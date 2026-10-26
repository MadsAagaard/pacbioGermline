#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * KG Vejle — PacBio LRS germline RNA (Kinnex, blood)          TEST PHASE
 * =============================================================================
 * Use case 2: blood DNA (already run with pacbioGermline) + blood RNA
 * Use case 3: blood RNA only
 *
 * RNA:  refine/cluster2 -> pbmm2 -> collapse -> pigeon      (transcriptome)
 *                                -> pbfusion -> fusion_annotate
 *       refined BAM -> [identity check -> whatshap haplotag]  (use case 2 only)
 *                   -> lrs_splice -> lrs_splice_report
 *
 * Run:
 *   nextflow run MadsAagaard/pacbioGermline -r <tag> \
 *       -main-script rna.nf -profile slurm \
 *       --rnaSS pairs.tsv [--outdir <dir>]
 *
 * --rnaSS: tab-separated, one row per RNA sample; '#' lines and a header row
 * starting with "rnaNPN" are ignored.
 *     rnaNPN   dnaNPN   [phasedVcf]
 *   - dnaNPN NA (or -, empty)  -> use case 3, no DNA
 *   - phasedVcf (optional)     -> this file is used as is, no lookup
 *
 * The germline DNA input is the hiPhase-phased DeepVariant VCF (+ .tbi beside it):
 *   <dnaNPN>.<genomeVersion>.<readSet tag>.hiphase.deepvariant.vcf.gz
 * Lookup at launch, in this order; each step must find exactly one file:
 *   1) phasedVcf column
 *   2) <germlineVcfStore>/<readSet>/deepVariant/vcfs/           (central store)
 *   3) the DNA case folders under <germlineCaseRoot>/<year>/perSampleAnalysis/
 *      (folder name contains _<dnaNPN>_, file in SNV_and_INDELs/)
 * Anything missing, ambiguous or without .tbi stops the run before any task.
 * =============================================================================
 */

if (params.help) {
    log.info """
    Germline RNA (TEST PHASE) — see the header of rna.nf.

    REQUIRED
      --rnaSS  [path]      rnaNPN <TAB> dnaNPN|NA [<TAB> phasedVcf]

    OPTIONS
      --outdir [path]      output root (default /lnx01_data3/pacBioAnalyzedSamples/<genome>/<year>/RNAdata)
      --inputRNA [path]    folder to search for RNA uBAMs (**/*.lima.*.bam) instead of the archive
      --inhouseSplicingGenelist [path]
                           splicing panel (default: none = all genes)
      --rnaIdentityCheck false
                           skip the RNA vs DNA identity check (not recommended)
    """.stripIndent()
    exit 0
}


// =============================================================================
// GUARDS
// =============================================================================

if (!params.rnaSS) {
    exit 1, "USER INPUT ERROR: --rnaSS is required (rnaNPN <TAB> dnaNPN|NA [<TAB> phasedVcf])."
}
if (params.samplesheet || params.familySS || params.jointSS || params.customSS) {
    exit 1, "USER INPUT ERROR: rna.nf takes --rnaSS only; --samplesheet/--familySS/--jointSS/--customSS belong to the DNA pipeline (main.nf)."
}
if (params.rnaIdentityCheck == null) {
    exit 1, "conf/germlineRNA.config was not loaded. nextflow.config includes it when --rnaSS is set; check that file."
}


include { RNA_ARM; RNA_PHASE_GERMLINE; RNA_SPLICING_GERMLINE } from './subworkflows/RNA.nf'


// =============================================================================
// SAMPLE SHEET  (read at launch, so every problem is reported before any task)
// =============================================================================

def isMissing = { v -> (v == null) || (v.toString().trim().toLowerCase() in ['', 'na', '-', '.', 'null']) }

def problems = []
def rows     = []

file(params.rnaSS, checkIfExists: true).readLines().eachWithIndex { line, i ->
    def l = line.trim()
    if (!l || l.startsWith('#')) return
    def f = l.split('\t', -1).collect { it.trim() }
    if (f[0].toLowerCase() == 'rnanpn') return                          // header row

    if (f.size() < 2) {
        problems << "line ${i + 1}: expected rnaNPN <TAB> dnaNPN|NA [<TAB> phasedVcf], got '${l}'"
        return
    }
    def rnaNPN = f[0]
    def dnaNPN = isMissing(f[1]) ? null : f[1]
    def vcfCol = (f.size() > 2 && !isMissing(f[2])) ? f[2] : null

    if (rnaNPN.contains('.') || rnaNPN.contains('_') || rnaNPN.contains(' '))
        problems << "line ${i + 1}: rnaNPN '${rnaNPN}' contains '.', '_' or a space"
    if (vcfCol && !dnaNPN)
        problems << "line ${i + 1}: phasedVcf given but dnaNPN is NA — give the dnaNPN too, so the pair is named"

    rows << [line: i + 1, rnaNPN: rnaNPN, dnaNPN: dnaNPN, vcfCol: vcfCol]
}

if (!rows && !problems) problems << "no sample rows in ${params.rnaSS}"

rows.groupBy { it.rnaNPN }.findAll { k, v -> v.size() > 1 }.each { k, v ->
    problems << "rnaNPN ${k} is listed ${v.size()} times (lines ${v*.line.join(', ')})"
}


// =============================================================================
// GERMLINE DNA VCF LOOKUP
// =============================================================================

def vcfName = { dnaNPN -> "${dnaNPN}.${params.genomeVersion}.*.hiphase.deepvariant.vcf.gz" }

// returns [path, source] or [null, problem]
def resolveVcf = { String dnaNPN, String explicit ->

    if (explicit) {
        def v = file(explicit)
        if (!v.exists())                     return [null, "phasedVcf not found: ${explicit}"]
        if (!file("${v}.tbi").exists())      return [null, "phasedVcf has no .tbi beside it: ${explicit}"]
        if (!v.name.startsWith("${dnaNPN}."))
            log.warn "phasedVcf for DNA ${dnaNPN} does not start with '${dnaNPN}.': ${v} — used as given (the identity check will test it)"
        if (!v.name.contains(".${params.genomeVersion}."))
            log.warn "phasedVcf for DNA ${dnaNPN} does not contain '.${params.genomeVersion}.' in its name: ${v} — check that it is the same reference build"
        return [v, 'samplesheet']
    }

    def steps = [
        [ 'centralStore', params.germlineReadSets.collect { rs -> "${params.germlineVcfStore}/${rs}/deepVariant/vcfs/${vcfName(dnaNPN)}" } ],
        [ 'caseFolder',   [ "${params.germlineCaseRoot}/*/perSampleAnalysis/*/*_${dnaNPN}_*/SNV_and_INDELs/${vcfName(dnaNPN)}" ] ],
    ]

    for (step in steps) {
        def (source, globs) = step
        def hits = globs.collectMany { g -> files(g) }.unique()
        if (hits.size() > 1)
            return [null, "DNA ${dnaNPN}: ${hits.size()} candidate VCFs in ${source} — put the right one in the phasedVcf column:\n      " + hits.join('\n      ')]
        if (hits.size() == 1) {
            def v = hits[0]
            if (!file("${v}.tbi").exists()) return [null, "DNA ${dnaNPN}: VCF has no .tbi beside it: ${v}"]
            return [v, source]
        }
    }
    return [null, "DNA ${dnaNPN}: no ${vcfName(dnaNPN)} found. Searched:\n      " + steps.collectMany { it[1] }.join('\n      ')]
}

def metas = rows.collect { r ->
    def vcf    = null
    def source = 'none'
    if (r.dnaNPN) {
        def (v, s) = resolveVcf(r.dnaNPN, r.vcfCol)
        if (v == null) problems << "line ${r.line}: ${s}"
        vcf    = v
        source = v ? s : 'NOT FOUND'
    }
    [
        id         : r.rnaNPN,                 // join key in all RNA subworkflows
        npn        : r.rnaNPN,                 // pbmm2 --sample
        npnRNA     : r.rnaNPN,
        dnaNPN     : r.dnaNPN,
        patient    : r.dnaNPN ?: r.rnaNPN,     // lrs_splice --patient (recurrence store key)
        useCase    : r.dnaNPN ? 'DNA+RNA' : 'RNA only',
        phasedVcf  : vcf?.toString(),
        vcfSource  : source,
        sampletype : 'rnaGermline',
        prefixRNA  : "${r.rnaNPN}.rnaGermline.${params.genomeVersion}",
        outName    : r.dnaNPN ? "${r.rnaNPN}_${r.dnaNPN}" : "${r.rnaNPN}_noDNA",
    ]
}

// =============================================================================
// RNA uBAMs  (one archive walk at launch; a sheet NPN without uBAMs is an input error)
// =============================================================================
// demultiplexed Kinnex uBAMs, e.g. 113125099265.m84328_..._s4.hifi_reads.lima.IsoSeqX_bc02_5p.SMfixed.bam

def rnaRoot   = params.inputRNA ?: params.dataArchiveRNA
def sheetNpns = metas*.npnRNA as Set
def ubamsByNpn = files("${rnaRoot}/**/*.lima.*.bam", followLinks: true)
    .findAll { bam -> bam.name.tokenize('.')[0] in sheetNpns }
    .groupBy { bam -> bam.name.tokenize('.')[0] }
    .collectEntries { npn, bams -> [npn, bams.sort { it.toString() }] }

metas.findAll { m -> !ubamsByNpn[m.npnRNA] }.each { m ->
    problems << "RNA ${m.npnRNA}: no ${m.npnRNA}.*.lima.*.bam under ${rnaRoot}"
}

if (problems) {
    exit 1, "\nGERMLINE RNA INPUT ERRORS (${params.rnaSS}) — nothing was started:\n  - " + problems.join('\n  - ') + "\n"
}


// =============================================================================
// RUN MANIFEST  (written at launch: which DNA file and uBAMs each RNA sample uses)
// =============================================================================

def runInfoDir = "${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}"
def manifest   = file("${runInfoDir}/germlineRNA_inputs.${params.traceStamp}.tsv")
manifest.parent.mkdirs()
manifest.text = (
    [ ['rnaNPN', 'dnaNPN', 'useCase', 'vcfSource', 'phasedVcf', 'outputFolder', 'ubams'].join('\t') ] +
    metas.collect { m ->
        [m.npnRNA, m.dnaNPN ?: 'NA', m.useCase, m.vcfSource, m.phasedVcf ?: 'NA',
         params.outBase(m), ubamsByNpn[m.npnRNA].join(',')].join('\t')
    }
).join('\n') + '\n'


log.info """\
======================================================
KG Vejle: PacBio LRS germline RNA (Kinnex)  TEST PHASE
======================================================
Genome          : ${params.genome} (${params.genomeVersion})
Genome FASTA    : ${params.genomeFasta}
RNA sample sheet: ${params.rnaSS}
Samples         : ${metas.size()}  (DNA+RNA: ${metas.count { it.dnaNPN }}, RNA only: ${metas.count { !it.dnaNPN }})
RNA uBAM search : ${params.inputRNA ?: params.dataArchiveRNA}/**/*.lima.*.bam
OutputDirBase   : ${params.outputDirBase}
Input manifest  : ${manifest}
Splicing panel  : ${params.inhouseSplicingGenelist ?: 'none (all genes)'}
Recurrence store: ${params.spliceRecurrenceDb ?: 'none'}
Identity check  : ${params.rnaIdentityCheck}
workDir         : ${workflow.workDir}
error mode      : ${params.errorMode}
lrsShared       : ${params.lrsSharedDir} (${params.lrsSharedVersion})
lrsCustomTools  : ${params.lrsCustomToolsBase}
"""


// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    Channel.fromList(metas)
        | map { meta -> tuple(meta, ubamsByNpn[meta.npnRNA]) }     // (meta, [uBAMs]); merged when > 1
        | set { rna_ch }

    RNA_ARM(rna_ch)

    RNA_PHASE_GERMLINE(RNA_ARM.out.refinedAlignedBam)

    RNA_SPLICING_GERMLINE(RNA_PHASE_GERMLINE.out.spliceBam)
}


// =============================================================================
// COMPLETION
// =============================================================================
workflow.onComplete {

    // errorMode 'cohort' (set in conf/germlineRNA.config) lets a failed task be
    // ignored so the other samples finish. Record it where a human will see it.
    def failed = (workflow.stats.failedCount ?: 0) + (workflow.stats.ignoredCount ?: 0)

    if (failed > 0) {
        try {
            def f = file("${runInfoDir}/FAILED_TASKS.germlineRNA.txt")
            f.parent.mkdirs()
            f.text = """\
                Run          : ${workflow.runName}
                Completed    : ${workflow.complete}
                Succeeded    : ${workflow.stats.succeedCount}
                Failed       : ${workflow.stats.failedCount}
                Ignored      : ${workflow.stats.ignoredCount}
                Error mode   : ${params.errorMode}
                Launch dir   : ${workflow.launchDir}

                Per-task detail is in the trace file under NextflowReports/.
                Output for the affected samples is INCOMPLETE.
                """.stripIndent()
            log.warn "*** ${failed} task(s) failed or were ignored — see ${f} ***"
        }
        catch (Exception e) {
            log.warn "Could not write FAILED_TASKS.germlineRNA.txt: ${e.message}"
        }
    }

    def idFail = file("${runInfoDir}/IDENTITY_CHECK_NOT_PASSED.txt")
    if (idFail.exists()) {
        log.warn "*** RNA/DNA identity check did not pass for some samples — see ${idFail} ***"
    }
}
