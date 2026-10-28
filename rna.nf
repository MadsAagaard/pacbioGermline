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
 *                   -> lrs_splice (blood recurrence store) -> lrs_splice_report
 *       all of the above -> germline_rna_summary (YAML / JSON / HTML per sample)
 *
 * Run:
 *   nextflow run MadsAagaard/pacbioGermline -r <tag> \
 *       -main-script rna.nf -profile slurm \
 *       --rnaSS pairs.tsv [--outdir <dir>]
 *
 * --rnaSS: tab-separated, one row per RNA sample; '#' lines and a header row
 * starting with "rnaNPN" are ignored.
 *     rnaNPN   dnaNPN   [phasedVcf]   [caseID]
 *   - dnaNPN NA (or -, empty)  -> use case 3, no DNA
 *   - phasedVcf (optional)     -> this file is used as is, no lookup (NA = look up)
 *   - caseID (optional)        -> family / case ID. Samples sharing it are one
 *                                 group in the blood recurrence store (relatives
 *                                 never count as recurrence for each other).
 *                                 Default group: dnaNPN, else rnaNPN.
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
      --rnaSS  [path]      rnaNPN <TAB> dnaNPN|NA [<TAB> phasedVcf|NA [<TAB> caseID]]

    OPTIONS
      --outdir [path]      output root (default /lnx01_data3/pacBioAnalyzedSamples/<genome>/<year>/RNAdata)
      --inputRNA [path]    folder to search for RNA uBAMs (**/*.lima.*.bam) instead of the archive
      --inhouseSplicingGenelist [path]
                           splicing panel (default: conf/germlineRNA.config;
                           null = all genes). Checked against the GTF at launch.
      --spliceGeneAliases [path]
                           legacy <TAB> current symbol, used when a panel symbol
                           is not in the GTF
      --rnaIdentityCheck false
                           skip the RNA vs DNA identity check (not recommended)
      --summaryFusions true
                           fusion cards in the summary (default false: counts only)
      --spliceRecurrenceDb [path|false]
                           blood recurrence store (default: conf/germlineRNA.config);
                           false = no store
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


include { RNA_ARM; RNA_PHASE_GERMLINE; RNA_SPLICING_GERMLINE; RNA_SUMMARY_GERMLINE } from './subworkflows/RNA.nf'


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
        problems << "line ${i + 1}: expected rnaNPN <TAB> dnaNPN|NA [<TAB> phasedVcf|NA [<TAB> caseID]], got '${l}'"
        return
    }
    def rnaNPN = f[0]
    def dnaNPN = isMissing(f[1]) ? null : f[1]
    def vcfCol = (f.size() > 2 && !isMissing(f[2])) ? f[2] : null
    def caseID = (f.size() > 3 && !isMissing(f[3])) ? f[3] : null

    if (rnaNPN.contains('.') || rnaNPN.contains('_') || rnaNPN.contains(' '))
        problems << "line ${i + 1}: rnaNPN '${rnaNPN}' contains '.', '_' or a space"
    if (vcfCol && !dnaNPN)
        problems << "line ${i + 1}: phasedVcf given but dnaNPN is NA — give the dnaNPN too, so the pair is named"
    if (caseID && !(caseID ==~ /[A-Za-z0-9._-]+/))
        problems << "line ${i + 1}: caseID '${caseID}' may only contain letters, digits, '.', '_' and '-'"
    if (f.size() > 4 && f[4..-1].any { it })
        problems << "line ${i + 1}: more than 4 columns (rnaNPN, dnaNPN, phasedVcf, caseID)"

    rows << [line: i + 1, rnaNPN: rnaNPN, dnaNPN: dnaNPN, vcfCol: vcfCol, caseID: caseID]
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
        // lrs_splice --patient = recurrence-store group key (leave-one-out): the
        // family / case when given, so relatives never count as recurrence
        patient    : r.caseID ?: r.dnaNPN ?: r.rnaNPN,
        useCase    : r.dnaNPN ? 'DNA+RNA' : 'RNA only',
        phasedVcf  : vcf?.toString(),
        vcfSource  : source,
        sampletype : 'rnaGermline',
        prefixRNA  : "${r.rnaNPN}.rnaGermline.${params.genomeVersion}",
        outName    : r.dnaNPN ? "${r.rnaNPN}_${r.dnaNPN}" : "${r.rnaNPN}_noDNA",
    ] + (r.caseID ? [caseID: r.caseID] : [:])   // key only when given: meta (and so every task hash) is unchanged for sheets without caseID
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


// =============================================================================
// SPLICING PANEL vs GTF  (launch check; same rule as lrs_splice load_panel)
// =============================================================================
// lrs_splice stops (exit 2) on any panel symbol it cannot find in the GTF, but
// only after preprocessing has run for hours. Resolve the panel here instead:
// a symbol is used as is when it is a gene_name in the GTF, else through the
// alias file (legacy <TAB> current). Anything left over stops the run now.

def panelResolved = []                  // [input, gtfSymbol, how]
def panelMissing  = []

if (params.inhouseSplicingGenelist) {
    def clean   = { String s -> s.replaceAll(/\r$/, '') }
    def wanted  = file(params.inhouseSplicingGenelist, checkIfExists: true).readLines()
        .collect { clean(it).split('\t')[0].trim() }
        .findAll { it && !it.startsWith('#') }
        .unique()
    def aliases = [:]
    if (params.spliceGeneAliases) {
        file(params.spliceGeneAliases, checkIfExists: true).readLines().each { l ->
            def f = clean(l).split('\t')
            if (f.size() >= 2 && f[0].trim() && f[1].trim() && !f[0].startsWith('#'))
                aliases[f[0].trim()] = f[1].trim()
        }
    }

    // gene_name of every 'gene' line in the GTF lrs_splice reads (plain or .gz)
    def gtf   = file(params.gencodeV49Gtf, checkIfExists: true)
    def known = new HashSet<String>()
    def input = gtf.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(gtf.newInputStream()) : gtf.newInputStream()
    input.withReader { rd ->
        String line
        while ((line = rd.readLine()) != null) {
            if (!line || line.charAt(0) == ('#' as char)) continue
            int t1 = line.indexOf('\t');      if (t1 < 0) continue
            int t2 = line.indexOf('\t', t1 + 1); if (t2 < 0) continue
            int t3 = line.indexOf('\t', t2 + 1); if (t3 < 0) continue
            if (line.substring(t2 + 1, t3) != 'gene') continue
            int i = line.indexOf('gene_name "', t3)
            if (i < 0) continue
            i += 11
            known << line.substring(i, line.indexOf('"', i))
        }
    }

    wanted.each { g ->
        if (g in known)                                 panelResolved << [g, g, 'GTF']
        else if (aliases[g] && aliases[g] in known)     panelResolved << [g, aliases[g], 'alias']
        else                                            panelMissing  << (aliases[g] ? "${g} (alias ${aliases[g]} also not in the GTF)" : g)
    }
    if (panelMissing) {
        problems << ("splicing panel ${params.inhouseSplicingGenelist}: ${panelMissing.size()} of ${wanted.size()} symbols not in ${gtf.name}" +
                     " — add legacy<TAB>current lines to ${params.spliceGeneAliases ?: 'an alias file (--spliceGeneAliases)'}:\n      " +
                     panelMissing.join('\n      '))
    }
}

// =============================================================================
// SPLICING RECURRENCE STORE  (blood; conf/germlineRNA.config)
// =============================================================================
if (params.spliceRecurrenceDb && !file(params.spliceRecurrenceDb.toString()).exists()) {
    problems << ("splicing recurrence store not found: ${params.spliceRecurrenceDb}\n" +
                 "      Build it with lrsCustomTools splicing/germline_recurrence_add.sh, or run without it: --spliceRecurrenceDb false")
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
    [ ['rnaNPN', 'dnaNPN', 'caseID', 'recurrenceKey', 'useCase', 'vcfSource', 'phasedVcf', 'outputFolder', 'ubams'].join('\t') ] +
    metas.collect { m ->
        [m.npnRNA, m.dnaNPN ?: 'NA', m.caseID ?: 'NA', m.patient, m.useCase, m.vcfSource, m.phasedVcf ?: 'NA',
         params.outBase(m), ubamsByNpn[m.npnRNA].join(',')].join('\t')
    }
).join('\n') + '\n'

if (panelResolved) {
    file("${runInfoDir}/splicingPanel_resolved.${params.traceStamp}.tsv").text =
        (['panelSymbol', 'gtfSymbol', 'resolvedBy'].join('\t') + '\n') +
        panelResolved.collect { it.join('\t') }.join('\n') + '\n'
}


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
Splicing panel  : ${params.inhouseSplicingGenelist ? "${params.inhouseSplicingGenelist} (${panelResolved.size()} genes, ${panelResolved.count { it[2] == 'alias' }} via alias)" : 'none (all genes)'}
Fusion gene list: ${params.inhouseFusionGenelist}
Summary fusions : ${params.summaryFusions ? 'cards' : 'counts only (--summaryFusions true for cards)'}
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

    RNA_SUMMARY_GERMLINE(
        RNA_ARM.out.isoseqForSummary,
        RNA_ARM.out.pigeonForSummary,
        RNA_ARM.out.fusionForSummary,
        RNA_ARM.out.fusionAnnot,
        RNA_PHASE_GERMLINE.out.identity,
        RNA_PHASE_GERMLINE.out.haplotagList,
        RNA_SPLICING_GERMLINE.out.json,
        RNA_SPLICING_GERMLINE.out.html
    )
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
