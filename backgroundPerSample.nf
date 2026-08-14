#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * KG Vejle — PacBio LRS background POOL producer
 * =============================================================================
 *
 * uBAM -> pbmm2 (HiFi + fail, separately) -> TRGT v5 (adotto + TRExplorer v2)
 *
 * Produces per-sample genome-wide TR calls into a flat pool. It does NOT build
 * a cohort: which samples and how many is decided downstream in the LPS
 * toolkit, which is also where the yymmdd+letter stamp belongs.
 *
 *   <bgRoot>/repeatExpansions/background/adotto/         <- one dir, all samples
 *   <bgRoot>/repeatExpansions/background/adotto_LPS/
 *   <bgRoot>/repeatExpansions/background/trexplorer/
 *   <bgRoot>/repeatExpansions/background/trexplorer_LPS/
 *   <bgRoot>/runInfo/<batch>/                            <- inputs.tsv, batchInfo.tsv
 *
 * Run:
 *   nextflow run MadsAagaard/pacbioGermline -r refactor \
 *       -main-script backgroundPerSample.nf \
 *       -profile slurm,background \
 *       --bgSheet batch01.tsv
 *
 * -main-script takes ONE dash. Two dashes makes it a pipeline parameter, which
 * is ignored, and main.nf runs instead.
 *
 * ALWAYS FROM RAW uBAM. There is no aligned-input mode: the old merged-allReads
 * alignments predate the hifi+fail split and cannot be reused, and anything
 * produced by the refactored clinical pipeline already has adotto calls.
 *
 * BOTH CATALOGS RUN UNCONDITIONALLY: the alignments are not published, so every
 * catalog that might ever be wanted must be genotyped during this one visit to
 * the BAM. A third catalog added later means realigning the whole pool.
 *
 * READ MODE IS FIXED at HiFi+fail, so every file in the pool carries the same
 * strTag ('AllReadsNew'). --hifiReads / --failedReads are refused: see GUARDS.
 *
 * SPANNING BAMs: not produced. Both TRGT processes run --disable-bam-output, so
 * there is nothing to publish and no sort/index cost. Requires that the *_bam
 * output declarations have been removed from dnaModules.nf, or every task fails
 * on a missing output file.
 *
 * NO PHASING, on purpose. TRGT does not need it and the existing background
 * (adotto_n86.260812a) was built unphased. Changing that means rebuilding the
 * pool — phased and unphased calls are indistinguishable once mixed.
 * =============================================================================
 */

if (params.help) {
    log.info """
    Background pool producer — per-sample TRGT v5 (adotto + TRExplorer v2).

    REQUIRED
      --bgSheet    [path]  Concatenated LabWare metadata, one line per sample:
                             rekv_npn_material_testlist_gender_proband_intRef[_count]
                           Tab-separated also accepted. Only npn and gender are used.

    OPTIONS
      --bgRoot     [path]  Pool root. Defaults to the production pool in the
                           'background' profile — don't pass it unless you mean it.
      --initPool           Permit creating a pool that does not yet exist.
      --input      [path]  uBAM search root (default: params.dataArchive)
      --batch      [str]   Label for runInfo/ (default: timestamp)
      --minHifiGB  [int]   Min combined HiFi uBAM size per sample (default 30)
      --force              Re-run samples already in the pool
    """.stripIndent()
    exit 0
}

include { sexFromGender; parseMetaLine; parseUbamName } from './modules/lrsFunctions.nf'
include { PREPROCESS }                                 from './subworkflows/PREPROCESS.nf'
include { trgt5_all_adotto; trgt5_all_TRexplorer }     from './modules/dnaModules.nf'


// =============================================================================
// GUARDS
// =============================================================================

if (!params.bgRoot)  exit 1, "params.bgRoot is null — add: includeConfig 'conf/backgroundPerSample.config'  as the last line of the repo root nextflow.config."
if (!params.bgSheet) exit 1, "USER INPUT ERROR: --bgSheet is required."

// READ MODE IS FIXED AT HiFi+fail. params.strTag tracks the read mode
// ('AllReadsNew' by default, 'HifiReadsOnly' under --hifiReads) and is baked
// into every output filename — but the pool is one flat directory per catalog,
// so files from both modes would sit side by side and the LPS toolkit would
// glob them into a single background. Samples genotyped with fail reads and
// samples genotyped without are not comparable, and nothing in the VCFs says
// which is which.
//
// A background must also match the clinical calls, which always use HiFi+fail.
// If a HiFi-only pool is ever wanted, it needs its own --bgRoot, not a flag.
if (params.hifiReads || params.failedReads) {
    exit 1, "USER INPUT ERROR: --hifiReads and --failedReads are not supported here. The pool is one flat directory per catalog, so a second read mode would mix incomparable samples under strTag '${params.strTag}'. Use a separate --bgRoot if a HiFi-only pool is genuinely wanted."
}

// A typo'd --bgRoot does not error on its own: it silently starts a second pool
// and re-aligns everything into it. Require the pool to exist, or to be declared.
if (!file(params.bgPool).exists() && !params.initPool) {
    exit 1, "POOL NOT FOUND: ${params.bgPool} does not exist. Check --bgRoot for a typo, or pass --initPool if a new pool is intended."
}

// Without -profile background, the process publishDir blocks still point at the
// clinical tree and this run would overwrite clinical TRGT output under
// identical filenames. Probing outBase catches a missing profile even though the
// background run publishes nothing per-sample. Also catches profile ORDER:
// -profile background,slurm lets the slurm profile win on workDir.
def probeOut = params.outBase([id: '__probe__']).toString()
if (!probeOut.startsWith(params.bgRoot.toString()) || !params.lrsStorage.toString().startsWith(params.bgRoot.toString())) {
    exit 1, "OUTPUT SAFETY ABORT: publish targets resolve outside --bgRoot.\n  outBase    -> ${probeOut}\n  lrsStorage -> ${params.lrsStorage}\nThe 'background' profile is not active. Use: -profile slurm,background"
}

def runInfo = "${params.bgRoot}/runInfo/${params.batchTag}"

log.info """\
======================================================
PacBio LRS — background POOL producer
======================================================
Pool root   : ${params.bgRoot}
Pool dir    : ${params.bgPool}
Batch       : ${params.batchTag}
Sheet       : ${params.bgSheet}
uBAM root   : ${params.input ?: params.dataArchive}
Read mode   : HiFi + fail (fixed)
strTag      : ${params.strTag}
adotto      : ${params.trAdottoCatalog}
TRExplorer  : ${params.trExplorerCatalog}
workDir     : ${workflow.workDir}
"""

// Run-level provenance. Per-batch rather than per-sample: a batch is homogeneous
// by construction, and inputs.tsv maps samples to their batch. To trace a pool
// sample back: grep <npn> runInfo/*/inputs.tsv
file(runInfo).mkdirs()
file("${runInfo}/batchInfo.tsv").text = [
    "batch\t${params.batchTag}",
    "readMode\thifiPlusFail",
    "phased\tfalse",
    "genomeVersion\t${params.genomeVersion}",
    "genomeFasta\t${params.genomeFasta}",
    "strTag\t${params.strTag}",
    "adottoCatalog\t${params.trAdottoCatalog}",
    "trexplorerCatalog\t${params.trExplorerCatalog}",
    "trgtEnv\t${params.condaEnvs.trgt51}",
    "trgtLps\t${params.trgtLps}",
    "revision\t${workflow.revision ?: 'local'}",
    "commitId\t${workflow.commitId ?: 'NA'}",
    "started\t${new Date().format('yyMMdd_HHmmss')}",
].join('\n') + '\n'


// =============================================================================
// SAMPLE SHEET
// =============================================================================
// Ordinary LabWare metadata, concatenated. Fed raw rather than pre-extracted:
// `cut -d_ -f2,5` is one off-by-one from putting material in the gender column,
// and nothing downstream would notice. Batch prep is: cat metadata/*.txt > batch01.tsv
//
// Only npn and gender are carried. Not rekv/testlist/groupKey — those would let
// the clinical 'auto' output layout look satisfied if the overlay went missing.

channel.fromPath(params.bgSheet, checkIfExists: true)
    .splitText()
    .map { line -> line.trim() }
    .filter { line -> line && !line.startsWith('#') }
    .map { line ->
        def m = parseMetaLine(line)
        tuple(m.id, [id: m.id, gender: m.gender, sex: m.sex])
    }
    .set { sheet_ch }


// =============================================================================
// uBAM DISCOVERY
// =============================================================================

def ubamRoot = params.input ?: params.dataArchive
def ubamGlob = params.input ? "${ubamRoot}/**/*.bam" : "${ubamRoot}/**/*_reads/*.bam"

channel.fromPath(ubamGlob, followLinks: true)
    .map { bam -> tuple(parseUbamName(bam), bam) }
    .filter { parsed, bam -> parsed != null }
    .map { parsed, bam -> tuple(parsed.id, parsed.genderFile, bam) }
    .groupTuple()
    .join(sheet_ch, failOnDuplicate: true)   // inner join: archive samples not in the sheet drop out here
    .map { id, genderFiles, bams, meta ->

        // GENDER CROSS-CHECK. The Revio encodes gender in field 4 of the uBAM
        // name, same M/K vocabulary as LabWare — two independent records of the
        // same fact. A disagreement means the sample is mis-tracked, and a
        // mis-sexed sample poisons the pool invisibly: --karyotype changes the
        // chrX calls and the VCF does not say which karyotype produced them.
        //
        // ONLY npn AND gender. Not testlist/intRef/rekv — those describe a
        // referral, not a sample, and change legitimately for the same NPN.
        def seen = genderFiles.findAll { g -> g?.toString()?.trim() }.collect { g -> g.toString().trim() }.unique()

        if (!seen) {
            log.warn "${id}: no gender in uBAM filenames — cannot verify sheet gender '${meta.gender}'."
        }
        else if (seen.size() > 1) {
            throw new IllegalStateException("${id}: uBAMs disagree on gender (${seen.join(', ')}) — reads from more than one individual may be merged under this NPN.")
        }
        else {
            def fileSex = null
            try { fileSex = sexFromGender(seen[0], id) }
            catch (IllegalArgumentException e) { log.warn "${id}: uBAM gender '${seen[0]}' unrecognised — cannot verify sheet gender '${meta.gender}'." }

            if (fileSex && fileSex != meta.sex) {
                throw new IllegalStateException("${id}: GENDER MISMATCH. Sheet says '${meta.gender}' (${meta.sex}), uBAMs say '${seen[0]}' (${fileSex}). Either the sheet row is wrong or the data is assigned to the wrong NPN.")
            }
        }

        // Sorted explicitly, not via groupTuple(sort:true): with two collected
        // lists, relying on the operator to keep them aligned is not worth the
        // risk. Order matters for FOFN stability and so for -resume hashing.
        def sorted = bams.sort { a, b -> a.name <=> b.name }
        def gb     = (sorted.sum { b -> b.size() } as long) / (1024.0 * 1024 * 1024)

        tuple(meta + [nBams: sorted.size(), totalsizeGB: gb], sorted)
    }
    .set { ubamInput_ch }


// =============================================================================
// INPUT SUMMARY
// =============================================================================
// One row per sample found, with why it was or was not run. This is the record
// of what went into the pool from this batch — diff it against the sheet to see
// what had no sequencing data at all.

def inPool = { id ->
    !params.force &&
    file("${params.bgPool}/adotto").exists() &&
    (file("${params.bgPool}/adotto").list() as List).any { n -> n.startsWith("${id}.") }
}

ubamInput_ch
    .map { meta, bams ->
        def status = inPool(meta.id)                          ? 'already_in_pool'
                   : (meta.totalsizeGB < params.minHifiGBeff) ? "below_minHifiGB_${params.minHifiGBeff}"
                   :                                            'run'
        String.format(java.util.Locale.US, "%s\t%s\t%d\t%.2f\t%s",
                      meta.id, meta.sex, meta.nBams, meta.totalsizeGB, status)
    }
    .collectFile(name: 'inputs.tsv', storeDir: runInfo, newLine: true, sort: true,
                 seed: "npn\tsex\tn_ubam\ttotal_gb\tstatus")

ubamInput_ch
    .filter { meta, bams -> !inPool(meta.id) && meta.totalsizeGB >= params.minHifiGBeff }
    .set { toRun_ch }


// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    PREPROCESS(toRun_ch)

    trgt5_all_adotto(PREPROCESS.out.alignedFinal)
    trgt5_all_TRexplorer(PREPROCESS.out.alignedFinal)
}


workflow.onComplete {
    log.info """
    Batch ${params.batchTag}: ${workflow.success ? 'OK' : 'FAILED'} — ${workflow.stats.succeedCount} ok, ${workflow.stats.failedCount} failed, ${workflow.stats.ignoredCount} ignored
    Pool     : ${params.bgPool}/
    Run info : ${runInfo}/
    """.stripIndent()
}
