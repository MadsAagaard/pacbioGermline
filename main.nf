#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import java.util.Locale

/*
 * =============================================================================
 * KG Vejle — PacBio HiFi LRS germline pipeline (hg38)
 * =============================================================================
 *
 * REFACTOR: separate HiFi / fail-read alignments
 * ----------------------------------------------
 * There is no merged all-reads alignment any more:
 *
 *   uBAM -> create_fofn -> pbmm2_align_hifi        -> HiFi BAM
 *                       -> pbmm2_align_failedOnly  -> fail BAM
 *
 *   TRGT consumes both (--reads / --fail-reads); everything else uses HiFi.
 *   hiPhaseTwoAln phases and emits BOTH alignments.
 *
 * Fixes vs. the previous version:
 *   - hiPhaseTwoAln was never imported, and the downstream joins referenced
 *     hiPhase.out.* (the superseded single-BAM process).
 *   - The --aligned re-entry branch was removed; the pipeline always starts
 *     from unmapped BAMs.
 *   - The summary/symlink processes referenced channels that only exist when
 *     --samplesheet is given, so --input alone crashed.
 *   - Destructuring assignments inside map closures had no `def`, so they wrote
 *     to the script-level binding from concurrently executing closures. `date2`
 *     in particular was clobbered by the filename parser.
 *   - `exit 0` on user input errors reported success to the caller.
 * =============================================================================
 */

def runDate      = new Date().format('yyMMdd')
def runDateTime  = new Date().format('yyMMdd HH:mm:ss')
def runUser      = System.getenv('USER') ?: 'unknown'
def runID        = "${runDate}.${runUser}"

def readModeText = params.failedReads ? 'fail reads only'
                 : params.hifiReads   ? 'HiFi only (TRGT runs without --fail-reads)'
                 : 'HiFi + fail reads, aligned separately'

log.info """\
======================================================
Clinical Genetics Vejle: PacBio LRS v4
======================================================
Genome        : ${params.genome}
GenomeDir     : ${params.refFilesDir}
Read mode     : ${readModeText}
Read set tag  : ${params.readSet}
STR name tag  : ${params.strTag}
RunID         : ${runID}
Script start  : ${runDateTime}
Genome FASTA  : ${params.genomeFasta}
Archive RAW   : ${params.dataArchive}
OutputDirBase : ${params.outputDirBase}
workDir       : ${workflow.workDir}
layout        : ${params.layoutMode}
error mode    : ${params.errorMode}
min HiFi GB   : ${params.minHifiGBeff}
"""


///////////////////////////////////////////////////
/////// ------- INPUT VALIDATION ------- //////////
///////////////////////////////////////////////////

if (!params.samplesheet && !params.input && !params.familySS) {
    exit 1, """
    USER INPUT ERROR: provide a samplesheet (--samplesheet) or an input folder
    containing the data to use (--input).
    """.stripIndent()
}

if (!params.samplesheet && params.hpo && !params.familySS) {
    exit 1, """
    USER INPUT ERROR: --hpo requires a samplesheet (--samplesheet) with the
    5 columns caseID, samplename, gender, relation, affection status.
    """.stripIndent()
}

if (params.aligned) {
    // --aligned was removed: this pipeline always starts from unmapped BAMs.
    // Analyses that begin from existing alignments are add-on work and belong
    // in a standalone script (cf. trgtGenomewide.nf), not in a re-entry branch
    // of main.nf that nothing exercises and nobody tests.
    // Kept as an explicit error rather than deleting the param, so the flag
    // fails loudly instead of being silently accepted and ignored.
    exit 1, """
    USER INPUT ERROR: --aligned is no longer supported.
    Run add-on analyses from existing alignments with a standalone script.
    """.stripIndent()
}

if (params.allReads) {
    log.warn """
    --allReads is deprecated. HiFi and fail reads are now aligned separately and
    TRGT consumes both via --fail-reads. The flag is ignored; the default
    behaviour is what you want. Use --hifiReads for HiFi-only.
    """.stripIndent()
}


///////////////////////////////////////////////////
/////// ------- SAMPLESHEET PARSING ------- ///////
///////////////////////////////////////////////////

def ssBaseName = params.ssBase

// -----------------------------------------------------------------------------
// Sex is validated once, here, where the raw gender field is actually read.
// A missing or malformed gender is a samplesheet data-entry error, not a state
// the pipeline should have a policy for: downstream code may assume
// meta.sex in ['male','female'].
//
// This replaces the old
//   (meta.sex=="male"||meta.sex=="M"||meta.genderFile=="M") ? "XY" : "XX"
// in Sawfish and every TRGT process, which silently genotyped anything it did
// not recognise as female.
// -----------------------------------------------------------------------------
def sexFromGender = { gender, sampleId ->
    switch ((gender ?: '').toString().trim().toLowerCase()) {
        case ['m', 'male']:
            return 'male'
        case ['k', 'f', 'female']:
            return 'female'
        default:
            throw new IllegalArgumentException(
                "Sample ${sampleId}: gender field is '${gender}', expected M or K. " +
                "Fix the samplesheet — the pipeline will not guess a karyotype.")
    }
}


if (params.samplesheet && !params.customSS && !params.jointSS) {

    channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(sep: '\t')
        .map { row ->
            def (rekv, npn, material, testlist, gender, proband, intRef) = row[0].tokenize('_')
            def groupKey = (intRef == 'noInfo') ? 'single' : intRef
            def sex      = sexFromGender(gender, npn)
            [ id       : npn,
              testlist : testlist,
              material : material,
              gender   : gender,
              sex      : sex,
              proband  : proband,
              intRef   : intRef,
              rekv     : rekv,
              groupKey : groupKey,
              ssBase   : ssBaseName ]
        }
        .set { samplesheet_full }

    samplesheet_full
        .branch { row ->
            singleSample : (row.groupKey == 'single')
            multiSample  : true
        }
        .set { samplesheetBranch }
}

if (params.samplesheet && (params.jointSS || params.familySS || params.customSS)) {

    channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(sep: '\t')
        .map { row ->
            def (rekv, npn, material, testlist, gender, proband, intRef) = row
            def sex = sexFromGender(gender, npn)

            def meta = [
                rekv     : rekv,
                id       : npn,
                material : material,
                testlist : testlist,
                gender   : gender,
                sex      : sex,
                proband  : proband,
                intRef   : intRef,
                ssBase   : ssBaseName,
                groupKey : intRef
            ]
            tuple(intRef, meta)
        }
        .groupTuple()
        .flatMap { intRef, metas ->

            def probands = metas.findAll { it.proband == 'T' }
            assert probands && probands.size() >= 1 : "No proband (T) found for intRef=${intRef}"

            def anchor = probands[0]
            def caseID = "${anchor.rekv}_${anchor.testlist}_${intRef}"

            metas.collect { m ->
                def relation = (m.proband == 'T') ? 'index'
                             : (m.gender  == 'M') ? 'pater'
                             : (m.gender  == 'K') ? 'mater'
                             : 'unknown_relation'

                m + [ caseID: caseID, relation: relation ]
            }
        }
        .set { samplesheet_full }
}


///////////////////////////////////////////////////
/////// ------- INPUT DISCOVERY ------- ///////////
///////////////////////////////////////////////////

// -------------------------------------------------------------------------
// uBAM discovery
// -------------------------------------------------------------------------
def inputBam
def root = params.input ?: params.dataArchive

if (params.hifiReads) {
    inputBam = params.input ? "${root}/**/*.hifi_reads.*.bam" : "${root}/**/hifi_reads/*.hifi_reads.*.bam"
}
else if (params.failedReads) {
    inputBam = params.input ? "${root}/**/*.fail_reads.*.bam" : "${root}/**/failed_reads/*.fail_reads.*.bam"
}
else {
    inputBam = params.input ? "${root}/**/*.bam" : "${root}/**/*_reads/*.bam"
}

channel.fromPath(inputBam, followLinks: true).set { rawBam }

if (params.samplesheet) {

    rawBam
        .map { bam ->
            def id = bam.baseName
            def (samplenameFull, pacbioID, readset, barcode) = id.tokenize('.')
            def (samplename, material, testlist, gender)     = samplenameFull.tokenize('_')
            tuple([id: samplename], bam)
        }
        .groupTuple(sort: true)
        .map { meta, bams ->
            long   totalBytes = (bams.sum { it.size() } as long)
            double totalGB    = totalBytes / (1024.0 * 1024 * 1024)
            tuple(meta + [nBams: bams.size(), totalsizeGB: totalGB], bams)
        }
        .branch { meta, bam ->
            UNASSIGNED : (meta.id =~ /UNASSIGNED/)
            samples    : true
        }
        .set { ubam_input }

    ubam_input.samples
        .map { meta, bam -> tuple(meta.id, meta, bam) }
        .set { ubam_input_samples }

    if (params.singleOnly) {
        samplesheetBranch.singleSample.map { row -> tuple(row.id, row) }.set { samplesheet_join }
    }
    else if (params.intrefOnly) {
        samplesheetBranch.multiSample.map { row -> tuple(row.id, row) }.set { samplesheet_join }
    }
    else {
        samplesheet_full.map { row -> tuple(row.id, row) }.set { samplesheet_join }
    }

    samplesheet_join
        .join(ubam_input_samples)
        .map { samplename, metaSS, metaData, bam -> tuple(metaSS + metaData, bam) }
        .set { ubam_ss_merged }

    // ---- summary lines -------------------------------------------------
    def summaryLine = { meta ->
        def gb = String.format(Locale.US, "%.2f", (meta.totalsizeGB as double))
        "${meta.id}\t${meta.nBams}\t${params.inputReadSet}\t${gb}\t${meta.testlist}"
    }
    def withHeader = { lines ->
        (["sample\tbamcount\treadSet\ttotal_gb\ttestlist"] + lines).join("\n")
    }

    ubam_ss_merged.map { meta, bams -> summaryLine(meta) }.collect()
        .map(withHeader).set { ubam_size_summary_ch }

    // ---- size gate -----------------------------------------------------
    ubam_ss_merged
        .branch { meta, bams ->
            keep : (meta.totalsizeGB as double) >= params.minHifiGBeff
            drop : true
        }
        .set { ubam_ss_merged_size_split }

    ubam_ss_merged_size_split.drop.map { meta, bams -> summaryLine(meta) }.collect()
        .map(withHeader).set { ubam_size_dropped_ch }

    ubam_ss_merged_size_split.keep.map { meta, bams -> summaryLine(meta) }.collect()
        .map(withHeader).set { ubam_size_keep_ch }

    // meta.sex was validated at parse time; nothing more to do here.
    ubam_ss_merged_size_split.keep.set { finalUbamInput }
}

if (!params.samplesheet) {

    rawBam
        .map { bam ->
            def id = bam.baseName
            def (samplenameFull, pacbioID, readset, barcode) = id.tokenize('.')
            def (instrument, runDateRaw, runTime)            = pacbioID.tokenize('_')
            def (samplename, material, testlist, gender)     = samplenameFull.tokenize('_')
            def meta = [
                id       : samplename,
                caseID   : "${runDateRaw}_${testlist}",
                gender   : gender,
                rundate  : runDateRaw,
                testlist : testlist,
                rekv     : 'noRekv',
                groupKey : 'single'
            ]
            tuple(meta, bam)
        }
        .groupTuple(sort: true)
        .map { meta, bams ->
            long   totalBytes = (bams.sum { it.size() } as long)
            double totalGB    = totalBytes / (1024.0 * 1024 * 1024)
            tuple(meta + [nBams: bams.size(), totalsizeGB: totalGB], bams)
        }
        .branch { meta, bam ->
            UNASSIGNED : (meta.id =~ /UNASSIGNED/)
            samples    : true
        }
        .set { ubam_input }

    ubam_input.samples.set { finalUbamInput }
}


///////////////////////////////////////////////////
/////// ------- MODULES / SUBWORKFLOWS ------- ////
///////////////////////////////////////////////////

include {
        write_input_summary;
        write_dropped_samples_summary;
        write_analyzed_samples_summary;
        symlinks_ubam_dropped;
        hiPhaseTwoAln;
        } from './modules/dnaModules.nf'

include { PREPROCESS }            from './subworkflows/PREPROCESS.nf'
include { PRE_PHASING }           from './subworkflows/PRE_PHASING.nf'
include { POST_PHASING }          from './subworkflows/POST_PHASING.nf'
include { FAMILY_ANALYSIS }       from './subworkflows/FAMILY_ANALYSIS.nf'
include { FAMILY_ANALYSIS_ENTRY } from './subworkflows/FAMILY_ANALYSIS.nf'


///////////////////////////////////////////////////
/////// ------- MAIN WORKFLOW ------- /////////////
///////////////////////////////////////////////////

workflow {

    if (params.samplesheet) {
        write_input_summary(ubam_size_summary_ch)
        write_analyzed_samples_summary(ubam_size_keep_ch)
        write_dropped_samples_summary(ubam_size_dropped_ch)
        symlinks_ubam_dropped(ubam_ss_merged_size_split.drop)
    }

    PREPROCESS(finalUbamInput)

    PRE_PHASING(PREPROCESS.out.alignedFinal)

    hiPhaseTwoAln(PRE_PHASING.out.hiphaseInput)

    hiPhaseTwoAln.out.hifi_bam
        .join(hiPhaseTwoAln.out.fail_bam, remainder: true)
        .map { meta, bam, bai, failBam, failBai ->
            tuple(meta, [bam: bam, bai: bai, failBam: failBam, failBai: failBai])
        }
        .join(hiPhaseTwoAln.out.dv_vcf)
        .join(hiPhaseTwoAln.out.sawfish_vcf)
        .join(PRE_PHASING.out.sawfish_supporting_reads)
        .map { meta, aln, dv_vcf, dv_idx, sv_vcf, sv_idx, sv_jsonReads ->
            tuple(meta, aln + [
                dv_vcf        : dv_vcf,
                dv_idx        : dv_idx,
                sawfish_vcf   : sv_vcf,
                sawfish_idx   : sv_idx,
                sawfish_reads : sv_jsonReads
            ])
        }
        .set { phasedAll }

    POST_PHASING(
        phasedAll,
        PRE_PHASING.out.sawfish_supporting_reads,
        PRE_PHASING.out.mosdepth,
        PRE_PHASING.out.nanoStat
    )

    if (params.hpo) {
        channel.fromPath(params.hpo, checkIfExists: true).set { hpo_ch }
    } else {
        Channel.empty().set { hpo_ch }
    }

    if (params.samplesheet) {
        channel.fromPath(params.samplesheet, checkIfExists: true).set { ss_ch }
    } else {
        Channel.empty().set { ss_ch }
    }

    if (params.jointCall || params.jointSS) {

        // Manifests are written to the run's runInfo dir rather than launchDir,
        // so they are recoverable as provenance after the run.
        def manifestDir = file("${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/manifests")
        manifestDir.mkdirs()

        PRE_PHASING.out.sawfish_discover_dir      // tuple(meta, discoverDir, hifiBamPath)
            .map { meta, dir, bam ->
                assert meta.caseID : "Joint calling requested but sample ${meta.id} has no caseID. Use --jointSS, or add caseID to the samplesheet parser."
                tuple(meta.caseID, tuple(meta, "${dir.toString()}, ${bam.toString()}"))
            }
            .groupTuple()
            .map { caseID, records ->
                def anchorMeta = records[0][0]
                def content    = records.collect { it[1] }.join('\n') + '\n'
                def mf         = manifestDir.resolve("${caseID}.sawFishJointCall.manifest.csv")
                mf.text        = content
                tuple(anchorMeta, mf)
            }
            .set { sawfish_jointCall_manifest_ch }

        PRE_PHASING.out.dv_gvcf
            .map { meta, gvcf, tbi ->
                assert meta.caseID : "Joint calling requested but sample ${meta.id} has no caseID."
                tuple(meta.caseID, tuple(meta, gvcf.toString()))
            }
            .groupTuple()
            .map { caseID, records ->
                def anchorMeta = records[0][0]
                def content    = records.collect { it[1] }.join('\n') + '\n'
                def mf         = manifestDir.resolve("${caseID}.glnexus.manifest")
                mf.text        = content
                tuple(anchorMeta, mf)
            }
            .set { glnexus_manifest_ch }

        FAMILY_ANALYSIS(
            glnexus_manifest_ch,
            sawfish_jointCall_manifest_ch,
            hpo_ch,
            ss_ch
        )
    }
}


///////////////////////////////////////////////////
/////// ------- COMPLETION ------- ////////////////
///////////////////////////////////////////////////

workflow.onComplete {

    // ---- failure manifest ---------------------------------------------------
    // params.errorMode = 'cohort' lets a bad sample be skipped so a large run
    // can finish. That is only acceptable if the skip is recorded somewhere a
    // human will see it.
    def failed = (workflow.stats.failedCount ?: 0) + (workflow.stats.ignoredCount ?: 0)

    if (failed > 0) {
        try {
            def f = file("${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/FAILED_TASKS.txt")
            f.parent.mkdirs()
            f.text = """\
                Run          : ${workflow.runName}
                Completed    : ${workflow.complete}
                Succeeded    : ${workflow.stats.succeedCount}
                Failed       : ${workflow.stats.failedCount}
                Ignored      : ${workflow.stats.ignoredCount}
                Error mode   : ${params.errorMode}

                Per-task detail is in the trace file under NextflowReports/.
                Output for the affected samples is INCOMPLETE and must not be
                reported without review.
                """.stripIndent()
            log.warn "*** ${failed} task(s) failed or were ignored — see ${f} ***"
        }
        catch (Exception e) {
            log.warn "Could not write FAILED_TASKS.txt: ${e.message}"
        }
    }

    // ---- symlink maintenance ------------------------------------------------
    if (!params.createSymlinks) {
        log.info "Symlink maintenance disabled by config."
        return
    }
    if (!workflow.success) {
        log.warn "Workflow failed — skipping symlink maintenance."
        return
    }

    def mirrorScript  = params.mirrorSampleData
    def collectScript = params.collectDataTypeSymlink

    if (!mirrorScript || !collectScript) {
        log.warn "Symlink script paths not defined in config — skipping."
        return
    }

    ["bash '${collectScript}'", "bash '${mirrorScript}'"].each { cmd ->
        log.info "onComplete: running: ${cmd}"
        try {
            def p = ["bash", "-lc", cmd].execute()
            p.waitForProcessOutput(System.out, System.err)
            if (p.exitValue() != 0) {
                log.warn "onComplete: command failed (exit ${p.exitValue()}): ${cmd}"
            } else {
                log.info "onComplete: finished OK: ${cmd}"
            }
        }
        catch (Exception e) {
            log.warn "onComplete: exception while running '${cmd}': ${e.message}"
        }
    }
}
