#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * FAMILY_ANALYSIS — joint calling + HPO prioritisation for a case/family.
 *
 * Changes:
 *   - params.outputDirTMP -> params.outputDirBase (single name for the output
 *     root; params.outBase() reads it).
 *   - Hard-coded "hg38v4LRS" / "HifiReads" in EXOMISER_ONLY_ENTRY replaced with
 *     params.genomeVersion / params.tagHifi, so a tag change in the config
 *     cannot silently desynchronise the retrospective path reconstruction.
 *   - Consistent dot notation throughout (mixing `|` with `.method()` causes
 *     Groovy precedence surprises).
 */

include { glNexus_jointCall;
          sawFish2_jointCall_caseID;
          svdb_sawFish2_jointCall_caseID;
          exo14_2508_exome;
          exo14_2508_genome;
          exo14_2508_SV } from "../modules/dnaModules.nf"


workflow FAMILY_ANALYSIS {

    take:
    glnexus_manifest_ch      // tuple(meta, manifestFile)
    sawfish_manifest_ch      // tuple(meta, manifestFile)
    hpo_ch                   // Channel.empty() when --hpo is not given
    ss_ch                    // samplesheet path channel

    main:
    glNexus_jointCall(glnexus_manifest_ch)
    sawFish2_jointCall_caseID(sawfish_manifest_ch)
    svdb_sawFish2_jointCall_caseID(sawFish2_jointCall_caseID.out.sv_jointCall_caseID_vcf)

    if (params.hpo) {

        glNexus_jointCall.out.glnexus_wes_roi_vcf
            .combine(hpo_ch)
            .combine(ss_ch)
            .set { exomiser_ch }

        glNexus_jointCall.out.glnexus_vcf
            .combine(hpo_ch)
            .combine(ss_ch)
            .set { genomiser_ch }

        svdb_sawFish2_jointCall_caseID.out.sawfish_caseID_AF10
            .combine(hpo_ch)
            .combine(ss_ch)
            .set { exomiserSV_ch }

        exo14_2508_exome(exomiser_ch)
        exo14_2508_genome(genomiser_ch)
        exo14_2508_SV(exomiserSV_ch)
    }
}


workflow FAMILY_ANALYSIS_ENTRY {

    // Load the family JSON written by pacbio_familyAnalysis_v3.sh (Step 5).
    def familyData = new groovy.json.JsonSlurper()
                         .parse(new File(params.familyJSON))

    // params.outBase(meta) with layoutMode=jointAnalysis resolves to
    //   ${params.outputDirBase}/jointAnalysis/${meta.caseID}_${params.readSet}
    // which is exactly what the shell script built in Step 3.
    def anchorMeta = [
        caseID     : familyData.caseID,
        id         : familyData.caseID,
        groupKey   : familyData.familyID,
        layoutMode : 'jointAnalysis',
        rekv       : '',
        testlist   : '',
    ]

    params.outputDirBase = params.familyDir
    params.layoutMode    = 'jointAnalysis'

    Channel.of( tuple(anchorMeta, file(params.gvcfManifest, checkIfExists: true)) )
        .set { glnexus_manifest_ch }

    Channel.of( tuple(anchorMeta, file(params.sawfishCSV, checkIfExists: true)) )
        .set { sawfish_manifest_ch }

    if (params.hpo) {
        channel.fromPath(params.hpo, checkIfExists: true).set { hpo_ch }
    }
    else {
        Channel.empty().set { hpo_ch }
    }

    channel.fromPath(params.familySS, checkIfExists: true).set { ss_ch }

    FAMILY_ANALYSIS(
        glnexus_manifest_ch,
        sawfish_manifest_ch,
        hpo_ch,
        ss_ch
    )
}


workflow EXOMISER_ONLY_ENTRY {

    /*
     * Re-runs the three Exomiser processes on a completed family analysis.
     * All joint-calling is skipped; VCFs are located by their publishDir paths
     * under jointCalls/.
     *
     * Required params (supplied by pacbio_familyAnalysis_v3.sh --exomiser-only):
     *   --familyJSON  --familyDir  --familySS  --hpo  [--genome]
     */

    def familyData = new groovy.json.JsonSlurper()
                         .parse(new File(params.familyJSON))

    def anchorMeta = [
        caseID     : familyData.caseID,
        id         : familyData.caseID,
        groupKey   : familyData.familyID,
        layoutMode : 'jointAnalysis',
        rekv       : '',
        testlist   : '',
    ]

    params.outputDirBase = params.familyDir
    params.layoutMode    = 'jointAnalysis'

    def jointCallsDir = "${familyData.jointOutdir}/jointCalls"
    def caseID        = familyData.caseID

    // Sourced from params so the reconstruction cannot drift from what the
    // modules actually wrote.
    def gv = params.genomeVersion   // 'hg38v4LRS'
    def rs = params.tagHifi         // 'HifiReads'

    Channel.of( tuple(
        anchorMeta,
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.WES_ROI.vcf.gz",     checkIfExists: true),
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.WES_ROI.vcf.gz.tbi", checkIfExists: true)
    ) ).set { wes_roi_vcf_ch }

    Channel.of( tuple(
        anchorMeta,
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.vcf.gz",     checkIfExists: true),
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.vcf.gz.tbi", checkIfExists: true)
    ) ).set { wgs_vcf_ch }

    Channel.of( tuple(
        anchorMeta,
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz",     checkIfExists: true),
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz.tbi", checkIfExists: true)
    ) ).set { sv_vcf_ch }

    channel.fromPath(params.hpo,      checkIfExists: true).set { hpo_ch }
    channel.fromPath(params.familySS, checkIfExists: true).set { ss_ch  }

    wes_roi_vcf_ch.combine(hpo_ch).combine(ss_ch).set { exomiser_ch }
    wgs_vcf_ch.combine(hpo_ch).combine(ss_ch).set    { genomiser_ch }
    sv_vcf_ch.combine(hpo_ch).combine(ss_ch).set     { exomiserSV_ch }

    exo14_2508_exome(exomiser_ch)
    exo14_2508_genome(genomiser_ch)
    exo14_2508_SV(exomiserSV_ch)
}
