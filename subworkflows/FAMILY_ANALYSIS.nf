#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


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
    hpo_ch                   // optional — Channel.empty() if not provided
    ss_ch                    // samplesheet path channel

    main:
    glNexus_jointCall(glnexus_manifest_ch)
    sawFish2_jointCall_caseID(sawfish_manifest_ch)
    svdb_sawFish2_jointCall_caseID(sawFish2_jointCall_caseID.out.sv_jointCall_caseID_vcf)

    if (params.hpo) {

        // Small variant Exomiser (WES ROI)
        glNexus_jointCall.out.glnexus_wes_roi_vcf
            .combine(hpo_ch)
            .combine(ss_ch)
            | exo14_2508_exome

        // Genomiser (whole genome)
        glNexus_jointCall.out.glnexus_vcf
            .combine(hpo_ch)
            .combine(ss_ch)
            | exo14_2508_genome

        // SV Exomiser
        svdb_sawFish2_jointCall_caseID.out.sawfish_caseID_AF10
            .combine(hpo_ch)
            .combine(ss_ch)
            | exo14_2508_SV
    }
}


workflow FAMILY_ANALYSIS_ENTRY {

    // -------------------------------------------------------------------------
    // Load family JSON written by pacbio.familyAnalysis.sh Step 5
    // -------------------------------------------------------------------------
    def familyData = new groovy.json.JsonSlurper()
                         .parse(new File(params.familyJSON))

    // -------------------------------------------------------------------------
    // Reconstruct anchorMeta
    //
    // params.outBase(meta) for layoutMode=jointAnalysis resolves to:
    //   "${params.outputDirTMP}/jointAnalysis/${meta.caseID}_${params.readSet}"
    //
    // We set params.outputDirTMP = params.familyDir below, so the full path
    // becomes:
    //   params.familyDir/jointAnalysis/<caseID>_AllAndHifi
    //
    // This matches exactly what the shell script built in Step 3.
    // -------------------------------------------------------------------------
    def anchorMeta = [
        caseID     : familyData.caseID,
        id         : familyData.caseID,   // used for process tags
        groupKey   : familyData.familyID,
        layoutMode : 'jointAnalysis',
        rekv       : '',
        testlist   : '',
    ]

    params.outputDirTMP = params.familyDir
    params.layoutMode   = 'jointAnalysis' 

    Channel.of( tuple(anchorMeta, file(params.gvcfManifest)) )
    | set { glnexus_manifest_ch }

    Channel.of( tuple(anchorMeta, file(params.sawfishCSV)) )
    | set { sawfish_manifest_ch }
    
    def hpo_ch = params.hpo ? channel.fromPath(params.hpo) : Channel.empty()
    def ss_ch  = channel.fromPath(params.familySS)

    FAMILY_ANALYSIS(
        glnexus_manifest_ch,
        sawfish_manifest_ch,
        hpo_ch,
        ss_ch
    )
}


workflow EXOMISER_ONLY_ENTRY {

    // -------------------------------------------------------------------------
    // Re-runs the three Exomiser processes on a completed family analysis.
    // All joint-calling steps are skipped — VCFs are located by their known
    // publishDir paths under jointCalls/.
    //
    // Required params (supplied by pacbio_familyAnalysis_v3.sh --exomiser-only):
    //   --familyJSON   path to <caseID>.family.json from the original run
    //   --familyDir    base family output dir
    //   --familySS     family samplesheet (for make_ped_and_family_v2.py)
    //   --hpo          path to HPO terms file
    //   --genome       genome build (default: hg38)
    //
    // VCF naming mirrors exactly what glNexus_jointCall and
    // svdb_sawFish2_jointCall_caseID publish:
    //   <caseID>.hg38v4LRS.HifiReads.deepVariant.jointCall.vcf.gz
    //   <caseID>.hg38v4LRS.HifiReads.deepVariant.jointCall.WES_ROI.vcf.gz
    //   <caseID>.hg38v4LRS.HifiReads.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    // -------------------------------------------------------------------------

    def familyData = new groovy.json.JsonSlurper()
                         .parse(new File(params.familyJSON))

    // anchorMeta — identical shape to FAMILY_ANALYSIS_ENTRY so publishDir
    // resolves to the same jointAnalysis/<caseID>_<readSet>/ folder.
    def anchorMeta = [
        caseID     : familyData.caseID,
        id         : familyData.caseID,
        groupKey   : familyData.familyID,
        layoutMode : 'jointAnalysis',
        rekv       : '',
        testlist   : '',
    ]

    params.outputDirTMP = params.familyDir
    params.layoutMode   = 'jointAnalysis'

    // -------------------------------------------------------------------------
    // Resolve the jointCalls/ directory from the family JSON's jointOutdir field.
    // This is the same path the shell script built and Nextflow published into.
    // -------------------------------------------------------------------------
    def jointCallsDir = "${familyData.jointOutdir}/jointCalls"

    def caseID        = familyData.caseID
    // genome_version and readSubset_hifiDefault are env vars set in the config;
    // their values are "hg38v4LRS" and "HifiReads" — hard-coded here to build
    // the file paths without needing the env block at workflow scope.
    def gv            = "hg38v4LRS"
    def rs            = "HifiReads"

    // -------------------------------------------------------------------------
    // Build input channels from existing on-disk files.
    // checkIfExists: true gives an early, clear error if a VCF is missing
    // (i.e. the original family analysis did not complete successfully).
    // -------------------------------------------------------------------------

    // Small-variant WES ROI VCF  →  exo14_2508_exome
    Channel.of( tuple(
        anchorMeta,
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.WES_ROI.vcf.gz",    checkIfExists: true),
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.WES_ROI.vcf.gz.tbi", checkIfExists: true)
    ) )
    | set { wes_roi_vcf_ch }

    // Whole-genome VCF  →  exo14_2508_genome (Genomiser)
    Channel.of( tuple(
        anchorMeta,
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.vcf.gz",     checkIfExists: true),
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.deepVariant.jointCall.vcf.gz.tbi", checkIfExists: true)
    ) )
    | set { wgs_vcf_ch }

    // SVDB-filtered Sawfish VCF  →  exo14_2508_SV
    Channel.of( tuple(
        anchorMeta,
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz",     checkIfExists: true),
        file("${jointCallsDir}/${caseID}.${gv}.${rs}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz.tbi", checkIfExists: true)
    ) )
    | set { sv_vcf_ch }

    // HPO and samplesheet channels — identical pattern to FAMILY_ANALYSIS_ENTRY
    channel.fromPath(params.hpo,      checkIfExists: true) | set { hpo_ch }
    channel.fromPath(params.familySS, checkIfExists: true) | set { ss_ch  }

    // -------------------------------------------------------------------------
    // Run the three Exomiser processes — same calls as in FAMILY_ANALYSIS,
    // but fed from file channels instead of process outputs.
    // -------------------------------------------------------------------------
    wes_roi_vcf_ch
        .combine(hpo_ch)
        .combine(ss_ch)
        | exo14_2508_exome

    wgs_vcf_ch
        .combine(hpo_ch)
        .combine(ss_ch)
        | exo14_2508_genome

    sv_vcf_ch
        .combine(hpo_ch)
        .combine(ss_ch)
        | exo14_2508_SV
}
