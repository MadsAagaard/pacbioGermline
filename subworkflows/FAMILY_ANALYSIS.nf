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