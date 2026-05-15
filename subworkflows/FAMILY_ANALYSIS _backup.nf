#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include { glNexus_jointCall;
          sawFish2_jointCall_caseID;
          svdb_sawFish2_jointCall_caseID;
          exo14_2508_exome;
          exo14_2508_genome;
          exo14_2508_SV } from "../modules/dnaModules.nf"



workflow FAMILY_ANALYSIS {

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

    // Point outputDirTMP at the family base so outBase() resolves correctly
    params.outputDirTMP = params.familyDir

    // -------------------------------------------------------------------------
    // GLNexus joint calling
    //
    // Input:  gvcfManifest — plain text, one gVCF path per line
    // Output: joint-called VCF + WES ROI VCF
    //         → jointOutdir/jointCalls/
    // -------------------------------------------------------------------------
    Channel.of( tuple(anchorMeta, file(params.gvcfManifest)) )
    | set { glnexus_manifest_ch }

    glNexus_jointCall(glnexus_manifest_ch)

    // -------------------------------------------------------------------------
    // Sawfish joint calling
    //
    // Input:  sawfishCSV — one "discoverDir, bamPath" per line
    // Output: joint-called SV VCF
    //         → jointOutdir/jointCalls/
    // -------------------------------------------------------------------------
    Channel.of( tuple(anchorMeta, file(params.sawfishCSV)) )
    | set { sawfish_manifest_ch }

    sawFish2_jointCall_caseID(sawfish_manifest_ch)

    // -------------------------------------------------------------------------
    // SVDB annotation of Sawfish joint-call output
    //
    // Output: SVDB-annotated VCF + AF-filtered VCF (<10%)
    //         → jointOutdir/jointCalls/
    // -------------------------------------------------------------------------
    svdb_sawFish2_jointCall_caseID(sawFish2_jointCall_caseID.out.sv_jointCall_caseID_vcf)

    // -------------------------------------------------------------------------
    // Exomiser — only when --hpo is provided
    //
    // make_ped_and_family_v2.py reads --familySS to build the pedigree (.ped),
    // using the proband/pater/mater fields — identical to --jointSS + --hpo
    // in the normal main workflow.
    //
    // Three runs:
    //   exo14_2508_exome   → small variants, WES ROI VCF
    //   exo14_2508_genome  → small variants, whole genome VCF (Genomiser)
    //   exo14_2508_SV      → structural variants, SVDB-filtered Sawfish VCF
    // -------------------------------------------------------------------------
    if (params.hpo) {
        channel.fromPath(params.hpo)      | set { hpo_ch }
        channel.fromPath(params.familySS) | set { ss_ch  }

        // Small variant Exomiser (WES ROI)
        glNexus_jointCall.out.glnexus_wes_roi_vcf
            .combine(hpo_ch)
            .combine(ss_ch)
            | set { exomiser_ch }

        // Genomiser (whole genome)
        glNexus_jointCall.out.glnexus_vcf
            .combine(hpo_ch)
            .combine(ss_ch)
            | set { genomiser_ch }

        // SV Exomiser
        svdb_sawFish2_jointCall_caseID.out.sawfish_caseID_AF10
            .combine(hpo_ch)
            .combine(ss_ch)
            | set { exomiserSV_ch }

        exo14_2508_exome(exomiser_ch)
        exo14_2508_genome(genomiser_ch)
        exo14_2508_SV(exomiserSV_ch)
    }
}
