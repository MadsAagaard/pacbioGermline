#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * POST_PHASING — everything that runs on the phased HiFi alignment.
 */

include {
        kivvi_d4z4;
        kivvi_d4z4_v11;
        kivvi_report;
        pbCPGtools;
        paraphase;
        paraphase4;
        starphase;
        methBat;
        methBatNEW_profile_single;
        methBatNEW_pileup;
        multiQC;
        cramino;
        nanoStat;
        whatsHap_stats;
        bcftools_stats;
        svTopo;
        svTopo_filtered;
        mitorsaw;
        svdb_SawFish;
        collect_germline_summary;
        germline_circos;
        } from "../modules/dnaModules.nf"


workflow POST_PHASING {

    take:
    phasedAll                 // tuple(meta, [bam, bai, failBam, failBai, dv_vcf, dv_idx, sawfish_vcf, sawfish_idx, sawfish_reads])
    sawfish_supporting_reads
    sawfish_copynum           // tuple(meta, copynum.bedgraph) — germline summary section 05 + germline_circos
    sawfish_depth_bw          // tuple(meta, gc_bias_corrected_depth.bw) — germline_circos
    mosdepth
    mosdepthSummary
    nanoStat

    main:

    pbCPGtools(phasedAll)
    methBat(pbCPGtools.out)

    methBatNEW_pileup(phasedAll)
    methBatNEW_profile_single(methBatNEW_pileup.out.mC5)

    cramino(phasedAll)
    mitorsaw(phasedAll)
    whatsHap_stats(phasedAll)
    bcftools_stats(phasedAll)
    paraphase(phasedAll)
    paraphase4(phasedAll)
    kivvi_d4z4(phasedAll)
    starphase(phasedAll)
    svTopo(phasedAll)
    svdb_SawFish(phasedAll)

    phasedAll
        .join(svdb_SawFish.out.sawfishAF10)
        .join(sawfish_supporting_reads)
        .map { meta, data, sv10_vcf, sv10_idx, sv_jsonReads ->
            tuple(meta, [
                bam           : data.bam,
                bai           : data.bai,
                sawfish10_vcf : sv10_vcf,
                sawfish10_idx : sv10_idx,
                sawfish_reads : sv_jsonReads
            ])
        }
        .set { phasedSawfishAF10 }

    svTopo_filtered(phasedSawfishAF10)

    //------ clinical summary (html channel) ------------------------------------------------------

    paraphase4.out.paraphase_json
    .join(methBatNEW_profile_single.out.icReport)
    .join(mosdepthSummary)
    .join(cramino.out.multiqc)
    .join(bcftools_stats.out.forSummary)
    .join(whatsHap_stats.out.multiqc)
    // Section 05 (PWS/AS): sawfish copy number + the phased sawfish SV VCF,
    // SVDB-annotated so every call carries its in-house frequency (OCC / FRQ).
    // No remainder: a sample only reaches POST_PHASING when sawFish2 has run
    // (phasedAll is built from its VCF in main.nf), and svdb_SawFish runs on
    // every sample of phasedAll, so both are always present.
    .join(sawfish_copynum)
    .join(svdb_SawFish.out.svdb_vcf.map { meta, vcf, tbi -> tuple(meta, vcf) })
    .map { meta, paraphase_json, icReport,mosDist,mosSummary,cramino,bcfStats,bcfFilters,whatshap,sawfishCopynum,sawfishVcf ->
        tuple(meta, [
            sawfishCopynum: sawfishCopynum,
            sawfishVcf: sawfishVcf,
            paraphase_json: paraphase_json,
            icReport: icReport,
            mosDist: mosDist,
            mosSummary: mosSummary,
            cramino: cramino,
            bcfStats: bcfStats,
            bcfFilters: bcfFilters,
            whatshap: whatshap
            ])
    }
    .set { clinical_summary_inputs_ch }

    collect_germline_summary(clinical_summary_inputs_ch)

    //------ germline circos (TEST PHASE) ----------------------------------------------------------
    // Same no-remainder reasoning as above: copynum, the depth bigWig and the SVDB VCF exist for
    // every sample of phasedAll. A failed task shows as a missing *.germlineCircos.html.
    if (!params.skipGermlineCircos) {
        phasedAll
        .join(sawfish_copynum)
        .join(sawfish_depth_bw)
        .join(svdb_SawFish.out.svdb_vcf.map { meta, vcf, tbi -> tuple(meta, vcf) })
        .map { meta, data, copynum, depthBw, svdbVcf ->
            tuple(meta, [
                depthBw : depthBw,
                copynum : copynum,
                svdbVcf : svdbVcf,
                dvVcf   : data.dv_vcf
            ])
        }
        .set { germline_circos_inputs_ch }

        germline_circos(germline_circos_inputs_ch)
    }


    //------ Kivvi 1.1 D4Z4 + FSHD report (TEST PHASE, next to Kivvi 1.0) -------------------------
    // Failures are ignored (conf/germline.config section 8): a missing *.kivviReport.html
    // in newToolsTest/customReports/ means Kivvi 1.1 or the report failed for that sample.
    if (!params.skipKivvi11) {
        kivvi_d4z4_v11(phasedAll)
        kivvi_report(kivvi_d4z4_v11.out.d4z4)
    }


    if (!params.skipQC) {
        Channel.empty()
            .mix(mosdepth)
            .mix(nanoStat)
            .mix(whatsHap_stats.out.multiqc)
            .map { meta, qcfile -> tuple(params.multiqcKey(meta), meta, qcfile) }
            .groupTuple(by: 0)
            .map { key, metas, qcfiles ->
                def meta0 = metas.find { it.relation == 'index' } ?: metas[0]
                tuple(meta0, qcfiles)
            }
            .set { multiqc_inputs_ch }

        multiQC(multiqc_inputs_ch)
    }
}
