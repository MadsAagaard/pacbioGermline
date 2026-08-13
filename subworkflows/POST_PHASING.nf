#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * POST_PHASING — everything that runs on the phased HiFi alignment.
 */

include {
        kivvi_d4z4;
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
        svTopo;
        svTopo_filtered;
        mitorsaw;
        svdb_SawFish;
        } from "../modules/dnaModules.nf"


workflow POST_PHASING {

    take:
    phasedAll                 // tuple(meta, [bam, bai, failBam, failBai, dv_vcf, dv_idx, sawfish_vcf, sawfish_idx, sawfish_reads])
    sawfish_supporting_reads
    mosdepth
    nanoStat

    main:

    pbCPGtools(phasedAll)
    methBat(pbCPGtools.out)

    methBatNEW_pileup(phasedAll)
    methBatNEW_profile_single(methBatNEW_pileup.out.mC5)

    cramino(phasedAll)
    mitorsaw(phasedAll)
    whatsHap_stats(phasedAll)
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
