#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { 
        kivvi_d4z4;
        kivvi05_d4z4;
        pbCPGtools;
        paraphase;
        paraphase35
        starphase;
        methBat;
        multiQC;
        multiQC_ALL;
        mosdepthROI;
        cramino;
        nanoStat;
        whatsHap_stats;
        svTopo;
        svTopo_filtered;
        mitorsaw;
        svdb_SawFish;
        sawFish2_jointCall_all;
        svdb_sawFish2_jointCall_all;
        sawFish2_jointCall_caseID;
        svdb_sawFish2_jointCall_caseID;
        } from "../modules/dnaModules.nf" 

workflow POST_PHASING {

    take:
    phasedAll
    sawfish_supporting_reads
    mosdepth
    nanoStat

    main:
        pbCPGtools(phasedAll)
        methBat(pbCPGtools.out)
        cramino(phasedAll)
        mitorsaw(phasedAll)
        whatsHap_stats(phasedAll)
        paraphase(phasedAll)
        paraphase35(phasedAll)
        kivvi_d4z4(phasedAll)
        kivvi05_d4z4(phasedAll)
        starphase(phasedAll)
        svTopo(phasedAll)
        svdb_SawFish(phasedAll)

        /*
        hiPhase_OUT.hiphase_bam
        .join(svdb_SawFish.out.sawfishAF10)
        .join(sawfish_supporting_reads)
        | map {meta,bam,bai,sv10_vcf,sv10_idx,sv_jsonReads -> 
        tuple(meta,[bam:bam,bai:bai,sawfish10_vcf:sv10_vcf,sawfish10_idx:sv10_idx,sawfish_reads:sv_jsonReads])}
        |set {phasedSawfishAF10}   
        */
        phasedAll
        .join(svdb_SawFish.out.sawfishAF10)
        .join(sawfish_supporting_reads)
        | map { meta, data, sv10_vcf, sv10_idx, sv_jsonReads ->
            tuple(meta, [
                bam:            data.bam,
                bai:            data.bai,
                sawfish10_vcf:  sv10_vcf,
                sawfish10_idx:  sv10_idx,
                sawfish_reads:  sv_jsonReads
            ])
        }
        | set { phasedSawfishAF10 }


        svTopo_filtered(phasedSawfishAF10)

        if (!params.skipQC) {
            Channel.empty()
            .mix(mosdepth)
            .mix(nanoStat)
            .mix(whatsHap_stats.out.multiqc)
            .map { meta, qcfile ->
                tuple(params.multiqcKey(meta), meta, qcfile)
            }
            .groupTuple(by: 0)
            .map { key, metas, qcfiles ->

                // pick one representative meta for publishDir + naming
                def meta0 = metas.find { it.relation == 'index' } ?: metas[0]

                tuple(meta0, qcfiles)
            }
            .set { multiqc_inputs_ch }
            multiQC(multiqc_inputs_ch)
        }
}