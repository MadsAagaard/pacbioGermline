#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
        sawFish2;
        deepvariant;
        trgt4_diseaseSTRs;
        trgt4_diseaseSTRs_plots;
        trgt4_diseaseSTRs_plots_meth;
        trgt5_all_adotto;
        trgt5_diseaseSTRs;
        trgt5_diseaseSTRs_plots;
        trgt5_diseaseSTRs_plots_meth;
        mosdepthROI;
        nanoStat;
        } from "../modules/dnaModules.nf" 


workflow PRE_PHASING {

    take:
    aligned     // val(meta),val(data) -> [hifiBam,hifiBai,failBam,failBai]

    main:

    dv_vcf_ch           = Channel.empty()
    dv_gvcf_ch          = Channel.empty()
    glnexus_manifest_ch = Channel.empty()
    sawfish_vcf_ch      = Channel.empty()
    str_vcf_ch          = Channel.empty()
    hiphase_input_ch    = Channel.empty()

    if (!params.skipVariants) {
        deepvariant(aligned)

        deepvariant.out.vcf
//            .map { meta, vcf, idx -> tuple(meta, [dvVcf: vcf, dvTbi: idx]) }
            .map { meta, vcf, idx -> tuple(meta, [dvVcf: vcf]) }

            .set { dv_vcf_ch }
    
    if (!params.skipSV) {
        sawFish2(aligned)
        
        sawFish2.out.sv_vcf
            //.map { meta, vcf, idx -> tuple(meta, [sawfishVcf: vcf, sawfishTbi: idx]) }
            .map { meta, vcf, idx -> tuple(meta, [sawfishVcf: vcf]) }            
            .set { sawfish_vcf_ch }

    }

    if (!params.skipSTR) {

        trgt4_diseaseSTRs(aligned)

        trgt4_diseaseSTRs.out.str4_vcf
            //.map { meta, vcf, idx -> tuple(meta, [str4Vcf: vcf, str4Tbi: idx]) }
            .map { meta, vcf, idx -> tuple(meta, [str4Vcf: vcf]) }
            .set { str_vcf_ch }

        trgt4_diseaseSTRs.out.trgt_full
            |map {meta,bam,bai,vcf,tbi -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi])}
            |set {trgt4_plot_ch}

        trgt4_diseaseSTRs_plots(trgt4_plot_ch)

        trgt4_diseaseSTRs_plots_meth(trgt4_plot_ch)

        trgt5_all_adotto(aligned)
        trgt5_all_adotto.out.adotto_vcf
            //.map { meta, vcf, idx -> tuple(meta, [str5AdottoVcf: vcf, str5AdottoIdx: idx]) }
            .map { meta, vcf, idx -> tuple(meta, [str5AdottoVcf: vcf]) }
            .set { str5_adotto_vcf_ch }

        trgt5_diseaseSTRs(aligned)

        trgt5_diseaseSTRs.out.trgt_full
            .map {meta,bam,bai,vcf,tbi -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi])}
            .set {trgt5_plot_ch}

        trgt5_diseaseSTRs_plots(trgt5_plot_ch)

        trgt5_diseaseSTRs.out.trgt_full
            .map {meta,bam,bai,vcf,tbi -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi])}
            .set {trgt5_plot_ch_meth}

        trgt5_diseaseSTRs_plots_meth(trgt5_plot_ch_meth)
        
    }

    if (!params.skipQC) {
        mosdepthROI(aligned)
        nanoStat(aligned)
    }

    // Assemble hiPhase input — all the pieces exist here already
    if (!params.skipVariants && !params.skipSV && !params.skipSTR) {
        aligned
        .join(dv_vcf_ch)
        .join(sawfish_vcf_ch)
        .join(str_vcf_ch)
        .join(str5_adotto_vcf_ch)
        .map{ meta, aln, dv, sv, str4, str5 -> 
            tuple(meta, aln + dv + sf + str4 +str5) //bamHifi,baiHifi, bamFail, baiFail, dv,sv,str4,str5
        }
        | set { hiphase_input_ch }
    }


    emit:
    dv_vcf                   = dv_vcf_ch
    dv_gvcf                  = deepvariant.out.dv_gvcf
    glnexus_manifest         = glnexus_manifest_ch
    sawfish_vcf              = sawfish_vcf_ch
    //sawfish_discover_dir     = params.skipSV ? Channel.empty() : sawFish2.out.sv_discover_dir
    sawfish_discover_dir    = params.skipSV ? Channel.empty() : sawFish2.out.sv_discover_dir2
    sawfish_supporting_reads = params.skipSV ? Channel.empty() : sawFish2.out.sv_supporting_reads
    str4_vcf                 = params.skipSTR ? Channel.empty() : trgt4_diseaseSTRs.out.str4_vcf
    mosdepth                 = params.skipQC  ? Channel.empty() : mosdepthROI.out.multiqc
    nanoStat                 = params.skipQC  ? Channel.empty() : nanoStat.out.multiqc
    hiphaseInput             = hiphase_input_ch
}


