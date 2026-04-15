#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
        sawFish2;
        deepvariant;
        trgt4_diseaseSTRs;
        trgt4_diseaseSTRs_plots;
        trgt4_diseaseSTRs_plots_meth;
        trgt4_all;
        trgt5_diseaseSTRs;
        trgt5_diseaseSTRs_plots;
        trgt5_diseaseSTRs_plots_meth;
        mosdepthROI;
        nanoStat;
        } from "../modules/dnaModules.nf" 


workflow PRE_PHASING {

    take:
    aligned

    main:

    dv_vcf_ch           = Channel.empty()
    dv_gvcf_ch          = Channel.empty()
    glnexus_manifest_ch = Channel.empty()
    sawfish_vcf_ch      = Channel.empty()
    str_vcf_ch          = Channel.empty()
    hiphase_input_ch    = Channel.empty()

    if (!params.skipVariants) {
        deepvariant(aligned)
        deepvariant.out.dv_vcf
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            | set {dv_vcf_ch}
    
        deepvariant.out.dv_gvcf
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            | set {dv_gvcf_ch}
        
        if (params.jointCall || params.jointSS) {
            deepvariant.out.dv_gvcf
            | map { meta, gvcf, tbi ->
                // store one record per sample: (caseID, meta, gvcfPath)
                tuple(meta.caseID, tuple(meta, gvcf.toString()))
            }
            .groupTuple()
            | map { caseID, records ->
                def anchorMeta = records[0][0]
                def content = records.collect { it[1] }.join('\n') + '\n'
                def mf = file("${caseID}.manifest")
                mf.text = content
                tuple(anchorMeta, mf)
            }
            | set { glnexus_manifest_ch }
        }
    }

    if (!params.skipSV) {
        sawFish2(aligned)
        sawFish2.out.sawfish_vcf //meta, vcf, idx
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            | set {sawfish_vcf_ch}
    }

    if (!params.skipSTR) {

        trgt4_all(aligned)
        trgt4_diseaseSTRs(aligned)

        trgt4_diseaseSTRs.out.str4_vcf
        | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
        | set {str_vcf_ch}

        trgt4_diseaseSTRs.out.trgt_full
            .combine(params.puretargetPlotGenes)
            |map {meta,bam,bai,vcf,tbi,genes -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi,strID:genes])}
            |set {trgt4_plot_ch}

        trgt4_diseaseSTRs_plots(trgt4_plot_ch)

        trgt4_diseaseSTRs.out.trgt_full
            |map {meta,bam,bai,vcf,tbi -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi])}
            |set {trgt4_plot_ch_meth}

        trgt4_diseaseSTRs_plots_meth(trgt4_plot_ch_meth)

        trgt5_diseaseSTRs(aligned)

        trgt5_diseaseSTRs.out.trgt_full
            .combine(params.puretargetPlotGenes)
            |map {meta,bam,bai,vcf,tbi,genes -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi,strID:genes])}
            |set {trgt5_plot_ch}

        trgt5_diseaseSTRs_plots(trgt5_plot_ch)

        trgt5_diseaseSTRs.out.trgt_full
            |map {meta,bam,bai,vcf,tbi -> 
            tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi])}
            |set {trgt5_plot_ch_meth}

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
        | set { hiphase_input_ch }
    }


    emit:
    dv_vcf                   = dv_vcf_ch
    dv_gvcf                  = dv_gvcf_ch
    glnexus_manifest         = glnexus_manifest_ch
    sawfish_vcf              = sawfish_vcf_ch
    sawfish_discover_dir     = params.skipSV ? Channel.empty() : sawFish2.out.sv_discover_dir
    sawfish_discover_dir2    = params.skipSV ? Channel.empty() : sawFish2.out.sv_discover_dir2
    sawfish_supporting_reads = params.skipSV ? Channel.empty() : sawFish2.out.sv_supporting_reads
    str4_vcf                 = params.skipSTR ? Channel.empty() : trgt4_diseaseSTRs.out.str4_vcf
    mosdepth                 = params.skipQC  ? Channel.empty() : mosdepthROI.out.multiqc
    nanoStat                 = params.skipQC  ? Channel.empty() : nanoStat.out.multiqc
    hiphaseInput             = hiphase_input_ch
}


