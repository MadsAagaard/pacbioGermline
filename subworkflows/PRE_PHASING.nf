#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * PRE_PHASING — everything that runs on the unphased alignments.
 */

include {
        sawFish2;
        deepvariant;
        trgt4_diseaseSTRs;
        trgt4_diseaseSTRs_plots;
        trgt4_diseaseSTRs_plots_meth;
        trgt5_all_adotto;
        trgt5_all_TRexplorer;
        trgt5_diseaseSTRs;
        trgt5_diseaseSTRs_plots;
        trgt5_diseaseSTRs_plots_meth;
        mosdepthROI;
        nanoStat;
        } from "../modules/dnaModules.nf"


workflow PRE_PHASING {

    take:
    aligned     // tuple(meta, [hifiBam, hifiBai, failBam, failBai])

    main:

    dv_vcf_ch            = Channel.empty()
    dv_gvcf_ch           = Channel.empty()
    sawfish_vcf_ch       = Channel.empty()
    sawfish_discover_ch  = Channel.empty()
    sawfish_reads_ch     = Channel.empty()
    str4_vcf_ch          = Channel.empty()
    str5_adotto_vcf_ch   = Channel.empty()
    mosdepth_ch          = Channel.empty()
    nanostat_ch          = Channel.empty()
    hiphase_input_ch     = Channel.empty()

    // -------------------------------------------------------------------------
    // Small variants
    // -------------------------------------------------------------------------
    if (!params.skipVariants) {

        deepvariant(aligned)

        deepvariant.out.dv_vcf
            .map { meta, vcf, idx -> tuple(meta, [dvVcf: vcf]) }
            .set { dv_vcf_ch }

        dv_gvcf_ch = deepvariant.out.dv_gvcf
    }

    // -------------------------------------------------------------------------
    // Structural variants  (sex-dependent: --expected-cn)
    // -------------------------------------------------------------------------
    if (!params.skipSV) {

        sawFish2(aligned)

        sawFish2.out.sv_vcf
            .map { meta, vcf, idx -> tuple(meta, [sawfishVcf: vcf]) }
            .set { sawfish_vcf_ch }

        sawfish_discover_ch = sawFish2.out.sv_discover_dir2
        sawfish_reads_ch    = sawFish2.out.sv_supporting_reads
    }

    // -------------------------------------------------------------------------
    // Repeat expansions  (sex-dependent: --karyotype)
    // -------------------------------------------------------------------------
    if (!params.skipSTR) {

        // --- TRGT v4, disease loci -------------------------------------------
        trgt4_diseaseSTRs(aligned)

        trgt4_diseaseSTRs.out.str4_vcf
            .map { meta, vcf, idx -> tuple(meta, [str4Vcf: vcf]) }
            .set { str4_vcf_ch }

        trgt4_diseaseSTRs.out.trgt_full
            .map { meta, bam, bai, vcf, tbi -> tuple(meta, [bam: bam, bai: bai, vcf: vcf, tbi: tbi]) }
            .set { trgt4_plot_ch }

        trgt4_diseaseSTRs_plots(trgt4_plot_ch)
        trgt4_diseaseSTRs_plots_meth(trgt4_plot_ch)

        // --- TRGT v5, genome-wide adotto catalog -----------------------------
        trgt5_all_adotto(aligned)

        trgt5_all_adotto.out.adotto_vcf
            .map { meta, vcf, idx -> tuple(meta, [str5AdottoVcf: vcf]) }
            .set { str5_adotto_vcf_ch }

        // --- TRGT v5, genome-wide TRexplorer catalog -----------------------------
        trgt5_all_TRexplorer(aligned)

        // --- TRGT v5, disease loci -------------------------------------------
        trgt5_diseaseSTRs(aligned)

        trgt5_diseaseSTRs.out.trgt_full
            .map { meta, bam, bai, vcf, tbi -> tuple(meta, [bam: bam, bai: bai, vcf: vcf, tbi: tbi]) }
            .set { trgt5_plot_ch }

        trgt5_diseaseSTRs_plots(trgt5_plot_ch)
        trgt5_diseaseSTRs_plots_meth(trgt5_plot_ch)
    }

    // -------------------------------------------------------------------------
    // QC
    // -------------------------------------------------------------------------
    if (!params.skipQC) {
        mosdepthROI(aligned)
        nanoStat(aligned)

        mosdepth_ch = mosdepthROI.out.multiqc
        nanostat_ch = nanoStat.out.multiqc
    }

    // -------------------------------------------------------------------------
    // hiPhase input: HiFi BAM + fail BAM + DeepVariant + Sawfish + both TRGT VCFs
    // -------------------------------------------------------------------------
    if (!params.skipVariants && !params.skipSV && !params.skipSTR) {

        aligned
            .join(dv_vcf_ch)
            .join(sawfish_vcf_ch)
            .join(str4_vcf_ch)
            .join(str5_adotto_vcf_ch)
            .map { meta, aln, dv, sv, str4, str5 ->
                // hifiBam, hifiBai, failBam, failBai, dvVcf, sawfishVcf,
                // str4Vcf, str5AdottoVcf
                tuple(meta, aln + dv + sv + str4 + str5)
            }
            .set { hiphase_input_ch }
    }
    else {
        log.warn "PRE_PHASING: phasing disabled — hiPhase requires DeepVariant, Sawfish and TRGT output."
    }

    emit:
    dv_gvcf                  = dv_gvcf_ch
    sawfish_discover_dir     = sawfish_discover_ch
    sawfish_supporting_reads = sawfish_reads_ch
    mosdepth                 = mosdepth_ch
    mosdepthSummary            = mosdepthROI.out.summary
    nanoStat                 = nanostat_ch
    hiphaseInput             = hiphase_input_ch
}
