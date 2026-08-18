#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * PREPROCESS — uBAM -> two separate alignments.
 */

include {
        create_fofn;
        inputFiles_symlinks_ubam;
        pbmm2_align_hifi;
        pbmm2_align_failedOnly;
        } from "../modules/dnaModules.nf"


workflow PREPROCESS {

    take:
    finalUbamInput          // tuple(meta, [uBAMs])

    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    create_fofn(finalUbamInput)

    pbmm2_align_hifi(create_fofn.out.hifiReads)

    if (params.hifiReads) {

        // HiFi-only mode: no fail alignment, and TRGT runs without --fail-reads.
        pbmm2_align_hifi.out.bam
            .map { meta, hifiBam, hifiBai ->
                tuple(meta, [
                    hifiBam : hifiBam,
                    hifiBai : hifiBai,
                    failBam : null,
                    failBai : null
                ])
            }
            .set { alignedFinal_ch }
    }
    else {

        pbmm2_align_failedOnly(create_fofn.out.failReads)

        pbmm2_align_hifi.out.bam
            .join(pbmm2_align_failedOnly.out.bam)
            .map { meta, hifiBam, hifiBai, failBam, failBai ->
                tuple(meta, [
                    hifiBam : hifiBam,
                    hifiBai : hifiBai,
                    failBam : failBam,
                    failBai : failBai
                ])
            }
            .set { alignedFinal_ch }

    }

    emit:
    alignedFinal = alignedFinal_ch
}
