#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * PREPROCESS — uBAM -> two separate alignments.
 *
 * Changes vs. the old version:
 *   - extractHifi and pbmm2_align (merged all-reads) are gone. HiFi and fail
 *     reads are aligned independently.
 *   - The map literal was missing a comma after hifiBai (syntax error).
 *   - --hifiReads now short-circuits the fail alignment instead of producing a
 *     failing task that errorStrategy would silently swallow.
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

        // NOTE: the join above is an inner join — a sample without a fail
        // alignment vanishes here. create_fofn hard-fails on an empty fail
        // FOFN so that condition is caught upstream with an attributable error
        // rather than a silent drop.
    }

    emit:
    alignedFinal = alignedFinal_ch
}
