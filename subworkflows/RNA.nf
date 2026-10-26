#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * TEMPORARY COPY — germline RNA subworkflows (used by rna.nf only)
 * =============================================================================
 * RNA_PREPROCESS, RNA_TRANSCRIPTOME, RNA_FUSION and RNA_ARM are copied unchanged
 * from pacbioSomatic dev main.nf at c9ef209 (2026-10-23).
 * RNA_PHASE_GERMLINE and RNA_SPLICING_GERMLINE are germline versions:
 *   - the DNA VCF comes per sample in meta.phasedVcf (null = use case 3, no DNA)
 *   - RNA vs DNA identity check before haplotagging
 *   - no variant linking in lrs_splice (somaticVcf is always [])
 * Join key everywhere: meta.id (= rnaNPN).
 * =============================================================================
 */

include {
    merge_ubams;
    isoseq_refine_cluster;
    pbmm2_align_clusteredFLNC;
    pbmm2_align_refinedFLNC;
    isoseq_collapse;
    pigeon_classify;
    pbfusion;
    fusion_annotate;
    isocallProfile;
    rna_dna_identity;
    whatshap_haplotag;
    lrs_splice;
    lrs_splice_report;
} from '../modules/rnaModules.nf'


/* =============================================================================
 *  COPIED FROM pacbioSomatic (unchanged)
 * ========================================================================== */

workflow RNA_PREPROCESS {
    take: rna_ch
    main:
        rna_ch
            | branch { meta, bams ->
                multi:  (bams instanceof List) && (bams.size() > 1)
                single: true
            }
            | set { rna_branched }

        merge_ubams(rna_branched.multi)

        rna_branched.single
            | map { meta, bams -> tuple(meta, (bams instanceof List) ? bams[0] : bams) }
            | mix(merge_ubams.out)
            | set { ubam_final }                       // (meta, singleBam)

        isoseq_refine_cluster(ubam_final)

        pbmm2_align_refinedFLNC(isoseq_refine_cluster.out.isoseq_flnc_refined)
        pbmm2_align_clusteredFLNC(isoseq_refine_cluster.out.isoseq_flnc_clustered)

        isocallProfile(pbmm2_align_refinedFLNC.out.bam)

        // ---- base map: ------------------
        isoseq_refine_cluster.out.isoseq_flnc_refined
            .join(isoseq_refine_cluster.out.isoseq_flnc_clustered)
            .join(isoseq_refine_cluster.out.refine_report_json)
            .join(pbmm2_align_clusteredFLNC.out.bam)
            .join(pbmm2_align_refinedFLNC.out.bam)
            | map { meta, refinedFLNC, refinedFLNCPbi, clusteredFLNC, clusteredFLNCPbi, refine_json, clusteredBam, clusteredBai, refinedBAM, refinedBAI ->
                tuple(meta, [
                    refinedFLNC:      refinedFLNC,
                    refinedFLNCPbi:   refinedFLNCPbi,
                    clusteredFLNC:    clusteredFLNC,
                    clusteredFLNCPbi: clusteredFLNCPbi,
                    refineReportJSON: refine_json,
                    clusteredBAM:     clusteredBam,
                    clusteredBAI:     clusteredBai,
                    refinedBAM:       refinedBAM,
                    refinedBAI:       refinedBAI
                ])
            }
            | set { isoseq_pbmm2_joined }

        isoseq_collapse(isoseq_pbmm2_joined)

        // ---- Extend base map - still as metaMap for data ------------------
        isoseq_pbmm2_joined
            .join(isoseq_collapse.out.collapsed_list)
            .join(isoseq_collapse.out.read_stat)
            .map { meta, base, gff, counts, abund, readStat ->
                tuple(meta, base + [
                    collapsedGFF:  gff,
                    flncCounts:    counts,
                    flncAbundance: abund,
                    readStat:      readStat
                ])
            }
            .set { preprocess_all_joined }

    emit:
        preprocessFullOutput = preprocess_all_joined
        isoseqForSummary     = isoseq_refine_cluster.out.refine_report_json
        readStat             = isoseq_collapse.out.read_stat
        refinedAlignedBam    = pbmm2_align_refinedFLNC.out.bam
}

workflow RNA_TRANSCRIPTOME {
    take: preprocessFullOutput
    main:
        pigeon_classify(preprocessFullOutput)

    emit:
        classification      = pigeon_classify.out.classification_unfiltered
        classification_lite = pigeon_classify.out.classification
        sortedGFF           = pigeon_classify.out.sortedGFF
        pigeon              = pigeon_classify.out.pigeon
        pigeonForSummary    = pigeon_classify.out.pigeon_reports_json
}

workflow RNA_FUSION {
    take: preprocessFullOutput
    main:
        pbfusion(preprocessFullOutput)

        // fusion annotation -- TEST PHASE: annotates every pbfusion call
        // (panel membership is a tag, not a filter)
        fusionAnnot = Channel.empty()
        if (!params.skipFusionAnnotate) {
            fusion_annotate(
                pbfusion.out.bedpe,
                file(params.fusionAnnotation, checkIfExists: true),
                file(params.fusionDomains,    checkIfExists: true),
                file(params.fusionGeneList,   checkIfExists: true)
            )
            fusionAnnot = fusion_annotate.out.json
        }
    emit:
        fusion           = pbfusion.out.fusion
        fusionForSummary = pbfusion.out.inhouse_fusion
        fusionAnnot      = fusionAnnot                      // (meta, *.fusion_annot.json)
}

workflow RNA_ARM {
    take: rna_ch
    main:
        RNA_PREPROCESS(rna_ch)
        RNA_TRANSCRIPTOME(RNA_PREPROCESS.out.preprocessFullOutput)
        RNA_FUSION(RNA_PREPROCESS.out.preprocessFullOutput)
    emit:
        isoseqForSummary     = RNA_PREPROCESS.out.isoseqForSummary
        pigeonForSummary     = RNA_TRANSCRIPTOME.out.pigeonForSummary
        fusionForSummary     = RNA_FUSION.out.fusionForSummary
        preprocessFullOutput = RNA_PREPROCESS.out.preprocessFullOutput
        refinedAlignedBam    = RNA_PREPROCESS.out.refinedAlignedBam
        classification       = RNA_TRANSCRIPTOME.out.classification
        sortedGFF            = RNA_TRANSCRIPTOME.out.sortedGFF
        fusion               = RNA_FUSION.out.fusion
        fusionAnnot          = RNA_FUSION.out.fusionAnnot
}


/* =============================================================================
 *  GERMLINE VERSIONS
 * ========================================================================== */

/*
 * Use case 2 (meta.phasedVcf set): identity check -> whatshap haplotag.
 * Use case 3 (meta.phasedVcf null): the plain refined BAM goes straight on.
 * A sample whose identity check is not PASS is not haplotagged and not spliced;
 * it is listed in runInfo/<run>/IDENTITY_CHECK_NOT_PASSED.txt.
 * Each sample leaves on exactly one route.
 */
workflow RNA_PHASE_GERMLINE {
    take:
        rna_refined_bam                 // (meta, bam, bai)   refined pbmm2 RNA BAM
    main:
        rna_refined_bam
            | branch { meta, bam, bai ->
                withDna: meta.phasedVcf != null
                noDna:   true
            }
            | set { route }

        identity = Channel.empty()

        if (params.rnaIdentityCheck) {
            rna_dna_identity(route.withDna)
            identity = rna_dna_identity.out.result         // (meta, tsv, verdict)

            identity
                | branch { meta, tsv, verdict ->
                    pass: verdict == 'PASS'
                    fail: true
                }
                | set { idv }

            idv.fail
                | map { meta, tsv, verdict ->
                    log.warn "RNA ${meta.npnRNA} vs DNA ${meta.dnaNPN}: identity check ${verdict} — not haplotagged, not spliced. See ${params.outBase(meta)}/toolsOutputRNA/identityCheck/${tsv.name}"
                    "${meta.npnRNA}\t${meta.dnaNPN}\t${verdict}\t${meta.phasedVcf}\t${params.outBase(meta)}/toolsOutputRNA/identityCheck/${tsv.name}"
                }
                | collectFile(
                    name:     'IDENTITY_CHECK_NOT_PASSED.txt',
                    seed:     "rnaNPN\tdnaNPN\tverdict\tdnaVcf\tdetails",
                    newLine:  true,
                    sort:     true,
                    storeDir: "${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}")

            route.withDna
                | map  { meta, bam, bai -> [meta.id, meta, bam] }
                | join ( idv.pass.map { meta, tsv, verdict -> [meta.id] } )
                | map  { id, meta, bam -> tuple(meta, [bam: bam, phasedVcf: meta.phasedVcf]) }
                | set  { haplotag_in }
        }
        else {
            route.withDna
                | map { meta, bam, bai -> tuple(meta, [bam: bam, phasedVcf: meta.phasedVcf]) }
                | set { haplotag_in }
        }

        whatshap_haplotag(haplotag_in)

        whatshap_haplotag.out.bam
            | mix(route.noDna)
            | set { splice_bam }                            // (meta, bam, bai)

    emit:
        spliceBam    = splice_bam
        haplotagList = whatshap_haplotag.out.haplotag_list
        identity     = identity
}

workflow RNA_SPLICING_GERMLINE {
    take:
        rna_bam                         // (meta, bam, bai) haplotagged (use case 2) or plain (use case 3)
    main:
        rna_bam
            | map { meta, bam, bai ->
                tuple(meta, [
                    bam:        bam,
                    bai:        bai,
                    somaticVcf: []      // no variant linking for germline (common hets would override the recurrence cap)
                ])
            }
            | set { splice_in }

        // optional inputs: a real path, or [] (Nextflow stages nothing)
        def recDb     = params.spliceRecurrenceDb     ? file(params.spliceRecurrenceDb,     checkIfExists: true) : []
        def recExempt = params.spliceRecurrenceExempt ? file(params.spliceRecurrenceExempt, checkIfExists: true) : []
        lrs_splice(splice_in, recDb, recExempt)

        lrs_splice.out.json
            | map { meta, json -> tuple(meta, [json: json]) }
            | set { splice_report_in }

        lrs_splice_report(splice_report_in)

    emit:
        events = lrs_splice.out.events
        panel  = lrs_splice.out.panel
        json   = lrs_splice.out.json
        html   = lrs_splice_report.out.html
}
