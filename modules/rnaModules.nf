#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * TEMPORARY COPY — germline RNA modules (pacbioGermline, used by rna.nf only)
 * =============================================================================
 * Copied from pacbioSomatic dev, modules/rnaModules.nf at c9ef209 (2026-10-23).
 * Replace with the shared RNA modules when the shared-modules setup lands.
 * Do not fix bugs here only: fix them in the somatic file too, or in the
 * shared version.
 *
 * Changes against the somatic file:
 *   - every publishDir goes through params.outBase(meta) (somatic: "${meta.id}/...")
 *   - TUMORBOARDFILES/RNA/ -> reportFiles/
 *   - merge_ubams: merged uBAM is not published (somatic copies it to the output)
 *   - pbfusion: thresholds from params (lrsShared v2 section 4b); the INHOUSE
 *     grep no longer fails the task when no call matches the gene list
 *   - lrs_splice: --patient ${meta.patient ?: meta.id}
 *   - NEW rna_dna_identity: RNA vs DNA genotype check before haplotagging
 *   - stub blocks, for -stub test runs
 *   - only the processes rna.nf uses (no ASE, OUTRIDER, ISA, oarfish,
 *     sqanti3, isoCall merge/call)
 *   - NEW germline_rna_summary (germline counterpart of collect_clinical_summaryRNA;
 *     not in the somatic file)
 * =============================================================================
 */


////////////////////////////////////////////
/////// ------- PREPROCESS + ALN ------- ///
////////////////////////////////////////////

process merge_ubams {
    label "high"
    tag "$meta.id"
    conda "${params.condaEnvs.isoseqPbtk}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.${meta.npnRNA}.merged.bam")

    script:
    if (bams instanceof List && bams.size() > 1)
        """
        samtools merge -f -@ ${task.cpus} ${meta.id}.${meta.npnRNA}.merged.bam ${bams.join(' ')}
        """
    else
        """
        cp ${bams} ${meta.id}.${meta.npnRNA}.merged.bam
        """

    stub:
    """
    touch ${meta.id}.${meta.npnRNA}.merged.bam
    """
}


process isoseq_refine_cluster {
    label "high"
    tag "$meta.id"
    conda "${params.condaEnvs.isoseqPbtk}"

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta),
        path("${meta.id}.${meta.npnRNA}.flnc.clustered.bam"),
        path("${meta.id}.${meta.npnRNA}.flnc.clustered.bam.pbi"),  emit: isoseq_flnc_clustered

    tuple val(meta),
        path("${meta.id}.${meta.npnRNA}.flnc.clustered.bam"),       emit: bam_for_fofn

    tuple val(meta),
        path("${meta.id}.${meta.npnRNA}.flnc.refined.bam"),
        path("${meta.id}.${meta.npnRNA}.flnc.refined.bam.pbi"),     emit: isoseq_flnc_refined

    tuple val(meta),
        path("${meta.id}.${meta.npnRNA}.flnc.refined.*.report.json"), emit: refine_report_json

    script:
    """
    isoseq refine \
    ${data[0]} \
    ${params.kinnexBarcodes12plex} \
    ${meta.id}.${meta.npnRNA}.flnc.refined.bam

    isoseq cluster2 \
    ${meta.id}.${meta.npnRNA}.flnc.refined.bam \
    ${meta.id}.${meta.npnRNA}.flnc.clustered.bam

    pbindex ${meta.id}.${meta.npnRNA}.flnc.clustered.bam
    pbindex ${meta.id}.${meta.npnRNA}.flnc.refined.bam
    """

    stub:
    """
    touch ${meta.id}.${meta.npnRNA}.flnc.clustered.bam ${meta.id}.${meta.npnRNA}.flnc.clustered.bam.pbi
    touch ${meta.id}.${meta.npnRNA}.flnc.refined.bam   ${meta.id}.${meta.npnRNA}.flnc.refined.bam.pbi
    touch ${meta.id}.${meta.npnRNA}.flnc.refined.filter_summary.report.json
    """
}

process pbmm2_align_clusteredFLNC {
    label "high"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    publishDir (
        path: { "${params.outBase(meta)}/alignments/" },
        mode: 'copy',
        pattern: '*.pbmm2.*')

    input:
    tuple val(meta), path(bam), path(pbi)

    output:
    tuple val(meta),
        path("${meta.prefixRNA}.clustered.pbmm2.bam"),
        path("${meta.prefixRNA}.clustered.pbmm2*bai"),  emit: bam

    script:
    """
    pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.npn} \
    ${params.genomeMmi} \
    ${bam} \
    ${meta.prefixRNA}.clustered.pbmm2.bam
    """

    stub:
    """
    touch ${meta.prefixRNA}.clustered.pbmm2.bam ${meta.prefixRNA}.clustered.pbmm2.bam.bai
    """
}

process pbmm2_align_refinedFLNC {
    label "high"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    publishDir (
        path: { "${params.outBase(meta)}/alignments/" },
        mode: 'copy',
        pattern: '*.pbmm2.*')

    input:
    tuple val(meta), path(bam), path(pbi)

    output:
    tuple val(meta),
        path("${meta.prefixRNA}.refined.pbmm2.bam"),
        path("${meta.prefixRNA}.refined.pbmm2*bai"),  emit: bam

    script:
    """
    pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.npn} \
    ${params.genomeMmi} \
    ${bam} \
    ${meta.prefixRNA}.refined.pbmm2.bam
    """

    stub:
    """
    touch ${meta.prefixRNA}.refined.pbmm2.bam ${meta.prefixRNA}.refined.pbmm2.bam.bai
    """
}


process isoseq_collapse {
    label "high"
    tag "$meta.id"
    conda "${params.condaEnvs.isoseqPbtk}"

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/isoseq/collapsed/" },
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.collapsed.*"), emit: all

    tuple val(meta),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.gff"),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.flnc_count.txt"),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.abundance.txt"), emit: collapsed_list

    tuple val(meta),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.read_stat.txt"), emit: read_stat

    script:
    """
    isoseq collapse \
    --do-not-collapse-extra-5exons \
    ${data.refinedBAM} \
    ${data.refinedFLNC} \
    ${meta.prefixRNA}.refinedFLNC.collapsed.gff
    """

    stub:
    """
    for x in gff flnc_count.txt abundance.txt read_stat.txt; do
        touch ${meta.prefixRNA}.refinedFLNC.collapsed.\$x
    done
    """
}

process pigeon_classify {
    label "high"
    tag "$meta.id"
    conda "${params.condaEnvs.pigeon}"

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/isoseq/pigeon/" },
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
        path("*.pigeon*"),                           emit: pigeon

    tuple val(meta),
        path("*.filtered_lite_classification.txt"),  emit: classification

    tuple val(meta),
        path("*.pigeon_classification.txt"),         emit: classification_unfiltered

    tuple val(meta),
        path("*.collapsed.sorted.gff"),              emit: sortedGFF

    tuple val(meta),
        path("*.filtered.report.json"),
        path("*.pigeon.report.json"),                emit: pigeon_reports_json

    script:
    """
    cp ${data.collapsedGFF} ${meta.prefixRNA}.pbmm2.collapsed.gff
    pigeon prepare ${meta.prefixRNA}.pbmm2.collapsed.gff

    pigeon classify \
    ${meta.prefixRNA}.pbmm2.collapsed.sorted.gff \
    ${params.gencodeV49SortedGtf} \
    ${params.genomeFasta} \
    --fl ${data.flncCounts} \
    --cage-peak ${params.pigeonRefTSS} \
    --poly-a ${params.pigeonPolyA} \
    -o ${meta.prefixRNA}.pigeon

    pigeon filter \
    ${meta.prefixRNA}.pigeon_classification.txt \
    ${meta.prefixRNA}.pbmm2.collapsed.sorted.gff

    pigeon report \
    --exclude-singletons \
    ${meta.prefixRNA}.pigeon_classification.filtered_lite_classification.txt \
    ${meta.prefixRNA}.pigeon.report.txt
    """

    stub:
    """
    touch ${meta.prefixRNA}.pbmm2.collapsed.sorted.gff
    touch ${meta.prefixRNA}.pigeon_classification.txt
    touch ${meta.prefixRNA}.pigeon_classification.filtered_lite_classification.txt
    touch ${meta.prefixRNA}.pigeon_classification.filtered.report.json
    touch ${meta.prefixRNA}.pigeon.report.json
    """
}


////////////////////////////////////////////
/////// ------- FUSION ------- /////////////
////////////////////////////////////////////

process pbfusion {
    label "medium"
    tag "$meta.id"
    conda "${params.condaEnvs.pbfusion}"

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/pbfusion/" },
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.{pdf,bed,txt,vcf,idx}"),  emit: fusion

    tuple val(meta), path("${meta.prefixRNA}.refinedBAM.PBfusion.INHOUSE.txt"),  emit: inhouse_fusion

    tuple val(meta), path("${meta.prefixRNA}.refinedBAM.fusion.breakpoints.groups.bed"), emit: bedpe

    script:
    // data.refinedBAM = pbmm2-aligned refined FLNC BAM (per-molecule reads, NOT cluster2 output)
    """
    pbfusion discover \
    -b ${data.refinedBAM} \
    --threads ${task.cpus} \
    --gtf ${params.gencodeV49Gtf} \
    --min-coverage ${params.pbfusionMinCoverage} \
    --min-fusion-fraction ${params.pbfusionMinFusionFraction} \
    --min-fusion-quality ${params.pbfusionMinFusionQuality} \
    -o ${meta.prefixRNA}.refinedBAM.fusion

    # gene-list view; exit 1 (no match) is a valid empty result, exit 2 (error) still fails
    grep -w -f ${params.inhouseFusionGenelist} \
        ${meta.prefixRNA}.refinedBAM.fusion.breakpoints.groups.bed \
        > ${meta.prefixRNA}.refinedBAM.PBfusion.INHOUSE.txt || [ \$? -eq 1 ]
    """

    stub:
    """
    touch ${meta.prefixRNA}.refinedBAM.fusion.breakpoints.groups.bed
    touch ${meta.prefixRNA}.refinedBAM.PBfusion.INHOUSE.txt
    """
}

process fusion_annotate {
    label 'low'
    tag "$meta.id"
    conda "${params.condaEnvs.somaticSummaryEnv}"            // any python3; the tool is stdlib-only

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/fusion_annotate/" },
        mode: 'copy')

    input:
    tuple val(meta), path(bedpe)    // pbfusion *.breakpoints.groups.bed
    path annotation                 // MANE Select genePred asset
    path domains                    // reconciled domain table (*.mane.tsv)
    path geneList                   // panel, current HGNC symbols -- tags in_panel, filters nothing

    output:
    tuple val(meta), path("${meta.prefixRNA}.fusion_annot.json"),         emit: json
    tuple val(meta), path("${meta.prefixRNA}.fusion_annot.txt"),          emit: report
    tuple val(meta), path("${meta.prefixRNA}.fusion_annot*.svg"),         emit: svg,     optional: true
    tuple val(meta), path("${meta.prefixRNA}.fusion_annot.chimera.fa"),   emit: protein, optional: true
    tuple val(meta), path("${meta.prefixRNA}.fusion_annot.versions.yml"), emit: versions

    script:
    def out = "${meta.prefixRNA}.fusion_annot"
    """
    set -euo pipefail

    python3 ${params.fusionAnnotatePy} \\
        --bedpe          ${bedpe} \\
        --annotation     ${annotation} \\
        --domains        ${domains} \\
        --gene-list      ${geneList} \\
        --fasta          ${params.genomeFasta} \\
        --min-confidence ${params.fusionMinConfidence} \\
        --json           ${out}.json \\
        --report         ${out}.txt \\
        --svg            ${out}.svg \\
        --protein-fasta  ${out}.chimera.fa

    # contract check: every pbfusion call is annotated or listed as skipped
    python3 -c "import json,sys; d=json.load(open('${out}.json')); sys.exit(0 if d['schema_version'] >= 2 and d['n_calls_in'] == len(d['fusions']) + len(d['skipped']) else 'fusion_annotate: JSON contract check failed')"

    cat <<-END_VERSIONS > ${out}.versions.yml
    "${task.process}":
        fusion_annotate: \$(python3 ${params.fusionAnnotatePy} --version | cut -d' ' -f2)
        annotation: "${annotation.name} \$(sha256sum ${annotation} | cut -c1-64)"
        domains: "${domains.name} \$(sha256sum ${domains} | cut -c1-64)"
        gene_list: "${geneList.name} \$(sha256sum ${geneList} | cut -c1-64)"
        fasta: "${params.genomeFasta}"
    END_VERSIONS
    """

    stub:
    def out = "${meta.prefixRNA}.fusion_annot"
    """
    touch ${out}.json ${out}.txt ${out}.versions.yml
    """
}


////////////////////////////////////////////
/////// ------- ISOCALL ------- ////////////
////////////////////////////////////////////

process isocallProfile {
    label "medium"
    tag "$meta.id"

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/isoCall/profile/" },
        mode: 'copy')
    publishDir (
        path: { "${params.lrsStorage}/RNA/isoCallProfiles/" },
        mode: 'copy')

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.prefixRNA}.isoCallProfile.gz"),  emit: profile

    script:
    """
    ${params.isocallBin} profile \
    --reads ${bam} \
    --output ${meta.prefixRNA}.isoCallProfile.gz
    """

    stub:
    """
    touch ${meta.prefixRNA}.isoCallProfile.gz
    """
}


////////////////////////////////////////////
/////// ------- PHASING ------- ////////////
////////////////////////////////////////////

/*
 * RNA vs DNA identity check (NEW). The DNA VCF path is set by hand or by lookup,
 * so a wrong pairing would silently give wrong haplotags. At DNA hom-alt SNPs
 * the RNA of the same person carries (almost) only the alt allele; another
 * person is often het or hom-ref there. Allele-specific expression cannot fake
 * a pass or a fail at hom-alt sites, which is why het sites are not used.
 * The task always succeeds and reports VERDICT (PASS | FAIL | INSUFFICIENT);
 * the subworkflow haplotags only PASS samples.
 */
process rna_dna_identity {
    label 'low'
    tag "$meta.id"
    conda "${params.condaEnvs.aseGatk}"                     // bcftools + htslib

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/identityCheck/" },
        mode: 'copy',
        pattern: "*.identityCheck.tsv")

    input:
    tuple val(meta), path(bam), path(bai)   // refined pbmm2 RNA BAM; DNA VCF = meta.phasedVcf (.tbi beside it)

    output:
    tuple val(meta), path("${meta.prefixRNA}.identityCheck.tsv"), env(VERDICT), emit: result

    script:
    def out = "${meta.prefixRNA}.identityCheck.tsv"
    """
    set -euo pipefail

    # DNA hom-alt biallelic PASS SNPs, chrM excluded (heteroplasmy)
    bcftools view -f PASS -v snps -m2 -M2 -i 'GT="AA"' -t ^chrM \\
        ${meta.phasedVcf} -Oz -o dna.homAlt.vcf.gz

    # RNA allele depths at those sites (-T streams the BAM; no BAQ; high max depth
    # for highly expressed genes; ref-skip (N) does not count as coverage)
    bcftools mpileup -f ${params.genomeFasta} -T dna.homAlt.vcf.gz \\
        -a FORMAT/AD -B -d 100000 --skip-indels --ignore-RG \\
        -q ${params.rnaIdentityMinMapQ} -Q ${params.rnaIdentityMinBaseQ} \\
        -Ou ${bam} \\
      | bcftools query -f '%CHROM\\t%POS\\t[%AD]\\n' > rna.AD.tsv

    awk -F'\\t' -v OFS='\\t' \\
        -v minDP=${params.rnaIdentityMinDepth} -v maxRef=${params.rnaIdentityMaxRefFrac} \\
        -v minSites=${params.rnaIdentityMinSites} -v maxDisc=${params.rnaIdentityMaxDiscordant} \\
        -v rna=${meta.npnRNA} -v dna=${meta.dnaNPN} -v vcf=${meta.phasedVcf} '
        { k = split(\$3, a, ","); t = 0; for (i = 1; i <= k; i++) t += a[i]
          if (t < minDP) next
          n++; if (a[1] / t > maxRef) d++ }
        END {
          n += 0; d += 0
          f = (n > 0) ? d / n : 1
          v = (n < minSites) ? "INSUFFICIENT" : ((f > maxDisc) ? "FAIL" : "PASS")
          print "rnaNPN", "dnaNPN", "sitesTested", "sitesDiscordant", "discordantFrac", "verdict", "minDepth", "maxRefFrac", "maxDiscordantFrac", "dnaVcf"
          printf "%s\\t%s\\t%d\\t%d\\t%.4f\\t%s\\t%d\\t%s\\t%s\\t%s\\n", rna, dna, n, d, f, v, minDP, maxRef, maxDisc, vcf
        }' rna.AD.tsv > ${out}

    VERDICT=\$(tail -n 1 ${out} | cut -f 6)
    """

    stub:
    def out = "${meta.prefixRNA}.identityCheck.tsv"
    // -stub test hook: an rnaNPN starting with "MISMATCH" gets FAIL
    def v = meta.npnRNA.toString().startsWith('MISMATCH') ? 'FAIL' : 'PASS'
    """
    printf 'rnaNPN\\tdnaNPN\\tsitesTested\\tsitesDiscordant\\tdiscordantFrac\\tverdict\\n${meta.npnRNA}\\t${meta.dnaNPN}\\t0\\t0\\t0\\t${v}\\n' > ${out}
    VERDICT=${v}
    """
}

process whatshap_haplotag {
    label 'medium'
    tag "$meta.id"
    conda "${params.condaEnvs.whatshap}"                    // NB: env must also provide samtools

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/haplotag/" },
        mode: 'copy',
        pattern: "*.haplotag.tsv")
    publishDir (
        path: { "${params.outBase(meta)}/alignments/" },
        mode: 'copy',
        pattern: "*.haplotagged.ba*")

    input:
    tuple val(meta), val(data)   // data: bam, phasedVcf   (phasedVcf .tbi sits beside it)

    output:
    tuple val(meta), path("${meta.prefixRNA}.refined.pbmm2.haplotagged.bam"),
                     path("${meta.prefixRNA}.refined.pbmm2.haplotagged.bam.bai"), emit: bam
    tuple val(meta), path("${meta.prefixRNA}.refined.pbmm2.haplotag.tsv"),        emit: haplotag_list

    script:
    """
    whatshap haplotag \
        --reference ${params.genomeFasta} \
        --ignore-read-groups \
        --output-haplotag-list ${meta.prefixRNA}.refined.pbmm2.haplotag.tsv \
        --output ${meta.prefixRNA}.refined.pbmm2.haplotagged.bam \
        ${data.phasedVcf} \
        ${data.bam}

    samtools index ${meta.prefixRNA}.refined.pbmm2.haplotagged.bam
    """

    stub:
    """
    touch ${meta.prefixRNA}.refined.pbmm2.haplotagged.bam ${meta.prefixRNA}.refined.pbmm2.haplotagged.bam.bai
    touch ${meta.prefixRNA}.refined.pbmm2.haplotag.tsv
    """
}


////////////////////////////////////////////
/////// ------- SPLICING ------- ///////////
////////////////////////////////////////////

process lrs_splice {
    label 'medium'
    tag "$meta.id"
    conda "${params.condaEnvs.somaticSummaryEnv}"                // python3 + pysam (nothing else)

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/splicing/aberrantSplicing/" },
        mode: 'copy',
        pattern: "*.aberrantSplicing.*")
    publishDir (
        path: { "${params.outBase(meta)}/reportFiles/" },
        mode: 'copy',
        pattern: "*.aberrantSplicing.panel.tsv")
    // collection point for a future blood recurrence store (germline RNA: inside the RNAdata tree)
    publishDir (
        path: { "${params.lrsStorage}/splicing/json/" },
        mode: 'copy',
        pattern: "*.aberrantSplicing.json")

    input:
    tuple val(meta), val(data)
    path recurrenceDb       // cohort recurrence store (tsv.gz), or [] when not used
    path recurrenceExempt   // known pathogenic junctions (tsv), or [] when not used

    output:
    tuple val(meta), path("${meta.prefixRNA}.aberrantSplicing.tsv"),        emit: events
    tuple val(meta), path("${meta.prefixRNA}.aberrantSplicing.json"),       emit: json
    tuple val(meta), path("${meta.prefixRNA}.aberrantSplicing.panel.tsv"),  emit: panel, optional: true

    script:
    def somaticVcf = (data.somaticVcf == null || data.somaticVcf instanceof List) ? '' : "--somatic-vcf ${data.somaticVcf}"
    def genelist2  = params.inhouseSplicingGenelist ? "--genes ${params.inhouseSplicingGenelist}" : ''
    def aliases    = params.spliceGeneAliases ? "--gene-aliases ${params.spliceGeneAliases}" : ''
    // patient key: must equal the patient column used when adding samples to the store
    def patient    = "--patient ${meta.patient ?: meta.id}"
    def tumourType = meta.pcgr ? "--tumour-type ${meta.pcgr}" : ''
    def recurrence = recurrenceDb     ? "--recurrence-db ${recurrenceDb} --recurrence-mode ${params.spliceRecurrenceMode}" : ''
    def exempt     = recurrenceExempt ? "--recurrence-exempt ${recurrenceExempt}" : ''
    """
    python3 ${params.lrsSplicePy} \
        --bam                 ${data.bam} \
        --gtf                 ${params.gencodeV49Gtf} \
        --fasta               ${params.genomeFasta} \
        ${genelist2} \
        ${aliases} \
        ${somaticVcf} \
        --sample              ${meta.prefixRNA} \
        ${patient} \
        ${tumourType} \
        ${recurrence} \
        ${exempt} \
        --high-min-reads      ${params.spliceHighMinReads} \
        --tandem-max-shift    ${params.spliceTandemMaxShift} \
        --annotated-cap       ${params.spliceAnnotatedCap} \
        --min-reads           ${params.spliceMinReads} \
        --min-usage           ${params.spliceMinUsage} \
        --min-anchor          ${params.spliceMinAnchor} \
        --min-cluster-depth   ${params.spliceMinDepth} \
        --max-fdr             ${params.spliceMaxFDR} \
        --ir-min-frac         ${params.spliceIrMinFrac} \
        --variant-window      ${params.spliceVariantWindow} \
        --out-tsv             ${meta.prefixRNA}.aberrantSplicing.tsv \
        --out-json            ${meta.prefixRNA}.aberrantSplicing.json

    # report view: PASS events only (column looked up by header name, not index)
    awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="filter") c=i; print; next} c && \$c=="PASS"' \\
        ${meta.prefixRNA}.aberrantSplicing.tsv > ${meta.prefixRNA}.aberrantSplicing.panel.tsv

    # drop a header-only panel file so the optional output stays truly optional
    if [ \$(wc -l < ${meta.prefixRNA}.aberrantSplicing.panel.tsv) -le 1 ]; then
        rm -f ${meta.prefixRNA}.aberrantSplicing.panel.tsv
    fi
    """

    stub:
    """
    touch ${meta.prefixRNA}.aberrantSplicing.tsv ${meta.prefixRNA}.aberrantSplicing.json
    echo "bam=${data.bam} patient=${meta.patient ?: meta.id} genes=${params.inhouseSplicingGenelist ?: 'ALL'}" > ${meta.prefixRNA}.aberrantSplicing.stubinfo.txt
    """
}

process lrs_splice_report {
    label 'low'
    tag "$meta.id"
    conda "${params.condaEnvs.somaticSummaryEnv}"                // stdlib only, but keep one env

    publishDir (
        path: { "${params.outBase(meta)}/toolsOutputRNA/splicing/aberrantSplicing/" },
        mode: 'copy',
        pattern: "*.aberrantSplicing*.html")
    publishDir (
        path: { "${params.outBase(meta)}/reportFiles/" },
        mode: 'copy',
        pattern: "*.aberrantSplicing.report.html")

    input:
    tuple val(meta), val(data)                    // data: json

    output:
    tuple val(meta), path("${meta.prefixRNA}.aberrantSplicing.section.html"), emit: section
    tuple val(meta), path("${meta.prefixRNA}.aberrantSplicing.report.html"),  emit: html

    script:
    """
    # embeddable fragment (for a later RNA summary)
    python3 ${params.lrsSpliceReportPy} \
        --json  ${data.json} \
        --out   ${meta.prefixRNA}.aberrantSplicing.section.html \
        --mode  fragment

    # standalone, for review / QC
    python3 ${params.lrsSpliceReportPy} \
        --json  ${data.json} \
        --out   ${meta.prefixRNA}.aberrantSplicing.report.html \
        --mode  standalone
    """

    stub:
    """
    touch ${meta.prefixRNA}.aberrantSplicing.section.html ${meta.prefixRNA}.aberrantSplicing.report.html
    """
}


////////////////////////////////////////////
/////// ------- SUMMARY ------- ////////////
////////////////////////////////////////////

/*
 * Germline RNA summary: YAML + JSON + self-contained HTML per RNA sample
 * (lrsCustomTools summaryScripts/germlineRNA). TEST PHASE banner in the HTML.
 * The HTML sits in reportFiles/ next to the lrs_splice report it links to.
 */
process germline_rna_summary {
    label 'low'
    tag "$meta.id"
    conda "${params.condaEnvs.somaticSummaryEnv}"                // python3 + pyyaml

    publishDir (
        path: { "${params.outBase(meta)}/reportFiles/" },
        mode: 'copy',
        pattern: "*.summaryRNA.html")
    publishDir (
        path: { "${params.outBase(meta)}/summaryFiles/" },
        mode: 'copy',
        pattern: "*.summaryRNA.*")
    publishDir (
        path: { "${params.lrsStorage}/summaries/rna/" },
        mode: 'copy',
        pattern: "*.summaryRNA.{json,yaml}")

    input:
    tuple val(meta), val(data)
    // data: refineReportJSON, pigeonRawJSON, pigeonFilteredJSON, fusionInhouse
    //       fusionAnnot    fusion_annotate JSON, or []
    //       identityTsv    rna_dna_identity TSV, or []        (use case 2 only)
    //       haplotagList   whatshap haplotag list, or []      (identity PASS only)
    //       spliceJSON     lrs_splice JSON, or []
    //       spliceReport   file name of the standalone lrs_splice report, or null
    //       spliceNote     why spliceJSON is [], or null

    output:
    tuple val(meta), path("${meta.prefixRNA}.summaryRNA.yaml"), emit: yaml
    tuple val(meta), path("${meta.prefixRNA}.summaryRNA.json"), emit: json
    tuple val(meta), path("${meta.prefixRNA}.summaryRNA.html"), emit: html

    script:
    def opt   = { flag, v -> (v == null || v instanceof List) ? '' : "${flag} ${v}" }
    def panel = params.inhouseSplicingGenelist ? "--panel-name ${file(params.inhouseSplicingGenelist).name}" : ''
    """
    python3 ${params.germlineRnaSummaryPy} \
        --rna-npn                ${meta.npnRNA} \
        --dna-npn                ${meta.dnaNPN ?: 'NA'} \
        --case-id                ${meta.caseID ?: 'NA'} \
        --use-case               '${meta.useCase}' \
        --vcf-source             ${meta.vcfSource ?: 'NA'} \
        --phased-vcf             ${meta.phasedVcf ?: 'NA'} \
        --genome-version         ${params.genomeVersion} \
        ${panel} \
        --refine-json            ${data.refineReportJSON} \
        --pigeon-raw             ${data.pigeonRawJSON} \
        --pigeon-filtered        ${data.pigeonFilteredJSON} \
        --pbfusion               ${data.fusionInhouse} \
        ${opt('--fusion-annot', data.fusionAnnot)} \
        --identity-check-enabled ${params.rnaIdentityCheck ? 'true' : 'false'} \
        --identity-min-sites     ${params.rnaIdentityMinSites} \
        ${opt('--identity-tsv', data.identityTsv)} \
        ${opt('--haplotag-list', data.haplotagList)} \
        ${opt('--splice-json', data.spliceJSON)} \
        ${data.spliceReport ? "--splice-report-href ${data.spliceReport}" : ''} \
        ${data.spliceNote ? "--splice-note '${data.spliceNote}'" : ''} \
        --html-template          ${params.germlineRnaSummaryHtml} \
        --output                 ${meta.prefixRNA}.summaryRNA
    """

    stub:
    """
    touch ${meta.prefixRNA}.summaryRNA.yaml ${meta.prefixRNA}.summaryRNA.json ${meta.prefixRNA}.summaryRNA.html
    echo "splice=${data.spliceJSON instanceof List ? 'none' : data.spliceJSON} note=${data.spliceNote} identity=${data.identityVerdict} hap=${data.haplotagList instanceof List ? 'none' : 'yes'} annot=${data.fusionAnnot instanceof List ? 'none' : 'yes'}" > ${meta.prefixRNA}.summaryRNA.stubinfo.txt
    """
}
