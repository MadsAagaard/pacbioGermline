#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


////////////////////////////////////////////
/////// ------- RUN BOOKKEEPING ------- ////
////////////////////////////////////////////

process write_input_summary {
    label "low"
    publishDir ( 
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"}, 
        mode: 'copy', 
        pattern: "*.txt")

    publishDir (
        path: {"${params.lrsDocuments}/summaryData/allSamples/"}, 
        mode: 'copy', 
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    cat > ${params.ssBase}.${params.inputReadSet}.input.allSamples.summary.txt << 'EOF'
    ${summary_ch}
    """
}

process write_dropped_samples_summary {
    label "low"
    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"}, 
        mode: 'copy', 
        pattern: "*.txt")
    publishDir (
        path: {"${params.lrsDocuments}/summaryData/droppedSamples/"}, 
        mode: 'copy', 
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    cat > ${params.ssBase}.${params.inputReadSet}.dropped.samples.summary.txt << 'EOF'
    ${summary_ch}
    """
}

process write_analyzed_samples_summary {
    label "low"
    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"}, 
        mode: 'copy', 
        pattern: "*.txt")
    publishDir (
        path: {"${params.lrsDocuments}/summaryData/analyzedSamples/"}, 
        mode: 'copy', 
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    cat > ${params.ssBase}.${params.inputReadSet}.analyzed.samples.summary.txt << 'EOF'
    ${summary_ch}
    """
}


////////////////////////////////////////////
/////// ------- PREPROCESS + ALN ------- ///
////////////////////////////////////////////

process create_fofn {
    label "low"
    tag "$meta.id"

    publishDir (
        path: {"${params.outBase(meta)}/documents/"}, 
        mode: 'copy', 
        pattern: "*.fofn", 
        overwrite: true)

    input:
    tuple val(meta), path(data)   // unmapped BAMs

    output:
    tuple val(meta), path("${meta.id}.fofn"),                emit: allReads
    tuple val(meta), path("${meta.id}.hifi_reads.fofn"),     emit: hifiReads
    tuple val(meta), path("${meta.id}.fail_reads.fofn"),     emit: failReads

    script:

    def failCheck = params.hifiReads
        ? "true"
        : """if [ ! -s ${meta.id}.fail_reads.fofn ]; then
                 echo "ERROR: no fail_reads uBAMs found for sample ${meta.id}." >&2
                 echo "       Use --hifiReads to run HiFi-only, or fix the input glob." >&2
                 exit 1
             fi"""
    """
    realpath ${data} > ${meta.id}.fofn

    grep 'fail_reads' ${meta.id}.fofn > ${meta.id}.fail_reads.fofn || true
    grep 'hifi_reads' ${meta.id}.fofn > ${meta.id}.hifi_reads.fofn || true

    if [ ! -s ${meta.id}.hifi_reads.fofn ]; then
        echo "ERROR: no hifi_reads uBAMs found for sample ${meta.id}." >&2
        exit 1
    fi

    ${failCheck}
    """
}

process inputFiles_symlinks_ubam {
    label "low"
    tag "$meta.id"
    publishDir (
        path: {"${params.outBase(meta)}/documents/inputSymlinks/"}, 
        mode: 'symlink', 
        pattern: '*.{bam,pbi}', 
        overwrite: true)

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path(data)

    script:
    """
    """
}

process symlinks_ubam_dropped {
    label "low"
    tag "$meta.id"

    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/dropped_samples_ubam_symlinks/"}, 
        mode: 'symlink', 
        pattern: '*.{bam,pbi}', 
        overwrite: true)

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path(data)

    script:
    """
    """
}

process pbmm2_align_hifi {
    label "veryHigh"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    input:
    tuple val(meta), path(fofn)

    output:
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.pbmm2.merged.bam"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.pbmm2.merged.bam.bai"), emit: bam

    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genomeMmi} \
    ${fofn} \
    ${meta.id}.${params.genomeVersion}.${params.tagHifi}.pbmm2.merged.bam
    """
}

process pbmm2_align_failedOnly {
    label "intermediate"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    input:
    tuple val(meta), path(fofn)

    output:
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagFail}.pbmm2.bam"),
          path("${meta.id}.${params.genomeVersion}.${params.tagFail}.pbmm2.bam.bai"), emit: bam

    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genomeMmi} \
    ${fofn} \
    ${meta.id}.${params.genomeVersion}.${params.tagFail}.pbmm2.bam
    """
}


////////////////////////////////////////////
/////// ------- SMALL VARIANTS ------- /////
////////////////////////////////////////////

process deepvariant {
    label "veryHigh"
    tag "$meta.id"

    publishDir (
        path: {"${params.lrsStorage}/deepVariant/gvcf/"}, 
        mode: 'copy',
        pattern: "*.deepVariant.g.vcf.*")
    
    publishDir ( 
        path: {"${params.outBase(meta)}/SNV_and_INDELs/gvcf/"},
        mode: 'copy',
        pattern: "*.deepVariant.g.vcf.*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.deepVariant.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.deepVariant.vcf.gz.tbi"),   emit: dv_vcf

    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.deepVariant.g.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.deepVariant.g.vcf.gz.tbi"), emit: dv_gvcf

    script:
    """
    singularity run -B ${params.sBind} ${params.deepvariantSif} /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref=${params.genomeFasta} \
    --reads=${data.hifiBam} \
    --output_vcf=${meta.id}.${params.genomeVersion}.${params.tagHifi}.deepVariant.vcf.gz \
    --output_gvcf=${meta.id}.${params.genomeVersion}.${params.tagHifi}.deepVariant.g.vcf.gz \
    --num_shards=${task.cpus}
    """
}

process glNexus_jointCall {
    label "high"
    tag "$meta.caseID"
    conda "${params.condaEnvs.glnexus}"

    publishDir (
        path: {"${params.outBase(meta)}/jointCalls/"}, 
        mode: 'copy', 
        pattern: "*.jointCall.*")
    
    publishDir (
        path: {"${params.outBase(meta)}/documents/"},   
        mode: 'copy', 
        pattern: "*.manifest")

    input:
    tuple val(meta), path(manifest)

    output:
    tuple val(meta),
          path("${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.vcf.gz"),
          path("${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.vcf.gz.tbi"),          emit: glnexus_vcf
    
    tuple val(meta),
          path("${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.WES_ROI.vcf.gz"),
          path("${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.WES_ROI.vcf.gz.tbi"),  emit: glnexus_wes_roi_vcf
    
    tuple val(meta), path("${manifest}"), emit: glnexus_manifest

    script:

    """
    glnexus_cli \
    --config DeepVariant \
    --threads ${task.cpus} \
    --list ${manifest} > ${meta.caseID}.glnexus.bcf

    bcftools view -Oz -o ${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.vcf.gz ${meta.caseID}.glnexus.bcf
    bcftools index -t ${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.vcf.gz

    bcftools view -R ${params.roiBed} \
    ${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.vcf.gz \
    -Oz -o ${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.WES_ROI.vcf.gz

    bcftools index -t ${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.deepVariant.jointCall.WES_ROI.vcf.gz
    """
}


///////////////////////////////////////////////////
////// ------- PHASING ------- ////////////////////
///////////////////////////////////////////////////

process hiPhaseTwoAln {
    tag "$meta.id"
    label "intermediate"
    conda "${params.condaEnvs.hiphase}"

    publishDir (
        path: {"${params.outBase(meta)}/alignments/"},                  
        mode: 'copy', 
        pattern: "*.${params.tagHifi}.hiphase.ba*")

    publishDir (
        path: {"${params.outBase(meta)}/alignments/failedReads/"},      
        mode: 'copy', 
        pattern: "*.${params.tagFail}.hiphase.ba*")

    publishDir (
        path: {"${params.outBase(meta)}/SNV_and_INDELs/"},              
        mode: 'copy', 
        pattern: "*.hiphase.deepvariant.*")

    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/"},          
        mode: 'copy', 
        pattern: "*.hiphase.sawfish.*")

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/diseaseSTRs/"}, 
        mode: 'copy', 
        pattern: "*.hiphase.trgt4.*")
    


    publishDir (
        path: {"${params.lrsStorage}/deepVariant/vcfs/"},                        
        mode: 'copy', 
        pattern: "*.hiphase.deepvariant.vcf.*")

/*

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT5_all/adotto/"}, 
        mode: 'copy', 
        pattern: "*.hiphase.trgt5.adotto.sorted.*")


    publishDir (
        path: {"${params.outBase(meta)}/QC/hiphaseStats/"},
        mode: 'copy',
        pattern: "*.hiphase.{stats,blocks,summary}.{csv,tsv}"
        )
*/

    input:
    tuple val(meta), val(data)
    // data keys: hifiBam, hifiBai, failBam, failBai, dvVcf, sawfishVcf, str4Vcf, str5AdottoVcf

    output:
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.bam"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.bam.bai"),                            emit: hifi_bam

    // optional: absent in --hifiReads mode, where there is no fail alignment.
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagFail}.hiphase.bam"),
          path("${meta.id}.${params.genomeVersion}.${params.tagFail}.hiphase.bam.bai"), optional: true,            emit: fail_bam

    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.deepvariant.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.deepvariant.vcf.gz.tbi"),             emit: dv_vcf

    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.deepvariant.WES_ROI.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.deepvariant.WES_ROI.vcf.gz.tbi"),     emit: dv_wes_roi_vcf

    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.sawfish.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.sawfish.vcf.gz.tbi"),                 emit: sawfish_vcf

    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.hiphase.trgt4.STRchive.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.hiphase.trgt4.STRchive.sorted.vcf.gz.tbi"),    emit: trgt4_disease_vcf

    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.hiphase.trgt5.adotto.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.hiphase.trgt5.adotto.sorted.vcf.gz.tbi"),      emit: trgt5_adotto_vcf
/*
    tuple val(meta),
          path("${prefix}.hiphase.stats.csv"),
          path("${prefix}.hiphase.blocks.tsv"),
          path("${prefix}.hiphase.summary.tsv"), emit: hiphase_stats
*/


    script:
    def prefix  = "${meta.id}.${params.genomeVersion}"

    def bamArgs = "--bam ${data.hifiBam} --output-bam ${prefix}.${params.tagHifi}.hiphase.bam"
    if (data.failBam) {
        bamArgs += " --bam ${data.failBam} --output-bam ${prefix}.${params.tagFail}.hiphase.bam"
    }
    """
    hiphase \
    ${bamArgs} \
    --vcf ${data.dvVcf} \
    --output-vcf ${prefix}.${params.tagHifi}.hiphase.deepvariant.vcf.gz \
    --vcf ${data.sawfishVcf} \
    --output-vcf ${prefix}.${params.tagHifi}.hiphase.sawfish.vcf.gz \
    --vcf ${data.str4Vcf} \
    --output-vcf ${prefix}.${params.strTag}.hiphase.trgt4.STRchive.sorted.vcf.gz \
    --vcf ${data.str5AdottoVcf} \
    --output-vcf ${prefix}.${params.strTag}.hiphase.trgt5.adotto.sorted.vcf.gz \
    --reference ${params.genomeFasta} \
    --threads ${task.cpus} \
   #--stats-file ${prefix}.hiphase.stats.csv \
    #--blocks-file ${prefix}.hiphase.blocks.tsv \
    #--summary-file ${prefix}.hiphase.summary.tsv \
    --io-threads ${task.cpus}

    bcftools index -t -f ${prefix}.${params.tagHifi}.hiphase.deepvariant.vcf.gz
    bcftools index -t -f ${prefix}.${params.tagHifi}.hiphase.sawfish.vcf.gz
    bcftools index -t -f ${prefix}.${params.strTag}.hiphase.trgt4.STRchive.sorted.vcf.gz
    bcftools index -t -f ${prefix}.${params.strTag}.hiphase.trgt5.adotto.sorted.vcf.gz

    ${params.gatkExec} SelectVariants \
    -R ${params.genomeFasta} \
    -V ${prefix}.${params.tagHifi}.hiphase.deepvariant.vcf.gz \
    -L ${params.roiBed} \
    -O ${prefix}.${params.tagHifi}.hiphase.deepvariant.WES_ROI.vcf.gz
    """
}


///////////////////////////////////////////////////
////// ------- CNV AND STRUCTURAL VARIANTS ------ /
///////////////////////////////////////////////////

process sawFish2 {
    tag "$meta.id"
    label "high"
    conda "${params.condaEnvs.sawfish2}"

    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/${meta.id}.sawfishSV/supportingFiles/"}, 
        mode: 'copy',
        pattern: "*.{bedgraph,bw,json,json.gz}")

    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/${meta.id}.sawfishSV/"},
        mode: 'copy',
        pattern: "${meta.id}.sawfishDiscover")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.sawfishSV.*")
    tuple val(meta), path("*.sawfishSV.vcf.gz"), path("*.sawfishSV.vcf.gz.tbi"),      emit: sv_vcf
    tuple val(meta), path("${meta.id}.sawfishDiscover"), val("${data.hifiBam}"),      emit: sv_discover_dir2
    tuple val(meta), path("*.sawfishSV.supporting_reads.json.gz"),                    emit: sv_supporting_reads
    tuple val(meta), path("${meta.id}.sawfishSV/"),                                   emit: sawfish_out_dir

    script:
    def exclude = params.genome == "hg38" ? "--cnv-excluded-regions ${params.sawfishExcludeRegions}" : ""
    // FAIL-CLOSED: throws rather than silently assuming XX.
    def sex     = params.sawfishExpectedCn(meta)
    def prefix  = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    sawfish discover \
    --threads ${task.cpus} \
    --ref ${params.genomeFasta} \
    --bam ${data.hifiBam} \
    $exclude \
    $sex \
    --output-dir ${meta.id}.sawfishDiscover

    sawfish joint-call \
    --threads ${task.cpus} \
    --report-supporting-reads \
    --sample ${meta.id}.sawfishDiscover \
    --output-dir ${meta.id}.sawfishSV

    mv ${meta.id}.sawfishSV/genotyped.sv.vcf.gz            ${prefix}.sawfishSV.vcf.gz
    mv ${meta.id}.sawfishSV/genotyped.sv.vcf.gz.tbi        ${prefix}.sawfishSV.vcf.gz.tbi
    mv ${meta.id}.sawfishSV/supporting_reads.json.gz       ${prefix}.sawfishSV.supporting_reads.json.gz
    mv ${meta.id}.sawfishSV/samples/*/gc_bias_corrected_depth.bw ${prefix}.sawfishSV.gc_bias_corrected_depth.bw
    mv ${meta.id}.sawfishSV/samples/*/depth.bw             ${prefix}.sawfishSV.depth.bw
    mv ${meta.id}.sawfishSV/samples/*/copynum.bedgraph     ${prefix}.sawfishSV.copynum.bedgraph
    mv ${meta.id}.sawfishSV/samples/*/copynum.summary.json ${prefix}.sawfishSV.copynum.summary.json
    """
}

process svdb_SawFish {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.svdb}"

    publishDir (
        path: "${params.lrsStorage}/structuralVariants/sawfish/", 
        mode: 'copy',
        pattern: "*.sawfishSV.hiphase.svdb.vcf*")
    
    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/"},
        mode: 'copy',
        pattern: "*.sawfishSV.hiphase.svdb.*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.sawfishSV.hiphase.svdb.*")
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz.tbi"), emit: sawfishAF10

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    svdb --query \
    --query_vcf ${data.sawfish_vcf} \
    --sqdb ${params.sawfishSvdb} > ${prefix}.sawfishSV.hiphase.svdb.vcf

    bgzip ${prefix}.sawfishSV.hiphase.svdb.vcf
    bcftools index -t ${prefix}.sawfishSV.hiphase.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${prefix}.sawfishSV.hiphase.svdb.vcf.gz \
    -Oz -o ${prefix}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz

    bcftools index -t ${prefix}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz
    """
}

process sawFish2_jointCall_all {
    label "high"
    conda "${params.condaEnvs.sawfish2}"

    input:
    val(x)

    output:
    tuple path("*.sawfishSV_jointCall.vcf.gz"), path("*.sawfishSV_jointCall.vcf.gz.tbi"), emit: sv_jointCall_vcf

    script:
    """
    sawfish joint-call \
    --threads ${task.cpus} \
    ${x} \
    --output-dir ${params.rundir}.sawfishSV_jointCall

    mv ${params.rundir}.sawfishSV_jointCall/genotyped.sv.vcf.gz     ${params.rundir}.${params.genomeVersion}.${params.tagHifi}.sawfishSV_jointCall.vcf.gz
    mv ${params.rundir}.sawfishSV_jointCall/genotyped.sv.vcf.gz.tbi ${params.rundir}.${params.genomeVersion}.${params.tagHifi}.sawfishSV_jointCall.vcf.gz.tbi
    """
}

process svdb_sawFish2_jointCall_all {
    label "low"
    conda "${params.condaEnvs.svdb}"

    publishDir (
        path: "${params.outputDirBase}/jointCalls_All/",
        mode: 'copy',
        pattern: "*_jointCall.svdb.*")

    input:
    tuple path(vcf), path(idx)

    output:
    path("*_jointCall.svdb.*")

    script:
    def prefix = "${params.rundir}.${params.genomeVersion}.${params.tagHifi}"
    """
    svdb --query \
    --query_vcf ${vcf} \
    --sqdb ${params.sawfishSvdb} > ${prefix}.sawfishSV_jointCall.svdb.vcf

    bgzip ${prefix}.sawfishSV_jointCall.svdb.vcf
    bcftools index -t ${prefix}.sawfishSV_jointCall.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${prefix}.sawfishSV_jointCall.svdb.vcf.gz \
    -Oz -o ${prefix}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz

    bcftools index -t ${prefix}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    """
}

process sawFish2_jointCall_caseID {
    tag "$meta.caseID"
    label "high"
    conda "${params.condaEnvs.sawfish2}"

    // FIX: pattern was "*.sawfishSV.hiphase.svdb.*", which never matched this
    // process's outputs — nothing was ever published.
    
    publishDir (
        path: {"${params.outBase(meta)}/jointCalls/"},
        mode: 'copy',
        pattern: "*.sawfishSV_jointCall.vcf.gz*")
    
    publishDir (
        path: {"${params.outBase(meta)}/documents/"}, 
        mode: 'copy',
        pattern: "*.csv")

    input:
    tuple val(meta), path(manifest)

    output:
    tuple val(meta), path("*.sawfishSV_jointCall.vcf.gz"), path("*.sawfishSV_jointCall.vcf.gz.tbi"), emit: sv_jointCall_caseID_vcf
    tuple val(meta), path("${manifest}")

    script:
    def prefix = "${meta.caseID}.${params.genomeVersion}.${params.tagHifi}"
    """
    sawfish joint-call \
    --threads ${task.cpus} \
    --sample-csv ${manifest} \
    --output-dir ${meta.caseID}.sawfishSV_jointCall

    mv ${meta.caseID}.sawfishSV_jointCall/genotyped.sv.vcf.gz     ${prefix}.sawfishSV_jointCall.vcf.gz
    mv ${meta.caseID}.sawfishSV_jointCall/genotyped.sv.vcf.gz.tbi ${prefix}.sawfishSV_jointCall.vcf.gz.tbi
    """
}

process svdb_sawFish2_jointCall_caseID {
    tag "$meta.caseID"
    label "low"
    conda "${params.condaEnvs.svdb}"

    publishDir (
        path: {"${params.outBase(meta)}/jointCalls/"},
        mode: 'copy',
        pattern: "*_jointCall.svdb.*")

    input:
    tuple val(meta), path(vcf), path(idx)

    output:
    path("*_jointCall.svdb.*")
    tuple val(meta),
          path("${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz"),
          path("${meta.caseID}.${params.genomeVersion}.${params.tagHifi}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz.tbi"), emit: sawfish_caseID_AF10

    script:
    def prefix = "${meta.caseID}.${params.genomeVersion}.${params.tagHifi}"
    """
    svdb --query \
    --query_vcf ${vcf} \
    --sqdb ${params.sawfishSvdb} > ${prefix}.sawfishSV_jointCall.svdb.vcf

    bgzip ${prefix}.sawfishSV_jointCall.svdb.vcf
    bcftools index -t ${prefix}.sawfishSV_jointCall.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${prefix}.sawfishSV_jointCall.svdb.vcf.gz \
    -Oz -o ${prefix}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz

    bcftools index -t ${prefix}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    """
}

process svTopo {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.svtopo}"

    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/SVtopo/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.svtopo_out/")

    script:
    def exclude = params.genome == "hg38" ? "--exclude-regions ${params.sawfishExcludeRegions}" : ""
    """
    mkdir ${meta.id}.svtopo_out

    svtopo \
    --bam ${data.bam} \
    --vcf ${data.sawfish_vcf} \
    --variant-readnames ${data.sawfish_reads} \
    --prefix ${meta.id} \
    $exclude \
    --svtopo-dir ${meta.id}.svtopo_out/

    svtopovz \
    --svtopo-dir ${meta.id}.svtopo_out/ \
    --genes ${params.gencodeGtf} \
    --image-type jpg

    mv ${meta.id}.svtopo_out/index.html ${meta.id}.svtopo_out/${meta.id}.sawfishSV.svtopo.html
    """
}

process svTopo_filtered {
    tag "$meta.id"
    label "high"
    conda "${params.condaEnvs.svtopo}"

    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/SVtopo_filtered/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.svtopo_out/")

    script:
    def exclude = params.genome == "hg38" ? "--exclude-regions ${params.sawfishExcludeRegions}" : ""
    """
    mkdir ${meta.id}.svtopo_out

    svtopo \
    --bam ${data.bam} \
    --vcf ${data.sawfish10_vcf} \
    --variant-readnames ${data.sawfish_reads} \
    --prefix ${meta.id} \
    $exclude \
    --svtopo-dir ${meta.id}.svtopo_out/

    svtopovz \
    --svtopo-dir ${meta.id}.svtopo_out/ \
    --genes ${params.gencodeGtf} \
    --image-type jpg

    mv ${meta.id}.svtopo_out/index.html ${meta.id}.svtopo_out/${meta.id}.sawfishSV.svtopo.html
    """
}


///////////////////////////////////////////////////
/////// --- PSEUDOGENES, REPEATS, MITO --- ////////
///////////////////////////////////////////////////

process mitorsaw {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.mitorsaw}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/mitochondrialVariants/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.mitorsaw.*")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    mitorsaw haplotype \
    --reference ${params.genomeFasta} \
    --bam ${data.bam} \
    --output-vcf ${prefix}.mitorsaw.vcf.gz \
    --output-hap-stats ${prefix}.mitorsaw.hapstats.json
    """
}


process trgt4_diseaseSTRs {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.trgt4}"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/bam"},
        mode: 'copy',
        pattern: "*.sorted.ba*")
    
    publishDir (
        path: "${params.lrsStorage}/repeatExpansions/TRGT/diseaseSTRs/",
        mode: 'copy',
        pattern: "*.sorted.vcf.*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive.sorted.vcf.gz.tbi"), emit: str4_vcf
    
    tuple val(meta), path("*.sorted.*")
   
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive.sorted.bam"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive.sorted.bam.bai"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive.sorted.vcf.gz.tbi"), emit: trgt_full

    script:
    def karyotype = params.karyotype(meta)
    def prefix    = "${meta.id}.${params.genomeVersion}.${params.strTag}.trgt4.STRchive"
    def merged    = "${meta.id}.${params.genomeVersion}.AllReads.pbmm2.merged.bam"

    """
    # ---- merge HiFi + failed reads into a single allReads BAM (TRGT4 has no --fail-reads)
    samtools merge \
    --threads ${task.cpus} \
    -c -p \
    -o ${merged} \
    ${data.hifiBam} ${data.failBam}

    samtools index -@ ${task.cpus} ${merged}

    trgt genotype \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --reads ${merged} \
    $karyotype \
    --min-read-quality 0 \
    --output-prefix ${prefix}

    bcftools sort -Oz -o ${prefix}.sorted.vcf.gz ${prefix}.vcf.gz
    bcftools index -t ${prefix}.sorted.vcf.gz

    samtools sort -o ${prefix}.sorted.bam ${prefix}.spanning.bam
    samtools index ${prefix}.sorted.bam
    """
}

process trgt4_diseaseSTRs_plots {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.trgt4}"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/Plots/"},
        mode: 'copy',
        pattern: "*.{pdf,png,svg}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.{pdf,png,svg}")

    script:
    def geneList = params.puretargetPlotGenes.join(' ')
    def prefix   = "${meta.id}.${params.genomeVersion}.${params.strTag}"
    """
    for gene in ${geneList}; do
      trgt plot \
      --genome ${params.genomeFasta} \
      --repeats ${params.trDiseaseCatalog} \
      --vcf ${data.vcf} \
      --spanning-reads ${data.bam} \
      --repeat-id \$gene \
      --squished \
      -o ${prefix}.\$gene.allele.pdf

      trgt plot \
      --genome ${params.genomeFasta} \
      --repeats ${params.trDiseaseCatalog} \
      --vcf ${data.vcf} \
      --spanning-reads ${data.bam} \
      --repeat-id \$gene \
      --plot-type waterfall \
      -o ${prefix}.\$gene.waterfall.pdf
    done
    """
}

process trgt4_diseaseSTRs_plots_meth {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/METHplots/"},
        mode: 'copy',
        pattern: "*.{pdf,png,svg}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.{pdf,png,svg}")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.strTag}"
    """
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --show meth \
    --squished \
    --max-allele-reads 75 \
    -o FXS_FMR1.${prefix}.METH.alleleSquished.pdf

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --plot-type waterfall \
    --show meth \
    --max-allele-reads 75 \
    -o FXS_FMR1.${prefix}.METH.waterfall.pdf
    """
}

process trgt5_all_adotto {
    tag "$meta.id"
    label "intermediate"
    conda "${params.condaEnvs.trgt51}"

    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/repeatExpansions/TRGT5_all/adotto/"},
        mode: 'copy',
        pattern: "*.adotto.LPS.txt")

    publishDir (
        path: "${params.lrsStorage}/repeatExpansions/TRGT5_all/adotto_LPS/",
        mode: 'copy',
        pattern: "*.adotto.LPS.txt")

    publishDir (
        path: "${params.lrsStorage}/repeatExpansions/TRGT5_all/adotto/",
        mode: 'copy',
        pattern: "*.trgt5.adotto.sorted.vcf.*")




    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.adotto.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.adotto.sorted.vcf.gz.tbi"), emit: adotto_vcf
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.adotto.LPS.txt"),           emit: adotto_LPS

    script:
    def karyotype = params.karyotype(meta)
    def failReads = data.failBam ? "--fail-reads ${data.failBam}" : ""
    def prefix    = "${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.adotto"
    """
    trgt genotype \
    --genome ${params.genomeFasta} \
    --repeats ${params.trAdottoCatalog} \
    --reads ${data.hifiBam} \
    $failReads \
    $karyotype \
    --threads ${task.cpus} \
    --disable-bam-output \
    --output-prefix ${prefix}

    bcftools sort -Oz -o ${prefix}.sorted.vcf.gz ${prefix}.vcf.gz
    bcftools index -t ${prefix}.sorted.vcf.gz

    ${params.trgtLps} \
    --vcf ${prefix}.sorted.vcf.gz \
    --threads ${task.cpus} > ${prefix}.LPS.txt
    """
}

process trgt5_all_TRexplorer {
    tag "$meta.id"
    label "intermediate"
    conda "${params.condaEnvs.trgt51}"


    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/repeatExpansions/TRGT5_all/trexplorer/"},
        mode: 'copy',
        pattern: "*.sorted.vcf.*")

    publishDir (
        path: "${params.lrsStorage}/repeatExpansions/TRGT5_all/trexplorer/",
        mode: 'copy',
        pattern: "*.sorted.vcf.*")

    publishDir (
        path: "${params.lrsStorage}/repeatExpansions/TRGT5_all/trexplorer_LPS/",
        mode: 'copy',
        pattern: "*.LPS.txt")

    input:
    tuple val(meta), val(data)

    output:
    
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.trexplorer.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.trexplorer.sorted.vcf.gz.tbi"), emit: trexplorer_vcf
    
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.trexplorer.LPS.txt"),           emit: trexplorer_LPS

    script:
    def karyotype = params.karyotype(meta)
    def failReads = data.failBam ? "--fail-reads ${data.failBam}" : ""
    def prefix    = "${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.trexplorer"
    """
    trgt genotype \
    --genome ${params.genomeFasta} \
    --repeats ${params.trExplorerCatalog} \
    --reads ${data.hifiBam} \
    $failReads \
    $karyotype \
    --threads ${task.cpus} \
    --disable-bam-output \
    --output-prefix ${prefix}

    bcftools sort -Oz -o ${prefix}.sorted.vcf.gz ${prefix}.vcf.gz
    bcftools index -t ${prefix}.sorted.vcf.gz

    ${params.trgtLps} \
    --vcf ${prefix}.sorted.vcf.gz \
    --threads ${task.cpus} > ${prefix}.LPS.txt
    """
}

process trgt5_diseaseSTRs {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.trgt51}"

    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/repeatExpansions/TRGT5/bam"},
        mode: 'copy',
        pattern: "*.sorted.ba*")

    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/repeatExpansions/TRGT5/diseaseSTRs"},
        mode: 'copy',
        pattern: "*.sorted.vcf*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.sorted.*")
    tuple val(meta),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.STRchive.sorted.bam"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.STRchive.sorted.bam.bai"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.STRchive.sorted.vcf.gz"),
          path("${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.STRchive.sorted.vcf.gz.tbi"), emit: trgt_full

    script:
    def karyotype = params.karyotype(meta)
    def failReads = data.failBam ? "--fail-reads ${data.failBam}" : ""
    def prefix    = "${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5.STRchive"
    """
    trgt genotype \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --reads ${data.hifiBam} \
    $failReads \
    $karyotype \
    --output-prefix ${prefix}

    bcftools sort -Oz -o ${prefix}.sorted.vcf.gz ${prefix}.vcf.gz
    bcftools index -t ${prefix}.sorted.vcf.gz

    samtools sort -o ${prefix}.sorted.bam ${prefix}.spanning.bam
    samtools index ${prefix}.sorted.bam
    """
}

process trgt5_diseaseSTRs_plots {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.trgt51}"

    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/repeatExpansions/TRGT5/Plots/"},
        mode: 'copy',
        pattern: "*.{pdf,png,svg}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.{pdf,png,svg}")

    script:
    def geneList = params.puretargetPlotGenes.join(' ')
    def prefix   = "${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5"
    """
    for gene in ${geneList}; do
      trgt plot \
      --genome ${params.genomeFasta} \
      --repeats ${params.trDiseaseCatalog} \
      --vcf ${data.vcf} \
      --spanning-reads ${data.bam} \
      --repeat-id \$gene \
      --squished \
      -o ${prefix}.\$gene.allele.pdf

      trgt plot \
      --genome ${params.genomeFasta} \
      --repeats ${params.trDiseaseCatalog} \
      --vcf ${data.vcf} \
      --spanning-reads ${data.bam} \
      --repeat-id \$gene \
      --plot-type waterfall \
      -o ${prefix}.\$gene.waterfall.pdf
    done
    """
}

process trgt5_diseaseSTRs_plots_meth {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt51}"

    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/repeatExpansions/TRGT5/METHplots/"}, 
        mode: 'copy',
        pattern: "*.{pdf,png,svg}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.{pdf,png,svg}")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.strTag}.trgt5"
    """
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --show meth \
    --squished \
    --max-allele-reads 75 \
    -o FXS_FMR1.${prefix}.METH.alleleSquished.pdf

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --plot-type waterfall \
    --show meth \
    --max-allele-reads 75 \
    -o FXS_FMR1.${prefix}.METH.waterfall.pdf
    """
}

process kivvi_d4z4 {
    tag "$meta.id"
    label "medium"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/Kivvi_D4Z4_contraction/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.kivviD4Z4")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    ${params.kivviDir}/kivvi \
    -r ${params.genomeFasta} \
    --bam ${data.bam} \
    -p ${prefix} \
    -o ${prefix}.kivviD4Z4 \
    d4z4
    """
}

process paraphase {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.paraphaseMinimap2}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/paraphase3/"}, 
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.paraphase/*")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    paraphase \
    -b ${data.bam} \
    --reference ${params.genomeFasta} \
    -t ${task.cpus} \
    -o ${prefix}.hiphase.paraphase

    python ${params.pythonScripts}/flatten_paraphaseSMN.py \
    --json ${prefix}.hiphase.paraphase/${meta.id}.paraphase.json \
    --out  ${prefix}.hiphase.paraphase/${prefix}.paraphase.flattened.tsv

    python ${params.pythonScripts}/flatten_paraphaseSMN.py \
    --json ${prefix}.hiphase.paraphase/${meta.id}.paraphase.json \
    --loci SMN1,PMS2,IKBKG \
    --out  ${prefix}.hiphase.paraphase/${prefix}.paraphase.flattened.SMN_PMS2_IKBKG.tsv
    """
}

process paraphase4 {
    tag "$meta.id"
    label "lowCPU"
    conda "${params.condaEnvs.paraphase40}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/paraphase4/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.paraphase/*")
    tuple val(meta), path("{prefix}.paraphase.json"), emit: paraphase_json

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    paraphase \
    -b ${data.bam} \
    --reference ${params.genomeFasta} \
    -t ${task.cpus} \
    -o ${prefix}.hiphase.paraphase

    cp ${prefix}.hiphase.paraphase/${meta.id}.paraphase.json {prefix}.paraphase.json
    """
}

process starphase {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.starphase}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/starphase/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.starphase.*")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    pbstarphase diplotype \
    --database ${params.starphaseDb} \
    --bam ${data.bam} \
    --reference ${params.genomeFasta} \
    --vcf ${data.dv_vcf} \
    --sv-vcf ${data.sawfish_vcf} \
    --pharmcat-tsv ${prefix}.starphase.diplotypes_for_pharmCAT.tsv \
    --output-calls ${prefix}.starphase.json

    java -jar ${params.pharmcatJar} \
    -po ${prefix}.starphase.diplotypes_for_pharmCAT.tsv \
    -vcf ${data.dv_vcf} \
    -bf ${prefix}.starphase.pharmCAT \
    -o .
    """
}


///////////////////////////////////////////////////
/// TRIO: variant prioritization — Exomiser 14.1.0 /
///////////////////////////////////////////////////

process exo14_2508_exome {
    label "medium"
    tag "$meta.caseID"

    publishDir (
        path: {"${params.outBase(meta)}/exomiser14_2508/exomiser/"},
        mode: 'copy')

    publishDir (
        path: {"${params.outBase(meta)}/documents/"},
        mode: 'copy',
        pattern: "*.{hpo.txt,yml,ped}")

    input:
    tuple val(meta), path(vcf), path(idx), path(hpo), path(samplesheet)

    output:
    path("*.{html,tsv,vcf,json,hpo.txt,yml,ped}")

    script:
    """
    python3 ${params.pythonScripts}/make_ped_and_family_v2.py \
    --samplesheet ${samplesheet} \
    --vcf ${vcf} \
    --caseid ${meta.caseID} \
    --hpo ${hpo}

    java -jar ${params.exomiserJar} \
    --sample ${meta.caseID}-family.yml \
    --analysis ${params.exomiserExomeYml} \
    --spring.config.location=${params.exomiserSpringDir}

    mv results/* .
    mv exomiser_tmp.html         ${meta.caseID}.exo14_2508.html
    mv exomiser_tmp.variants.tsv ${meta.caseID}.exo14_2508.Variants.tsv
    mv exomiser_tmp.genes.tsv    ${meta.caseID}.exo14_2508.Genes.tsv
    mv exomiser_tmp.json         ${meta.caseID}.exo14_2508.json
    """
}

process exo14_2508_genome {
    label "medium"
    tag "$meta.caseID"

    publishDir (
        path: {"${params.outBase(meta)}/exomiser14_2508/genomiser/"},
        mode: 'copy')

    input:
    tuple val(meta), path(vcf), path(idx), path(hpo), path(samplesheet)

    output:
    path("*.{html,tsv,vcf,json}")

    script:
    """
    python3 ${params.pythonScripts}/make_ped_and_family_v2.py \
    --samplesheet ${samplesheet} \
    --vcf ${vcf} \
    --caseid ${meta.caseID} \
    --hpo ${hpo}

    java -jar ${params.exomiserJar} \
    --sample ${meta.caseID}-family.yml \
    --analysis ${params.exomiserGenomeYml} \
    --spring.config.location=${params.exomiserSpringDir}

    mv results/* .
    mv exomiser_tmp.html         ${meta.caseID}.genomiser14_2508.html
    mv exomiser_tmp.variants.tsv ${meta.caseID}.genomiser14_2508.Variants.tsv
    mv exomiser_tmp.genes.tsv    ${meta.caseID}.genomiser14_2508.Genes.tsv
    mv exomiser_tmp.json         ${meta.caseID}.genomiser14_2508.json
    """
}

process exo14_2508_SV {
    label "medium"
    tag "$meta.caseID"

    publishDir (
        path: {"${params.outBase(meta)}/exomiser14_2508/exomiserStructuralVariants/"},
        mode: 'copy')

    input:
    tuple val(meta), path(vcf), path(idx), path(hpo), path(samplesheet)

    output:
    path("*.{html,tsv,vcf,json,hpo.txt,yml,ped}")

    script:
    """
    zcat ${vcf} | sed 's/^##fileformat=VCFv4\\.4/##fileformat=VCFv4.2/' | bgzip > ${meta.caseID}.sawfish.forExomiser.vcf.gz

    python3 ${params.pythonScripts}/make_ped_and_family_v2.py \
    --samplesheet ${samplesheet} \
    --vcf ${meta.caseID}.sawfish.forExomiser.vcf.gz \
    --caseid ${meta.caseID} \
    --hpo ${hpo}

    java -jar ${params.exomiserJar} \
    --sample ${meta.caseID}-family.yml \
    --analysis ${params.exomiserExomeYml} \
    --spring.config.location=${params.exomiserSpringDir}

    mv results/* .
    mv exomiser_tmp.html         ${meta.caseID}.SVs.exo14_2508.html
    mv exomiser_tmp.variants.tsv ${meta.caseID}.SVs.exo14_2508.Variants.tsv
    mv exomiser_tmp.genes.tsv    ${meta.caseID}.SVs.exo14_2508.Genes.tsv
    mv exomiser_tmp.json         ${meta.caseID}.SVs.exo14_2508.json
    """
}


///////////////////////////////////////////////////
/////// ------- METHYLATION ------- ///////////////
///////////////////////////////////////////////////

process pbCPGtools {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.pbCPGtools}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/methylation/BigWigBed/"},
        mode: 'copy',
        pattern: "*.methylation.{hap1,hap2,combined}.*")

    publishDir (
        path: "${params.lrsStorage}/methylation/pbCpGtools/${meta.id}/",
        mode: 'copy',
        pattern: "*.bed.*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.methylation*")

    script:
    """
    aligned_bam_to_cpg_scores \
    --bam ${data.bam} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.tagHifi}.hiphase.methylation
    """
}

process methBat {
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.methbat}"

    publishDir (
        path:{"${params.outBase(meta)}/specialAnalysis/methylation/"},
        mode: 'copy')

    publishDir (
        path: "${params.lrsStorage}/methylation/methBatProfiles/",
        mode: 'copy',
        pattern: "*.profile")

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("*.met.*")

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    methbat segment \
    --input-prefix ${prefix}.hiphase.methylation \
    --output-prefix ${prefix}.met.segments

    methbat profile \
    --input-prefix ${prefix}.hiphase.methylation \
    --input-regions ${params.methBackground} \
    --output-region-profile ${prefix}.met.profile

    methbat profile \
    --input-prefix ${prefix}.hiphase.methylation \
    --input-regions ${params.methBackgroundLocal} \
    --output-region-profile ${prefix}.met.profileLOCAL

    methbat report \
    --input-prefix ${prefix}.hiphase.methylation \
    --input-regions ${params.methICRegions} \
    --output-report ${prefix}.met.imprintingReport.tsv
    """
}

process methBatNEW_pileup {
    tag "$meta.id"
    label "intermediateCPU"
    conda "${params.condaEnvs.methbatV1}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/methylationNEW/5mC_pileup/"},
        mode: 'copy',
        pattern: "*.5mC.bed.*")

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/methylationNEW/5mC_bedgraphs/"},
        mode: 'copy',
        pattern: "*.5mC.bedgraph.*")

    publishDir (
        path: "${params.lrsStorage}/methylationNEW/5mC_pileup/",
        mode: 'copy',
        pattern: "*.5mC.bed.*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.met.*"), path("*.5mC.bedgraph.*")
    tuple val(meta), path("*.5mC.bed.gz"),  path("*.5mC.bed.gz.tbi"),  emit: mC5
    tuple val(meta), path("*.5hmC.bed.gz"), path("*.5hmC.bed.gz.tbi"), emit: hmC5
    tuple val(meta), path("*.6mA.bed.gz"),  path("*.6mA.bed.gz.tbi"),  emit: mA6

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    methbat pileup \
    --threads ${task.cpus} \
    --input-bam ${data.bam} \
    --output-prefix ${prefix}.met.pileup

    zgrep "Total" ${prefix}.met.pileup.5mC.bed.gz | \
    cut -f 1-3,7 | \
    bgzip > ${prefix}.5mC.bedgraph.gz

    tabix -p bed ${prefix}.5mC.bedgraph.gz
    """
}

process methBatNEW_profile_single {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.methbatV1}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/methylationNEW/"},
        mode: 'copy',
        pattern: "*.met.5mC.*")

    publishDir (
        path: "${params.lrsStorage}/methylationNEW/5mC_CGI_profiles/",
        mode: 'copy',
        pattern: "*.profile.tsv")

    input:
    tuple val(meta), path(data), path(tbi)

    output:
    tuple val(meta), 
        path("*.5mC.cpgIslands.profile.tsv"),       emit: cgiProfile
    
    tuple val(meta), 
        path("*.5mC.segments.*"),                   emit: segments
    
    tuple val(meta),
         path("*.met.5mC.imprintingReport.tsv"),    emit: icReport

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    methbat profile \
    --input-regions ${params.methCpGRegions} \
    --input-pileup ${data} \
    --output-region-profile ${prefix}.met.5mC.cpgIslands.profile.tsv

    methbat segment \
    --input-pileup ${data} \
    --output-prefix ${prefix}.met.5mC.segments

    methbat report \
    --input-pileup ${data} \
    --input-regions ${params.methICRegions} \
    --output-report ${prefix}.met.5mC.imprintingReport.tsv
    """
}


///////////////////////////////////////////////////
/////// ------- QUALITY CONTROL ------- ///////////
///////////////////////////////////////////////////

process mosdepthROI {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.mosdepth}"

    publishDir (
        path: {"${params.outBase(meta)}/QC/mosdepth/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}_roi.*"), emit: mosdepth_roi
    tuple val(meta), path("*.region.dist.txt"),                        emit: multiqc

    script:
    def callable = params.genome == "hg38" ? "--by ${params.callableRoiBed}" : "--by 1000"
    """
    mosdepth \
    -t ${task.cpus} \
    $callable \
    ${meta.id}.${params.genomeVersion}_roi \
    ${data.hifiBam}
    """
}

process whatsHap_stats {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.whatshap}"

    publishDir (
        path: {"${params.outBase(meta)}/QC/whatsHap/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.whatshap.stats.tsv"), emit: multiqc

    script:
    """
    whatshap stats \
    ${data.dv_vcf} \
    --tsv=${meta.id}.whatshap.stats.tsv
    """
}

process cramino {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.cramino}"

    publishDir (
        path: {"${params.outBase(meta)}/QC/cramino/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.craminoQC.txt")

    script:
    """
    cramino \
    -t ${task.cpus} \
    --karyotype \
    --phased \
    ${data.bam} > ${meta.id}.${params.genomeVersion}.${params.tagHifi}.craminoQC.txt
    """
}

process nanoStat {
    tag "$meta.id"
    label "low"
    conda "${params.condaEnvs.nanostats}"

    publishDir (
        path:{"${params.outBase(meta)}/QC/nanoStat/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.tagHifi}.nanostat.txt"), emit: multiqc

    script:
    """
    NanoStat \
    -t ${task.cpus} \
    -n ${meta.id}.${params.genomeVersion}.${params.tagHifi}.nanostat.txt \
    --bam ${data.hifiBam}
    """
}

process multiQC {
    tag { params.layoutMode == 'jointAnalysis' ? meta.caseID : meta.id }
    label "low"
    conda "${params.condaEnvs.multiqc}"

    publishDir (
        path: {"${params.outBase(meta)}/QC/"},
        mode: 'copy')

    input:
    tuple val(meta), path(qcfiles)

    output:
    path("*MultiQC*.html")

    script:
    def reportName = (params.layoutMode == 'jointAnalysis') ? "${meta.caseID}.MultiQC.DNA.html" : "${meta.id}.MultiQC.DNA.html"
    """
    mkdir -p qc_in
    cp -L ${qcfiles} qc_in/

    multiqc \
    -c ${params.multiqcConfig} \
    -f -q qc_in \
    -n ${reportName}
    """
}

process collect_germline_summary {
    label "low"
    tag "$meta.id"
    conda "${params.condaEnvs.somaticSummaryEnv}"        // needs pyyaml only - its in our somaticSummary conda env

    publishDir (
        path: {"${params.outBase(meta)}/newToolsTest/customReports/"},
        mode: 'copy',
        pattern: "*.{yaml,json,htmls}")

    publishDir (
        patH: "${params.lrsStorageBase}/clinicalSummaries/germline/json/",
        mode: 'copy', 
        pattern: "*.clinical_summaryGermline.json")

    publishDir (
        patH: "${params.lrsStorageBase}/clinicalSummaries/germline/yaml/",
        mode: 'copy',
        pattern: "*.clinical_summaryGermline.yaml")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.clinical_summaryGermline.yaml"), emit: yaml
    tuple val(meta), path("${meta.id}.clinical_summaryGermline.json"), emit: json
    tuple val(meta), path("${meta.id}.clinical_summaryGermline.html"), emit: html

    script:
    def prefix = "${meta.id}.${params.genomeVersion}.${params.tagHifi}"
    """
    python3 ${params.germlineSummaryPy} \
        --case-id             ${prefix} \
        --npn                 ${meta.id} \
        --gender              ${meta.gender} \
        --testlist            ${meta.testlist} \
        --genome-version      ${params.genomeVersion} \
        --paraphase-json      ${data.paraphase} \
        --methbat-imprinting  ${data.methbat_imprinting} \
        --regions             smn,pms2,rccx \
        --html-template       ${params.germlineSummaryHtml} \
        --output              ${prefix}.clinical_summaryGermline
    """
}

