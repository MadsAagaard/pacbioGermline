#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
date2=new Date().format( 'yyMMdd HH:mm:ss' )
user="$USER"
runID="${date}.${user}"


log.info """\
======================================================
Clinical Genetics Vejle: PacBio LRS v2
======================================================
Genome       : $params.genome
Read set     : $readSet
RunID        : $runID
Script start : $date2
Genome FASTA : ${genome_fasta}
Archive RAW  : ${dataArchive}
OutputDir    : ${outputDir}
sbind        : ${s_bind}
min input GB : $params.minGB
"""


////////////////////////////////////////////
/////// ------- PREPROCESS + ALN ------- ///
////////////////////////////////////////////
process check_tmpdir {
    label "low"
    script:
    """
    echo "TMPDIR is: \$TMPDIR"
    df -h \$TMPDIR
    """
}

process write_input_summary {
    publishDir "${outputDir}/runInfo/${date}/", mode: 'copy', pattern: "*.txt"
    publishDir "${lrsDocuments}/run_summaries/", mode: 'copy', pattern: "*.txt"

    //publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/documents/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/documents/"}, mode: 'copy', pattern: "*.txt"
    input:
    val(summary_ch)

    output:
    path("*.txt")


    when:
    !params.intSS
    script:
    """
    cat > ${ssBase}.${readSet}.input.allSamples.summary.txt << 'EOF'
    ${summary_ch}
    """
}

process write_dropped_samples_summary {
    publishDir "${outputDir}/runInfo/${date}/", mode: 'copy', pattern: "*.txt"
    publishDir "${lrsDocuments}/dropped_samples/", mode: 'copy', pattern: "*.txt"
    input:
    val(summary_ch)

    output:
    path("*.txt")

    when:
    !params.intSS
    script:
    """
    cat > ${ssBase}.${readSet}.dropped.samples.summary.txt << 'EOF'
    ${summary_ch}
    """
}


process create_fofn {
    label "low"
    
    publishDir {params.groupedOutput ? \
    "${outputDir}/${meta.caseID}/documents/" : \
    "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/documents/"}, \
    mode: 'copy',pattern: '*.fofn'

    cpus 4
    input:
    tuple val(meta), path(data) //ubam

    output:
    tuple val(meta), path("${meta.id}.fofn")
    script:
    """
    `realpath ${data} > ${meta.id}.fofn`
    """
} 

process inputFiles_symlinks_ubam{
    label "low"
    
    publishDir {params.groupedOutput ? \
    "${outputDir}/${meta.caseID}/documents/inputSymlinks/" : \
    "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/documents/inputSymlinks/"}, \
     mode: 'symlink', pattern: '*.{bam,pbi}'

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
    
    publishDir "${outputDir}/runInfo/dropped_samples_ubam_symlinks/", mode: 'symlink', pattern: '*.{bam,pbi}'

    input:
    tuple val(meta), path(data)   

    output:
    tuple val(meta), path(data)
 
    script:
    """
    """
}

process pbmm2_align {
    label "veryHigh"
    tag "$meta.id"
    conda "${params.pbmm2}"
    
    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.pbmm2.bam"), path("${meta.id}.${genome_version}.${readSet}.pbmm2*bai"),  emit: bam
 
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${genome_mmi} \
    ${data[0]} \
    ${meta.id}.${genome_version}.${readSet}.pbmm2.bam
    """
}

process pbmm2_align_mergedData {
    label "veryHigh"
    tag "$meta.id"
    conda "${params.pbmm2}"

    input:
    tuple val(meta), path(fofn)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.pbmm2.merged.bam"), path("${meta.id}.${genome_version}.${readSet}.pbmm2.merged*bai"),  emit: bam


    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${genome_mmi} \
    ${fofn} \
    ${meta.id}.${genome_version}.${readSet}.pbmm2.merged.bam
    """
}


////////////////////////////////////////////
/////// ------- SMALL VARIANTS ------- /////
////////////////////////////////////////////

process deepvariant{
    label "veryHigh"
    tag "$meta.id"

    publishDir "${lrsStorage}/deepVariant/gvcf/", mode: 'copy', pattern: "*.deepVariant.g.vcf.*"    
   
    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/SNV_and_INDELs/gvcf/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/SNV_and_INDELs/gvcf/"}, mode: 'copy', pattern: "*.deepVariant.g.vcf.*"


    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.deepVariant.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.deepVariant.vcf.gz.tbi"), emit: dv_vcf

    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.deepVariant.g.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.deepVariant.g.vcf.gz.tbi"), emit: dv_gvcf    
    //path("${meta.id}.deepvariant.vcf_stats_report.txt")
    """
    singularity run -B ${s_bind} ${simgpath}/deepvariant190.sif /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref=${genome_fasta} \
    --reads=${data[0]} \
    --output_vcf=${meta.id}.${genome_version}.${readSet}.deepVariant.vcf.gz \
    --output_gvcf=${meta.id}.${genome_version}.${readSet}.deepVariant.g.vcf.gz \
    --num_shards=${task.cpus}
    """    
}

process glNexus_jointCall { 
    label "high"
    tag "$caseID"
    conda "${params.glnexus}"

    publishDir {params.groupedOutput ? "${outputDir}/${caseID}/jointCalls/" : "${outputDir}/jointCalls/"}, mode: 'copy', pattern: "*.jointCall.*"

    publishDir {params.groupedOutput ? "${outputDir}/${caseID}/documents/" : "${outputDir}/documents/"}, mode: 'copy', pattern: "*.manifest"

    input:
    tuple val(caseID), path(manifest)

    output:
    tuple val(caseID), path("${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.vcf.gz"), path("${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.vcf.gz.tbi"), emit: glnexus_vcf
    tuple val(caseID), path("${manifest}")
    tuple val(caseID), path("${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.WES_ROI.vcf.gz"), path("${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.WES_ROI.vcf.gz.tbi"),emit:glnexus_wes_roi_vcf
    
    script:
    """
    glnexus_cli \
    --config DeepVariant \
    --threads ${task.cpus} \
    --list ${manifest} > ${caseID}.glnexus.bcf

    bcftools view -Oz -o ${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.vcf.gz ${caseID}.glnexus.bcf
    bcftools index -t ${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.vcf.gz

    bcftools view -R ${ROI} ${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.vcf.gz -Oz -o ${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.WES_ROI.vcf.gz
    bcftools index -t ${caseID}.${genome_version}.${readSet}.deepVariant.jointCall.WES_ROI.vcf.gz

    """
}


///////////////////////////////////////////////////
////// -------PHASING ------- /////////////////////
///////////////////////////////////////////////////


process hiPhase {
    
    tag "$meta.id"
    label "medium"
    conda "${params.hiphase}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/alignments/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/alignments/"}, mode: 'copy', pattern: "*.hiphase.ba*"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/SNV_and_INDELs/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/SNV_and_INDELs/"}, mode: 'copy', pattern: "*.hiphase.deepvariant.*"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/TRGT/diseaseSTRs/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/TRGT/diseaseSTRs/"}, mode: 'copy', pattern: "*.hiphase.trgt4.*"
    
    publishDir "${lrsStorage}/alignment/BAM/", mode: 'copy', pattern:"*.hiphase.ba*"
    publishDir "${lrsStorage}/alignment/CRAM/", mode: 'copy', pattern:"*.hiphase.cra*"
    publishDir "${lrsStorage}/deepVariant/vcfs/", mode: 'copy', pattern:"*.hiphase.deepvariant.vcf.*"


    input:
    tuple val(meta), path(aln), path(vcf), path(sv), path(str)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.bam"), path("${meta.id}.${genome_version}.${readSet}.hiphase.bam.bai"), emit: hiphase_bam                                                
   
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.cram"), path("${meta.id}.${genome_version}.${readSet}.hiphase.cram.crai"), emit: hiphase_cram       

    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.vcf.gz.tbi"), emit: hiphase_dv_vcf

    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.WES_ROI.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.WES_ROI.vcf.gz.tbi")

    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.sawfish.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.hiphase.sawfish.vcf.gz.tbi"), emit: hiphase_sawfish_vcf
   
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.trgt4.STRchive.sorted.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.hiphase.trgt4.STRchive.sorted.vcf.gz.tbi"), emit: hiphase_trgt_vcf

    //topic channel:
   // tuple(val(task.process), eval('hiphase --version |head -n1 ||true')),emit:versions_ch,topic:'versions'

    script:
    """
    hiphase \
    --bam ${aln[0]} \
    --output-bam ${meta.id}.${genome_version}.${readSet}.hiphase.bam \
    --vcf ${vcf[0]} \
    --output-vcf ${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.vcf.gz \
    --vcf ${sv[0]} \
    --output-vcf ${meta.id}.${genome_version}.${readSet}.hiphase.sawfish.vcf.gz \
    --vcf ${str[0]} \
    --output-vcf ${meta.id}.${genome_version}.${readSet}.hiphase.trgt4.STRchive.sorted.vcf.gz \
    --reference ${genome_fasta} \
    --threads ${task.cpus} \
    --io-threads ${task.cpus}

    bcftools index -t ${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.vcf.gz

    ${gatk_exec} SelectVariants \
    -R ${genome_fasta} \
    -V  ${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.vcf.gz \
    -L ${ROI} \
    -O  ${meta.id}.${genome_version}.${readSet}.hiphase.deepvariant.WES_ROI.vcf.gz

    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${genome_version}.${readSet}.hiphase.cram  ${meta.id}.${genome_version}.${readSet}.hiphase.bam 

    samtools index ${meta.id}.${genome_version}.${readSet}.hiphase.cram


    """
}



///////////////////////////////////////////////////
////// -------CNV AND STRUCTURAL VARIANTS ------- /
///////////////////////////////////////////////////

process sawFish2{
    tag "$meta.id"
    label "high"
    conda "${params.sawfish2}"

    publishDir "${lrsStorage}/structuralVariants/sawfish/raw/", mode: 'copy', pattern:"*.sawfishSV.vcf.*"
    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/structuralVariants/${meta.id}.sawfishSV/supportingFiles/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/structuralVariants/${meta.id}.sawfishSV/supportingFiles/"}, mode: 'copy', pattern: "*.{bedgraph,bw}"


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path(data), path("*.sawfishSV.*")
    tuple val(meta), path("*.sawfishSV.vcf.gz"),path("*.sawfishSV.vcf.gz.tbi"),emit: sv_vcf
    path("${meta.id}.sawfishDiscover"), emit: sv_discover_dir
    tuple val(meta), path("${meta.id}.sawfishDiscover"), path("${data[0]}"), emit: sv_discover_dir2
    tuple val(meta), path("*.sawfishSV.supporting_reads.json.gz"), emit: sv_supporting_reads
    tuple val(meta), path("${meta.id}.sawfishSV/"), emit: sawfish_out_dir

    script:
    def exclude=params.genome=="hg38" ? "--cnv-excluded-regions ${cnv_exclude_sawfish}" : ""
    def sex = (meta.sex=="male"||meta.sex=="M"||meta.genderFile=="M") ? "--expected-cn ${sawfishMaleExpectedCN}" : "--expected-cn ${sawfishFemaleExpectedCN}"
   
    """
    sawfish discover \
    --threads ${task.cpus} \
    --ref ${genome_fasta} \
    --bam ${data[0]} \
    $exclude \
    $sex \
    --output-dir ${meta.id}.sawfishDiscover 

    sawfish joint-call \
    --threads ${task.cpus} \
    --report-supporting-reads \
    --sample ${meta.id}.sawfishDiscover \
    --output-dir ${meta.id}.sawfishSV 
    
    mv ${meta.id}.sawfishSV/genotyped.sv.vcf.gz ${meta.id}.${genome_version}.${readSet}.sawfishSV.vcf.gz

    mv ${meta.id}.sawfishSV/genotyped.sv.vcf.gz.tbi ${meta.id}.${genome_version}.${readSet}.sawfishSV.vcf.gz.tbi

    mv ${meta.id}.sawfishSV/supporting_reads.json.gz ${meta.id}.${genome_version}.${readSet}.sawfishSV.supporting_reads.json.gz

    mv ${meta.id}.sawfishSV/samples/*/gc_bias_corrected_depth.bw ${meta.id}.${genome_version}.${readSet}.sawfishSV.gc_bias_corrected_depth.bw

    mv ${meta.id}.sawfishSV/samples/*/depth.bw ${meta.id}.${genome_version}.${readSet}.sawfishSV.depth.bw

    mv ${meta.id}.sawfishSV/samples/*/copynum.bedgraph ${meta.id}.${genome_version}.${readSet}.sawfishSV.copynum.bedgraph
    """
}


process svdb_SawFish {
    tag "$meta.id"
    label "low"
    conda "${params.svdb}"

    publishDir "${lrsStorage}/structuralVariants/sawfish/svdb/", mode: 'copy',pattern: "*.sawfishSV.hiphase.svdb.vcf*"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/structuralVariants/vcfs/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/structuralVariants/vcfs/"}, mode: 'copy', pattern: "*.sawfishSV.hiphase.svdb.*"



    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.sawfishSV.hiphase.svdb.*")
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz"),path("${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz.tbi"), emit: sawfishAF10
    script:
    """
    svdb --query \
    --query_vcf ${data.sawfish_vcf} \
    --sqdb ${sawfish_sqdb} > ${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.vcf
    
    bgzip ${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.vcf
    
    bcftools index -t ${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.vcf.gz -Oz -o ${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz

    bcftools index -t ${meta.id}.${genome_version}.${readSet}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz

    """
}

process sawFish2_jointCall_all{
    label "high"
    conda "${params.sawfish2}"

    //publishDir "${outputDir}/jointCalls_All/", mode: 'copy', pattern: "*.sawfishSV_jointCall.*"

    input:
    val(x)
    
    output:
    tuple path("*.sawfishSV_jointCall.vcf.gz"),path("*.sawfishSV_jointCall.vcf.gz.tbi"),emit: sv_jointCall_vcf

    script:
    """
    sawfish joint-call \
    --threads ${task.cpus} \
    ${x} \
    --output-dir ${params.rundir}.sawfishSV_jointCall 
    
    mv ${params.rundir}.sawfishSV_jointCall/genotyped.sv.vcf.gz ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.vcf.gz

    mv ${params.rundir}.sawfishSV_jointCall/genotyped.sv.vcf.gz.tbi ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.vcf.gz.tbi
    """
}

process svdb_sawFish2_jointCall_all {
    label "low"
    conda "${params.svdb}"
    
    publishDir "${outputDir}/jointCalls_All/", mode: 'copy', pattern: "*_jointCall.svdb.*"


    when:
    params.jointCall

    input:
    tuple path(vcf), path(idx)
    
    output:
    path("*_jointCall.svdb.*")

    script:
    """
    svdb --query \
    --query_vcf ${vcf} \
    --sqdb ${sawfish_sqdb} > ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf
    bgzip ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf
    bcftools index -t ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf.gz -Oz -o ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    bcftools index -t ${params.rundir}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    """
}

process sawFish2_jointCall_caseID{
    tag "$caseID"
    label "high"
    conda "${params.sawfish2}"

    publishDir {params.groupedOutput ? "${outputDir}/${caseID}/jointCalls/" : "${outputDir}/jointCalls/"}, mode: 'copy', pattern: "*.sawfishSV.hiphase.svdb.*"

    //publishDir "${outputDir}/jointCalls_CaseID/", mode: 'copy', pattern: "*.sawfishSV_jointCall.*"

    input:
    tuple val(caseID), path(manifest)
    
    output:
    tuple val(caseID), path("*.sawfishSV_jointCall.vcf.gz"),path("*.sawfishSV_jointCall.vcf.gz.tbi"),emit: sv_jointCall_caseID_vcf

    script:
    """
    sawfish joint-call \
    --threads ${task.cpus} \
    --sample-csv ${manifest} \
    --output-dir ${caseID}.sawfishSV_jointCall 
    
    mv ${caseID}.sawfishSV_jointCall/genotyped.sv.vcf.gz ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.vcf.gz

    mv ${caseID}.sawfishSV_jointCall/genotyped.sv.vcf.gz.tbi ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.vcf.gz.tbi
    """
}

process svdb_sawFish2_jointCall_caseID {
    label "low"
    conda "${params.svdb}"
    
    publishDir {params.groupedOutput ? "${outputDir}/${caseID}/jointCalls/" : "${outputDir}/jointCalls/"}, mode: 'copy', pattern: "*_jointCall.svdb.*"


    when:
    params.jointCall

    input:
    tuple val(caseID), path(vcf), path(idx)
    
    output:
    path("*_jointCall.svdb.*")
    tuple val(caseID), path("${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz"),path("${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz.tbi"), emit: sawfish_caseID_AF10
    script:
    """
    svdb --query \
    --query_vcf ${vcf} \
    --sqdb ${sawfish_sqdb} > ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf
    bgzip ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf
    bcftools index -t ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.vcf.gz -Oz -o ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    bcftools index -t ${caseID}.${genome_version}.${readSet}.sawfishSV_jointCall.svdb.AF_below10pct.vcf.gz
    """
}

process svTopo {
    tag "$meta.id"
    label "high"
    conda "${params.svtopo}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/structuralVariants/SVtopo/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/structuralVariants/SVtopo/"}, mode: 'copy'


    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("${meta.id}.svtopo_out/")
    script:
    def exclude=params.genome=="hg38" ? "--exclude-regions ${cnv_exclude_sawfish}" : ""
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
    --genes ${gencode_gtf} \
    --image-type jpg 

    mv ${meta.id}.svtopo_out/index.html ${meta.id}.svtopo_out/${meta.id}.sawfishSV.svtopo.html
    """
}

process svTopo_filtered {
    tag "$meta.id"
    label "high"
    conda "${params.svtopo}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/structuralVariants/SVtopo_filtered/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/structuralVariants/SVtopo_filtered/"}, mode: 'copy'


    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("${meta.id}.svtopo_out/")
    script:
    def exclude=params.genome=="hg38" ? "--exclude-regions ${cnv_exclude_sawfish}" : ""
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
    --genes ${gencode_gtf} \
    --image-type jpg 

    mv ${meta.id}.svtopo_out/index.html ${meta.id}.svtopo_out/${meta.id}.sawfishSV.svtopo.html
    """
}


///////////////////////////////////////////////////
/////// ------- PSEUDO, VNTR, REPEATS, MITO ------- //
///////////////////////////////////////////////////


process mitorsaw {
    tag "$meta.id"
    label "medium"
    conda "${params.mitorsaw}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/specialAnalysis/mitochondrialVariants/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/specialAnalysis/mitochondrialVariants/"}, mode: 'copy'



    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.mitorsaw.*")

    script:
    """
    mitorsaw haplotype \
    --reference ${genome_fasta} \
    --bam ${data.bam} \
    --output-vcf ${meta.id}.${genome_version}.${readSet}.mitorsaw.vcf.gz \
    --output-hap-stats ${meta.id}.${genome_version}.${readSet}.mitorsaw.hapstats.json 

    """
}

process trgt4_diseaseSTRs{
   
    tag "$meta.id"
    label "low"
    conda "${params.trgt4}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/TRGT/bam" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/TRGT/bam"}, mode: 'copy', pattern: "*.sorted.ba*"

    publishDir "${lrsStorage}/STRs/repeatExpansions/TRGT/diseaseSTRs/", mode: 'copy', pattern:"*.sorted.vcf.*"

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.vcf.gz.tbi"),emit: str4_vcf
    
    tuple val(meta),path ("*.sorted.*")

    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.bam"), path("${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.bam.bai"), path("${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.vcf.gz.tbi"),emit: trgt_full
    
    script:
    def karyotype=(meta.sex=="male"||meta.sex=="M"||meta.genderFile=="M") ? "--karyotype XY" : "--karyotype XX"

    """
    trgt genotype \
    --genome ${genome_fasta} \
    --repeats ${tr_pathogenic_v2} \
    --reads ${data[0]} \
    $karyotype \
    --output-prefix ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive

    bcftools sort -Ov -o ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.vcf.gz ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.vcf.gz

    samtools sort -o ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.bam ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.spanning.bam
    samtools index ${meta.id}.${genome_version}.${readSet}.trgt4.STRchive.sorted.bam
    """
}

process trgt4_diseaseSTRs_plots{
    tag "$meta.id"
    label "low"
    conda "${params.trgt4}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/TRGT/diseaseSTRs/Plots/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/TRGT/Plots/"}, mode: 'copy', pattern: "*.{pdf,png,svg}"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    """
    trgt plot \
    --genome ${genome_fasta} \
    --repeats ${tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id ${data.strID} \
    --squished \
    -o ${meta.id}.${genome_version}.${readSet}.${data.strID}.allele.pdf

    trgt plot \
    --genome ${genome_fasta} \
    --repeats ${tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id ${data.strID} \
    --plot-type waterfall \
    -o ${meta.id}.${genome_version}.${readSet}.${data.strID}.waterfall.pdf

    """
}


process trgt4_all {

    tag "$meta.id"
    label "high"
    conda "${params.trgt4}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/TRGT/bam" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/TRGT/bam"}, mode: 'copy', pattern: "*.sorted.ba*"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/TRGT/allSTRs/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/TRGT/allSTRs/"}, mode: 'copy', pattern: "*.sorted.vcf.*"

    publishDir "${lrsStorage}/STRs/repeatExpansions/TRGT/all/", mode: 'copy', pattern:"*.sorted.vcf.*"

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.bam"), path("${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.bam.bai"),emit: str_spanning_bam
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.vcf.gz"), path("${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.vcf.gz.tbi"),emit: str4All_vcf
    
    script:
    def karyotype=(meta.sex=="male"||meta.sex=="M"||meta.genderFile=="M")  ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --genome ${genome_fasta} \
    --repeats ${tr_all} \
    --reads ${data[0]} \
    $karyotype \
    --output-prefix ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR

    bcftools sort -Ov -o ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.vcf.gz ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.vcf.gz 
    bcftools index -t ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.vcf.gz

    samtools sort -o ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.bam ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.spanning.bam
    samtools index ${meta.id}.${genome_version}.${readSet}.trgt4.allSTR.sorted.bam
    """
}

process kivvi_d4z4{
    tag "$meta.id"
    label "medium"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/Kivvi_D4Z4/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/Kivvi_D4Z4/"}, mode: 'copy'


    input:
    tuple val(meta), val(data)
   //  tuple val(meta), path(data)   
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.kivviD4Z4")
    
    script:
    """
    ${params.kivvi_dir}/kivvi \
    -r ${genome_fasta} \
    --bam ${data.bam} \
    -p ${meta.id}.${genome_version} \
    -o ${meta.id}.${genome_version}.${readSet}.kivviD4Z4 \
    d4z4
    """
}

process kivvi05_d4z4{
    tag "$meta.id"
    label "medium"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/repeatExpansions/kivviD4Z4_05/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/repeatExpansions/kivviD4Z4_05/"}, mode: 'copy'


    input:
    tuple val(meta), val(data)
   //  tuple val(meta), path(data)   
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.kivviD4Z4_05")
    
    script:
    """
    ${params.kivvi_dir2}/kivvi \
    -r ${genome_fasta} \
    --bam ${data.bam} \
    -p ${meta.id}.${genome_version} \
    -o ${meta.id}.${genome_version}.${readSet}.kivviD4Z4_05 \
    d4z4
    """
}


process paraphase {

    tag "$meta.id"
    label "medium"
    conda "${params.paraphaseMinimap2}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/specialAnalysis/paraphase/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/specialAnalysis/paraphase/"},mode: 'copy'



    input:
 //   tuple val(meta), path(aln)
    tuple val(meta), val(data)
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.paraphase/*")
        //-b ${aln[0]} \
    script:
    """
    paraphase \
    -b ${data.bam} \
    --reference ${genome_fasta} \
    -t ${task.cpus} \
    -o ${meta.id}.${genome_version}.${readSet}.hiphase.paraphase

    python ${localPythonScripts}/flatten_paraphaseSMN.py \
    --json ${meta.id}.${genome_version}.${readSet}.hiphase.paraphase/${meta.id}.paraphase.json \
    --out ${meta.id}.${genome_version}.${readSet}.hiphase.paraphase/${meta.id}.${genome_version}.${readSet}.paraphase.flattened.tsv

    python ${localPythonScripts}/flatten_paraphaseSMN.py \
    --json ${meta.id}.${genome_version}.${readSet}.hiphase.paraphase/${meta.id}.paraphase.json \
    --loci SMN1,PMS2,IKBKG \
    --out ${meta.id}.${genome_version}.${readSet}.hiphase.paraphase/${meta.id}.${genome_version}.${readSet}.paraphase.flattened.SMN_PMS2_IKBKG.tsv
    """
}


process starphase {

    tag "$meta.id"
    label "medium"
    conda "${params.starphase}"
    
    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/specialAnalysis/starphase/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/specialAnalysis/starphase/"},mode: 'copy'


    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.starphase.*")
    
    script:
    """
    pbstarphase diplotype \
    --database ${starphase_db} \
    --bam ${data.bam} \
    --reference ${genome_fasta} \
    --vcf ${data.dv_vcf} \
    --sv-vcf ${data.sawfish_vcf} \
    --pharmcat-tsv ${meta.id}.${genome_version}.${readSet}.starphase.pharmcat.tsv \
    --output-calls ${meta.id}.${genome_version}.${readSet}.starphase.json

    """
}


process advntr {

    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 8
    publishDir "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/advntr/", mode: 'copy', pattern: "*.advntr.*"

    conda "${params.advntr15}"

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.advntr.*")
    
    script:
    """
    advntr genotype \
    -f ${data[0]} \
    --pacbio \
    -m ${vntr_defaultModel} \
    -o ${meta.id}.advntrDefault
    """
}

/////////////////////// TRIO: Variant prioritization with Exomiser 14.1.0 ///////////////////////

process exo14_2508_exome {
    label "medium"
    tag "$caseID"

    publishDir "${outputDir}/${caseID}/exomiser14_2508/exomiser/", mode: 'copy'
    publishDir "${outputDir}/${caseID}/documents/", mode: 'copy',pattern:"*.{hpo.txt,yml,ped}"
    input:
    tuple val(caseID), path(vcf), path(idx), path(hpo), path(samplesheet)

    output:
    path("*.{html,tsv,vcf,json,hpo.txt,yml,ped}")

    script:
    """
    python3 ${localPythonScripts}/make_ped_and_family.py \
    --samplesheet ${samplesheet} \
    --vcf ${vcf} \
    --hpo ${hpo}

    java -jar ${localProgramPath}/exomiser-cli-14.1.0/exomiser-cli-14.1.0.jar \
    --sample ${caseID}-family.yml \
    --analysis ${exome_yml} \
    --spring.config.location=${localProgramPath}/exomiser-cli-14.1.0/
    
    mv results/* .
    mv exomiser_tmp.html ${caseID}.exo14_2508.html
    mv exomiser_tmp.variants.tsv ${caseID}.exo14_2508.Variants.tsv
    mv exomiser_tmp.genes.tsv ${caseID}.exo14_2508.Genes.tsv
    mv exomiser_tmp.json ${caseID}.exo14_2508.json
    """

}

process exo14_2508_genome {
    label "medium"
    tag "$caseID"

    publishDir "${outputDir}/${caseID}/exomiser14_2508/genomiser/", mode: 'copy'

    input:
    tuple val(caseID), path(vcf), path(idx), path(hpo), path(samplesheet)

    output:
    path("*.{html,tsv,vcf,json}")

    script:
    """
    python3 ${localPythonScripts}/make_ped_and_family.py \
    --samplesheet ${samplesheet} \
    --vcf ${vcf} \
    --hpo ${hpo}
    
    java -jar ${localProgramPath}/exomiser-cli-14.1.0/exomiser-cli-14.1.0.jar \
    --sample ${caseID}-family.yml \
    --analysis ${genome_yml} \
    --spring.config.location=${localProgramPath}/exomiser-cli-14.1.0/
    
    mv results/* .
    mv exomiser_tmp.html ${caseID}.genomiser14_2508.html
    mv exomiser_tmp.variants.tsv ${caseID}.genomiser14_2508.Variants.tsv
    mv exomiser_tmp.genes.tsv ${caseID}.genomiser14_2508.Genes.tsv
    mv exomiser_tmp.json ${caseID}.genomiser14_2508.json
    """
}

process exo14_2508_SV {
    label "medium"
    tag "$caseID"

    publishDir "${outputDir}/${caseID}/exomiser14_2508/exomiserStructuralVariants/", mode: 'copy'

    input:
    tuple val(caseID), path(vcf), path(idx), path(hpo), path(samplesheet)

    output:
    path("*.{html,tsv,vcf,json,hpo.txt,yml,ped}")

    script:
    """
    zcat ${vcf} | sed 's/^##fileformat=VCFv4\\.4/##fileformat=VCFv4.2/'| bgzip > ${caseID}.sawfish.forExomiser.vcf.gz

    python3 ${localPythonScripts}/make_ped_and_family.py \
    --samplesheet ${samplesheet} \
    --vcf ${caseID}.sawfish.forExomiser.vcf.gz \
    --hpo ${hpo}

    java -jar ${localProgramPath}/exomiser-cli-14.1.0/exomiser-cli-14.1.0.jar \
    --sample ${caseID}-family.yml \
    --analysis ${exome_yml} \
    --spring.config.location=${localProgramPath}/exomiser-cli-14.1.0/
    
    mv results/* .
    mv exomiser_tmp.html ${caseID}.SVs.exo14_2508.html
    mv exomiser_tmp.variants.tsv ${caseID}.SVs.exo14_2508.Variants.tsv
    mv exomiser_tmp.genes.tsv ${caseID}.SVs.exo14_2508.Genes.tsv
    mv exomiser_tmp.json ${caseID}.SVs.exo14_2508.json
    """

}


///////////////////////////////////////////////////
/////// ------- METHYLATION ------- ///////////////
///////////////////////////////////////////////////
process methylationBW{
    
    tag "$meta.id"
    label "medium"
    conda "${params.pbCPGtools}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/specialAnalysis/methylation/BigWigBed/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/specialAnalysis/methylation/BigWigBed/"}, mode: 'copy', pattern: "*.methylation.{hap1,hap2,combined}.*"

    publishDir "${lrsStorage}/methylation/pbCpGtools/${meta.id}/", mode: 'copy', pattern:"*.bed.*"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.hiphase.methylation*")
    
    script:
    """
    aligned_bam_to_cpg_scores \
    --bam ${data.bam} \
    --output-prefix ${meta.id}.${genome_version}.${readSet}.hiphase.methylation
    """
}

process methylationSegm{
    tag "$meta.id"
    label "medium"
    conda "${params.methbat}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/specialAnalysis/methylation/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/specialAnalysis/methylation/"}, mode: 'copy'

    publishDir "${lrsStorage}/methylation/methBatProfiles/", mode: 'copy', pattern:"*.profile"


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.met.*")
    
    script:
    """
    methbat segment \
    --input-prefix ${meta.id}.${genome_version}.${readSet}.hiphase.methylation \
    --output-prefix ${meta.id}.met.segments

    methbat profile \
    --input-prefix ${meta.id}.${genome_version}.${readSet}.hiphase.methylation \
    --input-regions ${methylationBackground} \
    --output-region-profile ${meta.id}.met.profile

    methbat profile \
    --input-prefix ${meta.id}.${genome_version}.${readSet}.hiphase.methylation \
    --input-regions ${methylationBackgroundLocal} \
    --output-region-profile ${meta.id}.met.profileLOCAL
    """
}


///////////////////////////////////////////////////
/////// ------- QUALITY CONTROL ------- ///////////
///////////////////////////////////////////////////

process mosdepthROI {
    tag "$meta.id"
    label "low"
    conda "${params.mosdepth}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/QC/mosdepth/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/mosdepth/"}, mode: 'copy'

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.${genome_version}_roi.*"),emit: mosdepth_roi
    tuple val(meta), path("*.region.dist.txt"), emit:multiqc
    script:
    def callable=params.genome=="hg38" ? "--by ${CALLABLE_ROI}" : "--by 1000"
    """
    mosdepth \
    -t ${task.cpus} \
    $callable \
    ${meta.id}.${genome_version}_roi \
    ${data[0]}

    """
}


process whatsHap_stats {
    tag "$meta.id"
    label "low"
    conda "${params.whatshap}" 

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/QC/whatsHap/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/whatsHap/"}, mode: 'copy'

    input: 
    tuple val(meta), val(data)  
    output:
    tuple val(meta), path("${meta.id}.whatshap.stats.tsv"),emit:multiqc

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
    conda "${params.cramino}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/QC/cramino/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/cramino/"}, mode: 'copy'

    input: 
    tuple val(meta), val(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.craminoQC.txt")

    script:
    """
    cramino \
    -t ${task.cpus} \
    --karyotype \
    --phased \
    ${data.bam} > ${meta.id}.${genome_version}.${readSet}.craminoQC.txt
    """
}

process nanoStat {
    tag "$meta.id"
    label "low"
    conda "${params.nanostats}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/QC/nanoStat/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/nanoStat/"}, mode: 'copy'

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.${readSet}.nanostat.txt"),emit: multiqc
    path("${meta.id}.${genome_version}.${readSet}.nanostat.txt")
    script:
    """
    NanoStat \
    -t ${task.cpus} \
    -n ${meta.id}.${genome_version}.${readSet}.nanostat.txt \
    --bam ${data[0]}
    """
}

process multiQC {
    tag "$meta.id"
    label "low"
    conda "${params.multiqc}"

    publishDir {params.groupedOutput ? "${outputDir}/${meta.caseID}/QC/" : "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/"}, mode: 'copy'

    when:
    !params.groupedOutput

    input:
    tuple val(meta),  path(data)  

    output:
    path ("*MultiQC*.html")

    script:
    """
    multiqc \
    -c ${multiqc_config} \
    -f -q ${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/ \
    -n ${meta.id}.MultiQC.DNA.html
    """
}
//    -f -q ${launchDir}/${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/QC/ \


process multiQC_ALL {
    label "low"
    conda "${params.multiqc}"

    publishDir "${outputDir}/", mode: 'copy'

    when:
    params.groupedOutput

    input:
    tuple val(meta),path(data)  

    output:
    path ("${params.rundir}.MultiQC.ALL.html")

 //   def exclude=params.genome=="hg38" ? "--cnv-excluded-regions ${cnv_exclude_sawfish}" : ""
   // def sex = (meta.sex=="male"||meta.sex=="M"||meta.genderFile=="M") ? "--expected-cn ${sawfishMaleExpectedCN}" : "--expected-cn ${sawfishFemaleExpectedCN}"


    script:
    def qcdir = params.groupedOutput ? "${launchDir}/${outputDir}/*/QC/" : "${launchDir}/${outputDir}/*/*/QC/"
    
    """
    multiqc \
    -c ${multiqc_config} \
    -f -q $qcdir \
    -n ${params.rundir}.MultiQC.ALL.html
    """
}
//    -f -q ${launchDir}/${outputDir}/*/*/QC/ \
/////////////// TO DO /////////////////////
/*
process collect_versions {
    input:
    val versions from versions_ch.collect()

    output:
    path 'software_versions.tsv'

    script:
    """
    echo -e "process\\ttool\\tversion" > software_versions.tsv
    for row in "${versions[@]}"; do
        echo "$row" | tr -d '[]' | tr ',' '\\t' >> software_versions.tsv
    done
    """
}
*/
















process vntyper2 {
    errorStrategy 'ignore'
    publishDir "${outputDir}/MUC1-VNTR_kestrel/", mode: 'copy'
    cpus 16

    input:
    tuple val(meta), path(reads)

    output:
    //tuple val(meta), path("vntyper${meta.id}.vntyper/*")
    tuple val(meta), path("*/*.{tsv,vcf}")
    script:
    
    def reads_command = "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
    
    """
    singularity run -B ${s_bind} ${simgpath}/vntyper20.sif \
    -ref ${vntyperREF}/chr1.fa \
    --fastq1 ${r1} --fastq2 ${r2} \
    -t ${task.cpus} \
    -w vntyper \
    -m ${vntyperREF}/hg19_genic_VNTRs.db \
    -o ${meta.id} \
    -ref_VNTR ${vntyperREF}/MUC1-VNTR_NEW.fa \
    --fastq \
    --ignore_advntr \
    -p ${localProgramPath}/vntyper/VNtyper/
    """
}


process build_symlinks {
    tag "build_symlinks"

    // No input channels are required; Nextflow will wait until all upstream processes are done
    // before scheduling this process, if we make it the last step in the workflow.

    output:
    path "symlinks.done"

    script:
    """
    @{localBashScripts}/createSymlinks_after_nextflow.sh ${testLinksInput} ${testLinksOutput}
    touch symlinks.done
    """
}


////////////// END TO DO END ///////////////////


///////////////////////////////////////////////////
////// ------- DE NOVO ASSEMBLY ------- ///////////
///////////////////////////////////////////////////

process hifiasm {
    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 8
    publishDir "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/deNovoAssembly/hifiasm/", mode: 'copy', pattern: "*.advntr.*"

    conda "${params.hifiasm}"

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.advntr.*")
    
    script:
    """
    samtools fastq \
    -@ 12 \
    ${data[0]} > ${meta.id}.hifireads.fastq
    hifiasm \
    -t ${task.cpus} \
    -o ${meta.id}.asm \
    ${meta.id}.hifireads.fastq

    """
}

/*
process asm_to_fasta {
    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 4
    publishDir "${outputDir}/${meta.caseID}/${meta.rekv}_${meta.id}_${meta.groupKey}_${readSet}/deNovoAssembly/hifiasm/", mode: 'copy', pattern: "*.fasta"



    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}.asm.p_ctg.fasta"), emit: asm_fasta

    script:
    """
    awk '/^>/ {print ">"$0; next} {print}' ${data[0]} | sed 's/ /_/g' > ${meta.id}.asm.p_ctg.fasta
    """

}

*/




















