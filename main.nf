#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
import java.util.Locale
date=new Date().format( 'yyMMdd' )
date2=new Date().format( 'yyMMdd HH:mm:ss' )
user="$USER"
runID="${date}.${user}"

log.info """\
======================================================
Clinical Genetics Vejle: PacBio LRS v4
======================================================
Genome        : $params.genome
GenomeDir     : $refFilesDir
Input Readset : $inputReadSet_allDefault
read Subset   : $readSubset_hifiDefault
RunID         : $runID
Script start  : $date2
Genome FASTA  : ${genome_fasta}
Archive RAW   : ${dataArchive}
OutputDirBase : ${outputDirBase}
workDir       : ${workflow.workDir}
layout        : $params.layoutMode
min input GB  : $params.minGB
"""


/* ----- Changes:

    - deprecate obsolete versions of programs:
    - TRGT4, paraphase3, pb-cpgtools, kivvi05

    - Current versions:
        TRGT5, paraphase4, methBat profile, kivvi v1

    - Remove allReads bam/cram output


*/

//////////// DEFAULT INPUT ///////////////////////

def inputError() {
    log.info"""
    USER INPUT ERROR: The user should point to a samplesheet (--samplesheet parameter) or input folder containing all data to be used as input (--input parameter).
    """.stripIndent()
}

def hpoInputError() {
    log.info"""
    USER INPUT ERROR: A samplesheet (--samplesheet parameter) containing 5 columns (caseID, samplename, gender, relation and affection status) is required when usign --hpo.  
    """.stripIndent()
}


if (!params.samplesheet && !params.input && !params.familySS) exit 0, inputError() 
if (!params.samplesheet && params.hpo && !params.familySS) exit 0, hpoInputError() 


if (params.hpo) {
    channel.fromPath(params.hpo)
    |set { hpo_ch }
}


if (params.aligned) {

    inputBam="${params.input}/*.bam"
    inputBai="${params.input}/*.bai"

    channel.fromPath(inputBam, followLinks: true)
    |map { f -> tuple(f.baseName, f) }
    |map {id,bam -> 
            (samplename,genomeversion)      =id.tokenize(".")
            meta=[id:samplename,genomeversion:genomeversion,type:"aligned"]
            tuple(meta,bam)        
        }
    |set {bamInput}

    channel.fromPath(inputBai, followLinks: true)
    |map { f -> tuple(f.baseName, f) }
    |map {id,bai -> 
            (samplename,genomeversion)      =id.tokenize(".")
            meta=[id:samplename,genomeversion:genomeversion,type:"aligned"]
            tuple(meta,bai)        
        }
    |set {baiInput}  

    bamInput.join(baiInput)
    |map { meta,bam,bai -> tuple(meta,[bam,bai]) } 
    | set {alignedInput_tmp}


    if (params.samplesheet) {
        alignedInput_tmp.join(samplesheet_join)
        |map {metaData,metaSS,meta,bam -> tuple(metaData,bam)}
        |set {alignedFinal}
    }
    if (!params.samplesheet) {
        alignedInput_tmp
        |set {alignedFinal}
    }
}


if (!params.aligned) {

    if (params.input) {
        if (params.hifiReads){
            inputBam="${params.input}/**/*.hifi_reads.*.bam"
        }
        if (params.failedReads){
            inputBam="${params.input}/**/*.fail_reads.*.bam"
        }
        if (!params.hifiReads && !params.failedReads) {
            inputBam="${params.input}/**/*.bam"
        }
    }
    
    if (!params.input) {
        if (params.hifiReads){
            inputBam="${params.dataArchive}/**/hifi_reads/*.hifi_reads.*.bam"
        }
        if (params.failedReads){
            inputBam="${params.dataArchive}/**/failed_reads/*.fail_reads.*.bam"
        }
        if (!params.hifiReads && !params.failedReads) {
            inputBam="${params.dataArchive}/**/*_reads/*.bam"
        }
    }

    if (params.samplesheet && !params.customSS && !params.jointSS) {
              
        def ssBase = params.samplesheet
                    .toString()
                    .tokenize('/')
                    .last()
                    .replaceFirst(/_metadata$/, '')

        channel.fromPath(params.samplesheet)
        | splitCsv(sep:'\t')
        |map { row ->
            (rekv, npn,material,testlist,gender,proband,intRef) = row[0].tokenize("_")
            def groupKey    = (intRef == 'noInfo')  ? "single" : intRef
            def sex         = (gender =="K")        ? "female" : "male"
            meta=[  id          :npn,
                    testlist    :testlist,
                    sex         :sex,
                    proband     :proband,
                    intRef      :intRef,
                    rekv        :rekv,
                    groupKey    :groupKey,
                    ssBase      :ssBase]
            meta
            }
        | set {samplesheet_full}
        samplesheet_full
        |branch {row ->
            singleSample: (row.groupKey== "single")
                return row
            multiSample: true
                return row
        }
        |set {samplesheetBranch}
    }
 
   /*
    if (params.samplesheet && params.customSS) {

        def ssBase = params.samplesheet
                    .toString()
                    .tokenize('/')
                    .last()
                    .replaceFirst(/\.txt$/, '')

        channel.fromPath(params.samplesheet)
        | splitCsv(sep:'\t')
        |map { row -> 
            (caseID, samplename, sex,outKey) =tuple(row)
            meta=[caseID:caseID,id:samplename,sex:sex,groupKey:"customSampleSheet",outKey:caseID,ssBase:ssBase,rekv:outKey] // edit back to normal, if needed.
            meta
        }
        | set {samplesheet_full}
    }
*/
    if (params.samplesheet && (params.jointSS || params.familySS|| params.customSS) ) {
        def ssBase = params.samplesheet
                    .toString()
                    .tokenize('/')
                    .last()
                    .replaceFirst(/\.txt$/, '')

    channel.fromPath(params.samplesheet)
    .splitCsv(sep: '\t')
    .map { row ->
        def (rekv, npn, material, testlist, gender, proband, intRef) = row
        def sex = (gender == 'K') ? 'female' : 'male'

        def meta = [
        rekv     : rekv,
        id       : npn,
        material : material,
        testlist : testlist,
        gender   : gender,
        sex      : sex,
        proband  : proband,
        intRef   : intRef,
        ssBase   : ssBase,
        groupKey : intRef
        ]

        tuple(intRef, meta)
    }
    .groupTuple()
    .flatMap { intRef, metas ->

        def probands = metas.findAll { it.proband == 'T' }
        assert probands && probands.size() >= 1 : "No proband (T) found for intRef=${intRef}"

        def anchor = probands[0]
        def caseID = "${anchor.rekv}_${anchor.testlist}_${intRef}"

        metas.collect { m ->
            def relation
            if( m.proband == 'T' ) {
                relation = 'index'
            } else if( m.gender == 'M' ) {
                relation = 'pater'
            } else if( m.gender == 'K' ) {
                relation = 'mater'
            } else {
                relation = 'unknown_relation'
            }

            // Return NEW map (don’t mutate original)
            m + [
                caseID: caseID,
                relation: relation
            ]
        }
    }
    | set {samplesheet_full}
    }


    if (params.samplesheet) {
        channel.fromPath(inputBam, followLinks: true)
        |map { f -> tuple(f.baseName, f) }
        |map {id,bam -> 
                (samplenameFull,pacbioID,readset,barcode)   =id.tokenize(".")
                (instrument,date,time)                      =pacbioID.tokenize("_")     
                (samplename,material,testlist,gender)       =samplenameFull.tokenize("_")
                //meta=[id:samplename,genderFile:gender,testlistFile:testlist]
                meta=[id:samplename]
                tuple(meta,bam)   
            }
        |groupTuple(sort:true)
        | map { meta, bams ->
            long totalBytes = (bams.sum { it.size() } as long)
            double totalGB  = totalBytes / (1024.0 * 1024 * 1024)
            def meta2 = meta + [
                nBams       : bams.size(),
                totalsizeGB : totalGB
            ]
            tuple(meta2, bams)
        }
        |branch  {meta,bam -> 
            UNASSIGNED: (meta.id=~/UNASSIGNED/)
                return [meta,bam]
            samples: true
                return [meta,bam]
        }
        | set { ubam_input }

        ubam_input.samples
            | map { meta, bam -> tuple(meta.id,meta,bam) }
        |set {ubam_input_samples}    

        if (!params.singleOnly && !params.intrefOnly) {
            samplesheet_full
            |map {row -> meta2=[row.id,row]}
            |set {samplesheet_join}
        }
        if (params.singleOnly) {
            samplesheetBranch.singleSample
            |map {row -> meta2=[row.id,row]}
            |set {samplesheet_join}
        }
        if (params.intrefOnly) {
            samplesheetBranch.multiSample
            |map {row -> meta2=[row.id,row]}
            |set {samplesheet_join}
        }
        samplesheet_join.join(ubam_input_samples)
            |map {samplename, metaSS, metaData, bam -> tuple(metaSS+metaData,bam)}
        |set {ubam_ss_merged} // full unfiltered set

        //write info of full set to summary file:

        ubam_ss_merged
        .map { meta, bams ->
            def gb = String.format(Locale.US, "%.2f", (meta.totalsizeGB as double))
            "${meta.id}\t${meta.nBams}\t${inputReadSet_allDefault}\t${gb}\t${meta.testlist}"
        }
        .collect()
        | map { lines ->
            def header  ="sample\tbamcount\treadSet\ttotal_gb\ttestlist"
            ([header] + lines).join("\n")
        }
        |set {ubam_size_summary_ch}

        //Branch by total input size (i.e. drop all samples with combined ubam size < e.g. 30GB)
        ubam_ss_merged
            |branch { meta, bams ->
                keep:   (meta.totalsizeGB as double) >= params.minGB 
                    return [meta, bams]
                drop:   true
                    return [meta, bams]
            }
        |set { ubam_ss_merged_size_split }

        //write out dropped samples info
        ubam_ss_merged_size_split.drop
        .map { meta, bams ->
            def gb = String.format(Locale.US, "%.2f", (meta.totalsizeGB as double))
            "${meta.id}\t${meta.nBams}\t${inputReadSet_allDefault}\t${gb}\t${meta.testlist}"
        }
        .collect()
        | map { lines ->
            def header  ="sample\tbamcount\treadSet\ttotal_gb\ttestlist"
            ([header] + lines).join("\n")
        }
        |set {ubam_size_dropped_ch}

        ubam_ss_merged_size_split.keep 
        .map { meta, bams ->
            def gb = String.format(Locale.US, "%.2f", (meta.totalsizeGB as double))
            "${meta.id}\t${meta.nBams}\t${inputReadSet_allDefault}\t${gb}\t${meta.testlist}"
        }
        .collect()
        | map { lines ->
            def header  ="sample\tbamcount\treadSet\ttotal_gb\ttestlist"
            ([header] + lines).join("\n")
        }
        |set {ubam_size_keep_ch}

        ubam_ss_merged_size_split.keep      // All data passing size limit - ready for downstream
        |set {finalUbamInput}
        
        channel.fromPath(params.samplesheet)
        |set {samplesheet_path_ch}
    }

    if (!params.samplesheet) {
        channel.fromPath(inputBam, followLinks: true)
        |map { f -> tuple(f.baseName, f) }

        |map {id,bam -> 
                (samplenameFull,pacbioID,readset,barcode)   =id.tokenize(".")
                (instrument,date2,time)                      =pacbioID.tokenize("_")     
                (samplename,material,testlist,gender)       =samplenameFull.tokenize("_")
                meta=[id:samplename,caseID:date+"_"+testlist, gender:gender,rundate:date,testlist:testlist]
                tuple(meta,bam)        
            }

        |groupTuple(sort:true)
        | map { meta, bams ->
            long totalBytes = (bams.sum { it.size() } as long)
            double totalGB  = totalBytes / (1024.0 * 1024 * 1024)
            def meta2 = meta + [
                nBams       : bams.size(),
                totalsizeGB : totalGB
            ]
            tuple(meta2, bams)
        }
        |branch  {meta,bam -> 
            UNASSIGNED: (meta.id=~/UNASSIGNED/)
                        return [meta,bam]
            samples: true
                        return [meta,bam]
        }
        | set {ubam_input }
        
        ubam_input.samples
        |set {finalUbamInput}
    }

}



/////////////////// MODULES ///////////////////////


 include {
        write_input_summary;
        write_dropped_samples_summary;
        symlinks_ubam_dropped;
        write_analyzed_samples_summary;
        hiPhase;
        } from "./modules/dnaModules.nf" 
///////////////// SUBWORKFLOWS ///////////////////////

include { PREPROCESS }              from './subworkflows/PREPROCESS.nf'
include { PRE_PHASING }             from './subworkflows/PRE_PHASING.nf'
include { POST_PHASING }            from './subworkflows/POST_PHASING.nf'
include { FAMILY_ANALYSIS }         from './subworkflows/FAMILY_ANALYSIS.nf'
include { FAMILY_ANALYSIS_ENTRY }   from './subworkflows/FAMILY_ANALYSIS.nf'
////////////////// WORKFLOWS AND PROCESSES ///////////////////////


workflow {

    if (!params.aligned) {
        write_input_summary(ubam_size_summary_ch)
        write_analyzed_samples_summary(ubam_size_keep_ch)
        write_dropped_samples_summary(ubam_size_dropped_ch)
        symlinks_ubam_dropped(ubam_ss_merged_size_split.drop)
        PREPROCESS(finalUbamInput)
    }
    PRE_PHASING(PREPROCESS.out.alignedFinal)

   // hiPhase(PRE_PHASING.out.hiphaseInput)
    hiPhaseTwoAln(PRE_PHASING.out.hiphaseInput)
    
    hiPhaseTwoAln.out.hifi_bam

        .join(hiPhase.out.dv_vcf)
        .join(hiPhase.out.sawfish_vcf)
        .join(PRE_PHASING.out.sawfish_supporting_reads)
        | map { meta, bam, bai, dv_vcf, dv_idx, sv_vcf, sv_idx, sv_jsonReads ->
            tuple(meta, [
                bam:           bam,
                bai:           bai,
                dv_vcf:        dv_vcf,
                dv_idx:        dv_idx,
                sawfish_vcf:   sv_vcf,
                sawfish_idx:   sv_idx,
                sawfish_reads: sv_jsonReads
            ])
        }
    | set { phasedAll }  // use for val(data) instead of path(data) setup in modules 

    POST_PHASING(
                phasedAll,
                PRE_PHASING.out.sawfish_supporting_reads,
                PRE_PHASING.out.mosdepth,
                PRE_PHASING.out.nanoStat
                )

    def hpo_ch = params.hpo        
        ? channel.fromPath(params.hpo)
        : channel.empty()

    def ss_ch  = params.samplesheet 
        ? channel.fromPath(params.samplesheet) 
        : channel.empty()


    if (params.jointCall || params.jointSS) {

        PRE_PHASING.out.sawfish_discover_dir   // tuple(meta), path(dir), val(bam)
        | map { meta, dir, bam ->
            tuple(
            meta.caseID, tuple(meta, "${dir.toString()}, ${bam.toString()}")
            )
        }
        | groupTuple()   // -> caseID, [ (meta,line), (meta,line), ... ]
        | map { caseID, records ->

            def anchorMeta = records[0][0]

            // build manifest file content
            def content = records.collect { it[1] }.join("\n") + "\n"

            // write the manifest to a file in the work dir
            def mf = file("${caseID}.sawFishJoinCall.manifest.csv")
            mf.text = content

            // emit tuple(meta, manifest)
            tuple(anchorMeta, mf)
        }
        | set { sawfish_jointCall_manifest_ch }

        PRE_PHASING.out.dv_gvcf
        .map { meta, gvcf, tbi ->
            // store one record per sample: (caseID, meta, gvcfPath)
            tuple(meta.caseID, tuple(meta, gvcf.toString()))
        }
        .groupTuple()
        .map { caseID, records ->
            def anchorMeta = records[0][0]
            def content = records.collect { it[1] }.join('\n') + '\n'
            def mf = file("${caseID}.manifest")
            mf.text = content
            tuple(anchorMeta, mf)
        }
        .set { glnexus_manifest_ch }

        FAMILY_ANALYSIS(
                        glnexus_manifest_ch,
                        sawfish_jointCall_manifest_ch,
                        hpo_ch,
                        ss_ch
                        )

    }
}

// Virker ikke lige pt.:
workflow.onComplete {

    if( !params.createSymlinks ) {
        log.info "Symlink maintenance disabled by config."
        return
    }

    if( !workflow.success ) {
        log.warn "Workflow failed – skipping symlink maintenance."
        return
    }

    def mirrorScript  = params.mirrorSampleData
    def collectScript = params.collectDataTypeSymlink

    if( !mirrorScript || !collectScript ) {
        log.warn "Symlink script paths not defined in config – skipping."
        return
    }

    def cmds = [
        "bash '${collectScript}'",
        "bash '${mirrorScript}'"
        
    ]

    cmds.each { cmd ->
        log.info "onComplete: running: ${cmd}"

        try {
            def p = ["bash", "-lc", cmd].execute()
            p.waitForProcessOutput(System.out, System.err)

            if( p.exitValue() != 0 ) {
                log.warn "onComplete: command failed (exit ${p.exitValue()}): ${cmd}"
            } else {
                log.info "onComplete: finished OK: ${cmd}"
            }
        }
        catch(Exception e) {
            log.warn "onComplete: exception while running '${cmd}': ${e.message}"
        }
    }
}

