#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
import java.util.Locale
date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"

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



if (!params.samplesheet && !params.input) exit 0, inputError() 
if (!params.samplesheet && params.hpo) exit 0, hpoInputError() 




if (params.hpo) {
    channel.fromPath(params.hpo)
    |set { hpo_ch }
}






if (params.aligned) {

    inputBam="${params.input}/*.bam"
    inputBai="${params.input}/*.bai"

    Channel.fromPath(inputBam, followLinks: true)
    |map { tuple(it.baseName,it) }
    |map {id,bam -> 
            (samplename,genomeversion)      =id.tokenize(".")
            meta=[id:samplename,genomeversion:genomeversion,type:"aligned"]
            tuple(meta,bam)        
        }
    |set {bamInput}

    Channel.fromPath(inputBai, followLinks: true)
    |map { tuple(it.baseName,it) }
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
            inputBam="${params.dataArchive}/**/*.bam"
        }
    }
    
    if (!params.input) {
        if (params.hifiReads){
            inputBam="${params.dataArchive}/**/*.hifi_reads.*.bam"
        }
        if (params.failedReads){
            inputBam="${params.dataArchive}/**/*.fail_reads.*.bam"
        }
        if (!params.hifiReads && !params.failedReads) {
            inputBam="${params.dataArchive}/**/*.bam"
        }
    }

    if (params.samplesheet && !params.intSS && !params.jointSS) {
              
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
            def outKey      = (intRef == 'noInfo')  ? "singleSampleAnalysis" : "multiSampleAnalysis"
            def sex         = (gender =="K")        ? "female" : "male"
            meta=[  id          :npn,
                    testlist    :testlist,
                    sex         :sex,
                    proband     :proband,
                    intRef      :intRef,
                    rekv        :rekv,
                    groupKey    :groupKey,
                    outKey      :outKey,
                    ssBase      :ssBase]
            meta
            }

        | set {samplesheet_full}
        samplesheet_full
        |branch {row ->
            singleSample: (row.groupKey=~/single/)
                return row
            multiSample: true
                return row
        }
        |set {samplesheetBranch}
    }

    // intermediate naming scheme:
    if (params.samplesheet && params.intSS) {

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

    if (params.samplesheet && params.jointSS) {
        def ssBase = params.samplesheet
                    .toString()
                    .tokenize('/')
                    .last()
                    .replaceFirst(/\.txt$/, '')

    Channel
    .fromPath(params.samplesheet)
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
        outKey   : 'multiSampleAnalysis',
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
        Channel.fromPath(inputBam, followLinks: true)
        |map { tuple(it.baseName,it) }
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

        if (!params.singleOnly) {
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
        Channel.fromPath(inputBam, followLinks: true)
        |map { tuple(it.baseName,it) }

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
        |view
        |set {finalUbamInput}
    }

}



/////////////////// MODULES ///////////////////////
include {pbmm2_align;
        create_fofn;
        pbmm2_align_mergedData;
        extractHifi;
        inputFiles_symlinks_ubam;
        sawFish2;
        svdb_SawFish;
        sawFish2_jointCall_all;
        svdb_sawFish2_jointCall_all;
        sawFish2_jointCall_caseID;
        svdb_sawFish2_jointCall_caseID;
        deepvariant;
        glNexus_jointCall;
        trgt4_diseaseSTRs;
        trgt4_diseaseSTRs_plots;
        trgt4_all;
        kivvi_d4z4;
        methylationBW;
        paraphase;
        starphase;
        methylationSegm;
        multiQC;
        multiQC_ALL;
        mosdepthROI;
        cramino;
        nanoStat;
        whatsHap_stats;
        hiPhase;
        build_symlinks;
        check_tmpdir;
        svTopo;
        svTopo_filtered;
        mitorsaw;
        exo14_2508_exome;
        exo14_2508_genome;
        exo14_2508_SV;
        kivvi05_d4z4;
        write_input_summary;
        write_dropped_samples_summary;
        symlinks_ubam_dropped;
        write_analyzed_samples_summary;
        } from "./modules/dnaModules.nf" 


puretargetPlotGenes=[
        "DRPLA_ATN1",
        "SCA1_ATXN1",
        "SCA10_ATXN10",
        "SCA2_ATXN2",
        "SCA3_ATXN3",
        "SCA7_ATXN7",
        "SCA8_ATXN8OS",
        "SCA31_BEAN1",
        "SCA6_CACNA1A",
        "SCA37_DAB1",
        "SCA27B_FGF14",
        "SCA36_NOP56",
        "SCA12_PPP2R2B",
        "SCA17_TBP",
        "SCA4_ZFHX3",
        "FXS_FMR1",
        "FRAXE_AFF2",
        "FTDALS1_C9orf72",
        "FRDA_FXN",
        "CANVAS_RFC1",
        "NIID_NOTCH2NLC",
        "DM1_DMPK",
        "DM2_CNBP",
        "HD_HTT",
        "HDL2_JPH3",
        "SBMA_AR",
        "OPMD_PABPN1",
        "OPDM5_ABCD3",
        "OPDM2_GIPC1",
        "OPDM1_LRP12",
        "OPDM4_RILPL1",
        "CCHS_PHOX2B",
        "CJD_PRNP",
        "FAME1_SAMD12",
        "EPM1_CSTB",]


////////////////// WORKFLOWS AND PROCESSES ///////////////////////

workflow PREPROCESS {

    take:
    finalUbamInput     
   
    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    create_fofn(finalUbamInput)
    pbmm2_align_mergedData(create_fofn.out)

    emit:
    alignedAll=pbmm2_align_mergedData.out.bamAll


    
}


workflow VARIANTS {

    take:
    aligned    
    main:
    deepvariant(aligned)

    emit:
    dv_vcf=deepvariant.out.dv_vcf
    dv_gvcf=deepvariant.out.dv_gvcf
}

workflow STRUCTURALVARIANTS {

    take:
    aligned

    main:

    sawFish2(aligned)

    emit:
    sawfish_vcf=sawFish2.out.sv_vcf
    sawfish_discover_dir=sawFish2.out.sv_discover_dir
    sawfish_discover_dir2=sawFish2.out.sv_discover_dir2
    sawfish_supporting_reads=sawFish2.out.sv_supporting_reads
}

workflow STR {
    take:
    aligned

    main:
    if (params.genome=="hg38") {
        trgt4_all(aligned)
    }

    trgt4_diseaseSTRs(aligned)
    trgt4_diseaseSTRs.out.trgt_full.combine(puretargetPlotGenes)
    |map {meta,bam,bai,vcf,tbi,genes -> 
    tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi,strID:genes])}
    |set {trgt_plot_ch}
    trgt4_diseaseSTRs_plots(trgt_plot_ch)

    emit:
    str4_vcf=trgt4_diseaseSTRs.out.str4_vcf
}


workflow QC {
    take:
    aligned

    main:
    mosdepthROI(aligned)
    nanoStat(aligned)

    emit:
    mosdepth=mosdepthROI.out.multiqc
    nanoStat=nanoStat.out.multiqc
}

workflow {
    if (params.test ||params.summary) {
        finalUbamInput.view()
        samplesheet_full.view()
        write_input_summary(ubam_size_summary_ch)
        write_analyzed_samples_summary(ubam_size_keep_ch)
        write_dropped_samples_summary(ubam_size_dropped_ch)
        symlinks_ubam_dropped(ubam_ss_merged_size_split.drop)
    }

    if (!params.test && !params.summary) {
        if (!params.aligned) {
            write_input_summary(ubam_size_summary_ch)
            write_analyzed_samples_summary(ubam_size_keep_ch)
            write_dropped_samples_summary(ubam_size_dropped_ch)
            symlinks_ubam_dropped(ubam_ss_merged_size_split.drop)
            PREPROCESS(finalUbamInput)

            if (!params.failedReads && !params.allReads && !params.hifiReads) {
                extractHifi(PREPROCESS.out.alignedAll)
                extractHifi.out.alignedHifi.join(PREPROCESS.out.alignedAll)
                | map {meta,bamHifi,baiHifi,bamAll,baiAll ->
                tuple(meta, [mainBamFile:bamHifi, mainBaiFile:baiHifi, bamAll:bamAll, baiAll:baiAll])}
                |set {alignedFinal}
            }
            if (params.allReads || params.hifiReads || params.failedReads) {
                PREPROCESS.out.alignedAll
                | map {meta,bamAll,baiAll ->
                tuple(meta,[mainBamFile:bamAll,mainBaiFile:baiAll])}
                |set {alignedFinal}
            }
        }

        if (!params.skipQC) {
            QC(alignedFinal)
        }
        
        if (!params.skipVariants) {
            VARIANTS(alignedFinal)
            VARIANTS.out.dv_vcf     //meta, vcf, idx
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            |set {dv_vcf}
            VARIANTS.out.dv_gvcf
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            |set {dv_gvcf}
        }

        if (!params.skipSV) {
            STRUCTURALVARIANTS(alignedFinal)
            STRUCTURALVARIANTS.out.sawfish_vcf //meta, vcf, idx
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            |set {sawfish_ch}
        }

        if (!params.skipSTR) {
            STR(alignedFinal)
        }

        if (!params.skipVariants && !params.skipSV && !params.skipSTR) {

            STR.out.str4_vcf
            | map {meta,vcf,idx -> tuple(meta,[vcf,idx])}
            | set {strchannel}

            alignedFinal.join(dv_vcf).join(sawfish_ch).join(strchannel)
            |set {hiphaseInput}

            hiPhase(hiphaseInput)
            
            hiPhase.out.hiphase_bam
            .join(hiPhase.out.hiphase_dv_vcf)
            .join(hiPhase.out.hiphase_sawfish_vcf)
            .join(STRUCTURALVARIANTS.out.sawfish_supporting_reads)
            | map {meta,bam,bai,dv_vcf,dv_idx,sv_vcf,sv_idx,sv_jsonReads -> 
            tuple(meta,[bam:bam,bai:bai,dv_vcf:dv_vcf,dv_idx:dv_idx,sawfish_vcf:sv_vcf,sawfish_idx:sv_idx,sawfish_reads:sv_jsonReads])}
            |set {phasedAll}    // use for val(data) instead of path(data) setup in modules 

            if (params.jointCall || params.jointSS) {
                STRUCTURALVARIANTS.out.sawfish_discover_dir
                | map {" --sample "+it}
                |collectFile(name: "sawfish_discover_dir_list.csv", newLine: false)
                |map {it.text.trim()}
                |set {sawfish_discover_bam_list_ch}
    
                STRUCTURALVARIANTS.out.sawfish_discover_dir2   // tuple(meta), path(dir), val(bam)
                | map { meta, dir, bam ->
                    // Emit key + a record line we’ll write into the manifest
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

                VARIANTS.out.dv_gvcf
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

                glNexus_jointCall(glnexus_manifest_ch)
                sawFish2_jointCall_caseID(sawfish_jointCall_manifest_ch)
                svdb_sawFish2_jointCall_caseID(sawFish2_jointCall_caseID.out.sv_jointCall_caseID_vcf)
            }

            methylationBW(phasedAll)
            methylationSegm(methylationBW.out)
            cramino(phasedAll)
            mitorsaw(phasedAll)
            whatsHap_stats(phasedAll)
            
            if (params.genome=="hg38") {
                paraphase(phasedAll)
                //kivvi_d4z4(phasedAll)
                kivvi05_d4z4(phasedAll)
                starphase(phasedAll)
                svTopo(phasedAll)
                svdb_SawFish(phasedAll)
            }

            hiPhase.out.hiphase_bam
            .join(svdb_SawFish.out.sawfishAF10)
            .join(STRUCTURALVARIANTS.out.sawfish_supporting_reads)
            | map {meta,bam,bai,sv10_vcf,sv10_idx,sv_jsonReads -> 
            tuple(meta,[bam:bam,bai:bai,sawfish10_vcf:sv10_vcf,sawfish10_idx:sv10_idx,sawfish_reads:sv_jsonReads])}
            |set {phasedSawfishAF10}   

            svTopo_filtered(phasedSawfishAF10)

            // trio specific analysis. 
            //NB: Currently only works for single-family or single-trio analysis!

            if (params.hpo && params.samplesheet && (params.jointCall || params.jointSS)) {
            glNexus_jointCall.out.glnexus_vcf.combine(hpo_ch).combine(samplesheet_path_ch)
            |set {genomiser_ch}
            glNexus_jointCall.out.glnexus_wes_roi_vcf.combine(hpo_ch).combine(samplesheet_path_ch)
            |set {exomiser_ch}
            svdb_sawFish2_jointCall_caseID.out.sawfish_caseID_AF10.combine(hpo_ch).combine(samplesheet_path_ch)
            |set {exomiserSV_ch}
                //above structure: caseID, vcf, idx, hpoFile,samplesheet
                exo14_2508_exome(exomiser_ch)
                exo14_2508_genome(genomiser_ch)
                exo14_2508_SV(exomiserSV_ch)
            }

            if (!params.skipQC) {

                Channel.empty()
                .mix(QC.out.mosdepth)
                .mix(QC.out.nanoStat)
                .mix(whatsHap_stats.out.multiqc)
                .map { meta, qcfile ->
                    tuple(params.multiqcKey(meta), meta, qcfile)
                }
                .groupTuple(by: 0)
                .map { key, metas, qcfiles ->

                    // pick one representative meta for publishDir + naming
                    def meta0 = metas.find { it.relation == 'index' } ?: metas[0]

                    tuple(meta0, qcfiles)
                }
                .set { multiqc_inputs_ch }
                multiQC(multiqc_inputs_ch)
            }
        }
    }
}

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

