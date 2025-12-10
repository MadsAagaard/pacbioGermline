#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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
        if (!params.allReads && !params.failedReads){
            inputBam="${params.input}/**/*.hifi_reads.*.bam"
        }
        if (params.allReads){
            inputBam="${params.input}/**/*.bam"
        }
        if (params.failedReads){
            inputBam="${params.input}/**/*.fail_reads.*.bam"
        }
    }
    
    if (!params.input) {
        if (!params.allReads && !params.failedReads){
            inputBam="${params.dataArchive}/**/*.hifi_reads.*.bam"
        }
        if (params.allReads){
            inputBam="${params.dataArchive}/**/*.bam"
        }
        if (params.failedReads){
            inputBam="${params.dataArchive}/**/*.fail_reads.*.bam"
        }
    }

    /*
        Different naming schemes during imnplementation

        params.oldSS (initial setup during testing): 
        - samplesheet contains 3 cols per default: caseid, npn, gender
        - input files are named npn.runinfo.readset.barcode.bam

        From 251107 to dec 2025:
        params.intSS (intermediate samplesetup):
        - samplesheet contains same 3 cols, or 5 cols for trio/family
        - input files are named with more comprehensive sampleInfo:
        Input file structure:
        sampleInfo.runinfo.readset.barcode.bam, where sampleInfo contains:
        npn_material_testlist_gender


        From december 2025 (run 251205 and onwards):
        default going forward, if params.oldSS or intSS not set:

        - samplesheet extracted directly from run samplesheets:
        Strucutre:
        rekv_npn_material_testlist_gender_proband(T or F)_internalRef (number or noInfo)

        Examples:
        0000093123_113720841930_78_SL-NGC-SJAELDNE_K_F_noInfo
        0000093645_113697182321_78_SL-NGC-NYRESVIGT_M_F_113632562421
        0000093646_113697182186_78_SL-NGC-NYRESVIGT_K_F_113632562421
        Input bam:
        same as for intSS:
        sampleInfo.runinfo.readset.barcode.bam, where sampleInfo contains:
        npn_material_testlist_gender
        Example:
        113618066447_78_NGC-NEUROGENETIK_K.m84313_251205_131758_s1.hifi_reads.bc2049.bam
    */
   
    // default from dec. 5th, 2025:

    if (params.samplesheet && !params.oldSS && !params.intSS) {

        // new samplesheet - directly from metadata extracted from LabWare:
        channel.fromPath(params.samplesheet)
        | splitCsv(sep:'\t')
        |map { row ->
             (rekv, npn,material,testlist,gender,proband,intRef) = row[0].tokenize("_")
            meta=[id:npn,caseID:testlist, sex:gender, proband:proband,intRef:intRef, rekv:rekv]
            meta
            }
        | set {samplesheet_full}



        Channel.fromPath(inputBam, followLinks: true)
        |map { tuple(it.baseName,it) }
        |map {id,bam -> 
                (samplenameFull,pacbioID,readset,barcode)   =id.tokenize(".")
                (instrument,date,time)                      =pacbioID.tokenize("_")     
                (samplename,material,testlist,gender)       =samplenameFull.tokenize("_")
                meta=[id:samplename,gender:gender]
                tuple(meta,bam)        
            }
        |groupTuple(sort:true)
        |branch  {meta,bam -> 
            UNASSIGNED: (meta.id=~/UNASSIGNED/)
                        return [meta,bam]
            samples: true
                        return [meta,bam]
        }
        | set {ubam_input }
    }

    // intermediate naming scheme:
    if (params.samplesheet && !params.oldSS && params.intSS) {
        channel.fromPath(params.samplesheet)
        | splitCsv(sep:'\t')
        |map { row -> 
                (caseID, samplenameFull) =tuple(row)
                (samplename,material,testlist,gender)       =samplenameFull.tokenize("_")
            meta=[id:samplename,caseID:caseID, sex:gender, testlist:testlist]
            meta
            }
        | set {samplesheet_full}

        Channel.fromPath(inputBam, followLinks: true)
        |map { tuple(it.baseName,it) }
        |map {id,bam -> 
                (samplenameFull,pacbioID,readset,barcode)   =id.tokenize(".")
                (instrument,date,time)                      =pacbioID.tokenize("_")     
                (samplename,material,testlist,gender)       =samplenameFull.tokenize("_")
                meta=[id:samplename,gender:gender,rundate:date,testlistFile:testlist]
                tuple(meta,bam)        
            }
        |groupTuple(sort:true)
        |branch  {meta,bam -> 
            UNASSIGNED: (meta.id=~/UNASSIGNED/)
                        return [meta,bam]
            samples: true
                        return [meta,bam]
        }
        | set {ubam_input }
    }

    // Old (initial) naming scheme:
    if (params.samplesheet && params.oldSS && !params.intSS) {

        channel.fromPath(params.samplesheet)
            | splitCsv(sep:'\t')
            |map { row -> 
                (caseID, samplename, sex) =tuple(row)

                meta=[caseID:caseID,id:samplename,sex:sex]
                meta
                }
            | set {samplesheet_full}
  

        Channel.fromPath(inputBam, followLinks: true)
        |map { tuple(it.baseName,it) }
        |map {id,bam -> 
                (samplename,pacbioID,hifi,barcode)      =id.tokenize(".")
                (instrument,date,time)                  =pacbioID.tokenize("_")     
                meta=[id:samplename,gender:"NA"]
                tuple(meta,bam)        
            }
        |groupTuple(sort:true)
        |branch  {meta,bam -> 
            UNASSIGNED: (meta.id=~/UNASSIGNED/)
                        return [meta,bam]
            samples: true
                        return [meta,bam]
        }
        | set {ubam_input }
    }


    if (params.samplesheet) {

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
        |branch  {meta,bam -> 
            UNASSIGNED: (meta.id=~/UNASSIGNED/)
                        return [meta,bam]
            samples: true
                        return [meta,bam]
        }
        | set {ubam_input }
    }

    ubam_input.samples
        | map { meta, bam -> tuple(meta.id,meta,bam) }
        | set {ubam_input_samples}    


    if (params.samplesheet) {

        samplesheet_full
        |map {row -> meta2=[row.id,row]}
        |set {samplesheet_join}

        samplesheet_join.join(ubam_input_samples)
        |map {samplename, metaSS, metaData, bam -> tuple(metaSS+metaData,bam)}
        |view
        |set {finalUbamInput}
    }

    if (!params.samplesheet) {
        ubam_input.samples
        |set {finalUbamInput}
    }
}



/////////////////// MODULES ///////////////////////
include {pbmm2_align;
        create_fofn;
        pbmm2_align_mergedData;
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
        //collect_versions;
        } from "./modules/dnaModules.nf" 


puretargetPlotGenes=["SCA1_ATXN1",
                     "SCA2_ATXN2",
                     "SCA3_ATXN3",
                     "SCA6_CACNA1A",
                     "SCA7_ATXN7",
                     "CANVAS_RFC1",
                     "DM1_DMPK",
                     "DM2_CNBP",
                     "FTDALS1_C9orf72",
                     "FXS_FMR1",
                     "FRDA_FXN",
                     "HD_HTT"]


////////////////// WORKFLOWS AND PROCESSES ///////////////////////

workflow PREPROCESS {

    take:
    finalUbamInput     
   
    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    if (params.noMerge) {
        pbmm2_align(finalUbamInput)
        
        pbmm2_align.out.bam
        |set {alignedTMP}
    }
    if (!params.noMerge) {
        create_fofn(finalUbamInput)
        pbmm2_align_mergedData(create_fofn.out)
        pbmm2_align_mergedData.out.bam
        |set {alignedTMP}
    }

    emit:
    aligned=alignedTMP
    
}

workflow VARIANTS {

    take:
    aligned     //tuple(meta,[bam,bai])
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
    //tuple(meta,bam,genes)}
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


/*
workflow PHASED {

    take:
    phasedAll

    main:

    methylationBW(phasedAll)
    methylationSegm(methylationBW.out)
    cramino(phasedAll)
    mitorsaw(phasedAll)
    whatsHap_stats(phasedAll)

    if (params.genome=="hg38") {
        paraphase(phasedAll)
        kivvi_d4z4(phasedAll)
        starphase(phasedAll)
        svTopo(phasedAll)
        svdb_SawFish(phasedAll)
    }

    emit:


}
*/



//Channel.topic('versions') as versions_ch
workflow {
    if (params.test) {
        finalUbamInput.view()
        samplesheet_full.view()
    }

    if (!params.test) {
        if (!params.aligned) {

            PREPROCESS(finalUbamInput)

            PREPROCESS.out.aligned
            | map {meta,bam,bai -> tuple(meta,[bam,bai])}
            |set {alignedFinal}
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




            if (params.jointCall) {
                STRUCTURALVARIANTS.out.sawfish_discover_dir
                | map {" --sample "+it}
                |collectFile(name: "sawfish_discover_dir_list.csv", newLine: false)
                |map {it.text.trim()}
                |set {sawfish_discover_bam_list_ch}
    
                STRUCTURALVARIANTS.out.sawfish_discover_dir2 //meta, sawfishDir,bam,bai
                | map {meta, dir, bam ->
                    def dirPath = dir.toString()
                    def bamPath = bam.toString()
                    return [ meta.caseID, dirPath+", "+ bamPath ]
                }
                | collectFile(newLine: true) { item  ->
                    def caseID = item[0]
                    def line = item[1]
                    return [ "${caseID}.sawFishJoinCall.manifest.csv", line ]
                }
                | map { manifestFile -> 
                    def caseID = manifestFile.getName().tokenize(".")[0]
                    return tuple(caseID, [manifestFile])
                    }
                | set { sawfish_jointCall_manifest_ch }


                manifestChannel =dv_gvcf
                | map { meta, files ->
                    def vcfPath = files[0].toString()
                    return [ meta.caseID, "${vcfPath}" ]
                }
                | collectFile(newLine: true) { item ->
                    def caseID = item[0]
                    def line   = item[1]
                    return [ "${caseID}.manifest", line ]
                }
                
                manifestChannel
                | map { manifestFile -> manifestFile
                    def caseID = manifestFile.getName().tokenize(".")[0]
                    return tuple(caseID, [manifestFile])
                }
                | set { glnexus_manifest_ch }

                if (!params.groupedOutput) {           
                    sawFish2_jointCall_all(sawfish_discover_bam_list_ch)   
                    svdb_sawFish2_jointCall_all(sawFish2_jointCall_all.out.sv_jointCall_vcf)
                }

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
                kivvi_d4z4(phasedAll)
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

            if (params.hpo && params.samplesheet && params.jointCall && params.groupedOutput) {
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




            // input channel for multiQC
            if (!params.skipQC) {
                QC.out.mosdepth.join(QC.out.nanoStat).join(whatsHap_stats.out.multiqc)
                | map {meta,mosdepth,nanoStat,whatshap -> tuple(meta,[mosdepth,nanoStat,whatshap])}
                |set {multiqcSingleInput}   
                multiQC(multiqcSingleInput)

                def allOutputs = Channel.empty()
                allOutputs = allOutputs.mix(QC.out.mosdepth)    
                allOutputs = allOutputs.mix(QC.out.nanoStat)          
                allOutputs = allOutputs.mix(whatsHap_stats.out.multiqc)    

                allOutputs
                |groupTuple
                |set {multiqcAllInput}
                multiQC_ALL(multiqcAllInput)
            }
        }
    }
}


