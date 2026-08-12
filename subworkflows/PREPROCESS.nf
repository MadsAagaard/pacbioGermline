#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {create_fofn;
        inputFiles_symlinks_ubam;
        pbmm2_align_hifi;
        pbmm2_align_failedOnly;
        extractHifi;
        } from "../modules/dnaModules.nf" 



workflow PREPROCESS {

    take:
    finalUbamInput     
   
    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    create_fofn(finalUbamInput)

    pbmm2_align_hifi(create_fofn.out.hifiReads)
    pbmm2_align_failedOnly(create_fofn.out.failReads)
    
    pbmm2_align_hifi.out.bam //val(meta), path(bam),path(bai)
    .join(pbmm2_align_hifi.out.bam)
    | map {meta,hifiBam,hifiBai,failBam,failBai -> 
        tuple(meta, [
                hifiBam:hifiBam,
                hifiBai:hifiBai
                failBam:failBam,
                failBai:failBai                    
            ])
    }
    |set {alignedFinal_ch}

    emit:
    alignedFinal=alignedFinal_ch

}



/*
workflow PREPROCESS {

    take:
    finalUbamInput     
   
    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    create_fofn(finalUbamInput)

    pbmm2_align_hifi(create_fofn.out.hifiReads)
    pbmm2_align_failedOnly(create_fofn.out.failReads)
    
    if (!params.failedReads && !params.allReads && !params.hifiReads) {
        extractHifi(pbmm2_align_mergedData.out.bamAll)
        extractHifi.out.alignedHifi
            .join(pbmm2_align_mergedData.out.bamAll)
            | map {meta,bamHifi,baiHifi,bamAll,baiAll ->
            tuple(meta, [mainBamFile:bamHifi, mainBaiFile:baiHifi, bamAll:bamAll, baiAll:baiAll])}
            | set {alignedFinal_ch}
    }
    if (params.allReads || params.hifiReads || params.failedReads) {
        pbmm2_align_mergedData.out.bamAll
            | map {meta,bamAll,baiAll ->
            tuple(meta,[mainBamFile:bamAll,mainBaiFile:baiAll])}
            | set {alignedFinal_ch}
    }

    emit:
    alignedFinal=alignedFinal_ch

}
*/