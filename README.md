# KG Vejle Germline PacBio LRS pipeline

## General info:
This pipeline is used for PacBio germline WGS at Clinical Genetics, Vejle


## Default analysis steps and tools used:

- Alignment (pbmm2)
- Small variants (DeepVariant)
- Structural variants (Sawfish)
- Repeat expansions (TRGT)
- Repeat contraction (Kivvi)
- Phasing (hiPhase)
- Pseudogenes (Paraphase)
- Pharmacogenomics (Starphase)
- Mitochondrial variants (mitorsaw)
- Methylation profiles (pb-cpg-tools and methBat)
- QC module (nanostat, mosdepth, cramino, whatsHap, multiQC)

## Additional user-defined output
- Exomiser will be run, if the user provides a file with hpo terms (e.g. for rare disease trio analysis). Three exomiser runs based on smallvariants (jointGenotyped DeepVariant vcf), structural variants (jointGenotyped sawfish vcf) 



# Usage

The tools used and output generated depends on how the pipeline is run. See below for instructions.
