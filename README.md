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
- Exomiser will be included based on smallvariants (jointGenotyped DeepVariant vcf), structural variants (jointGenotyped sawfish vcf), if the user provides a file with hpo terms (e.g. for rare disease trio analysis).
- JointGenotyping is disabled by default, but can be activated with --jointCall (see parameter section)
- Grouped output (e.g. collect data for all samples for each tool in a single outputfolder) is disabled by default (i.e. output data is collected per sample by default). Grouped output can be activated with --groupedOutput (see parameter section).
- Tools and modules can be disabled using e.g. --skipQC, --skipVariants, --skipSV, --skipSTR (see parameter section)

# Usage

The tools used and output generated depends on how the pipeline is run. See below for instructions.

## Options:
  --help
  --input
  --allReads
  --noMerge
  --samplesheet
  --skipQC
  --skipVariants
  --skipSV
  --skipSTR
  --groupedOutput
  --jointCall
  --hpo

## Usage examples:

#### Default: Analyze all samples in samplesheet. Use all unmapped bam files available (across multiple SMRTcells) for each sample. Run all default analysis steps:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt

#### Analyze all samples in samplesheet. Group output for all samples. Run joint genotyping for DeepVariant and Sawfish:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt --joiuntCall --groupedOutput

#### Trio analysis: Run exomiser in addition to all default analysis steps (requires the use of --hpo, and that the samplesheet contains the trio samples with gender, relation and affection status set):
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/TrioSamplesheet.txt --joiuntCall --groupedOutput --hpo /path/to/TrioHPOfile.txt
