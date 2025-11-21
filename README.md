# KG Vejle Germline PacBio LRS pipeline

## General info:
This pipeline is used for PacBio germline WGS at Clinical Genetics, Vejle


## Default analysis steps and tools used

- Alignment (pbmm2)
- Small variants (DeepVariant)
- Structural variants (Sawfish)
- Inhouse allele frequency annotation of structural variants (SVDB)
- Repeat expansions (TRGT)
- Repeat contraction (Kivvi)
- Phasing (hiPhase)
- Pseudogenes (Paraphase)
- Pharmacogenomics (Starphase)
- Mitochondrial variants (mitorsaw)
- Methylation profiles (pb-cpg-tools and methBat)
- QC module (nanostat, mosdepth, cramino, whatsHap, multiQC)

## Additional user-defined output
- Exomiser will be included based on smallvariants (jointGenotyped DeepVariant vcf) and structural variants (jointGenotyped sawfish vcf), if the user provides a file with hpo terms (e.g. for rare disease trio analysis).
- JointGenotyping is disabled by default, but can be activated with --jointCall (see parameter section)
- Grouped output (e.g. collect data for all samples for each tool in a single outputfolder based on the used samplesheet) is disabled by default (i.e. output data is collected per sample by default). Grouped output can be activated with --groupedOutput (see parameter section).
- Tools and modules can be disabled using e.g. --skipQC, --skipVariants, --skipSV, --skipSTR (see parameter section)

# Usage

The tools used and output generated depends on how the pipeline is run. See below for instructions.
The script requires a samplesheet as input:

## Samplesheet format, unrelated samples
The most basic samplesheet contains 3 tab-separated columns in this specific order:

    CASE_GROUP NPN GENDER

Where CASE_GROUP can be either the NPN for unrelated samples, or e.g. contain a groupID for samples that should be analyzed together, e.g. "WCS_CNV", "TRIO_NAME" etc.

Example: Unrelated samples, separate output for each sample, use NPN as CASE_GROUP, so each sampleoutput is stored in an output folder named {NPN}

    123456789012  123456789012   female
    234567890123  234567890123   male
    345678901234  345678901234   male

Example: Unrelated samples, but collect sampleoutput per group based on values in CASE_GROUP

    WGS_CNV      123456789012    female
    WGS_CNV      234567890123    male
    Pseudogene   345678901234    male
    Pseudogene   456789012345    female
    ManualKey    567890123456    male
    ManualKey    678901234567    female


When using the above samplesheet with the --groupedOutput option, the output will be separated into WGS_CNV, Pseudogene and ManualKey.

## Samplesheet format, trios

    CASEID  NPN  GENDER  RELATION  AFFECTED_STATUS

Example:

    trio_name	113648565123	female	mater	normal
    trio_name	345678965123	female	index	affected
    trio_name	123456789123	male	pater	normal

For trios, if --hpo is used, the script will generate a pedigree file (.ped) and run exomiser for the trio, using the information in the samplesheet. Make sure to have each field set correctly!

Note: GENDER should be male/female, RELATION should be mater/index/pater and AFFECTED_STATUS should be normal/affected/unknown

## Options and parameters:
    --help                  Show this help menu with available options
    --input [path]:         Path to data to use as input. Default: Search KG Vejle archive for input unmapped bams (search across all previous PacBio runs)
    --allReads [bool]:      Use allreads, i.e. HiFi reads and failed reads as input. Default: Only uses HiFi reads as input.
    --noMerge [bool]:       Do not use all available unmapped bam files as input. Only use files in folder set with --input). Requires that the --input parameter is set
    --samplesheet [path]:   path to samplesheet to use. Required
    --skipQC [bool]:        Do not run QC module
    --skipVariants [bool]:  Do not call small variants (i.e. skip DeepVariant)
    --skipSV [bool]:        Do not call structural variants (i.e. skip Sawfish)
    --skipSTR [bool]:       Do not call repeat expansions (i.e. skip TRGT and Kivvi)
    --groupedOutput [bool]: Merge output based on value in first column of samplesheet
    --jointCall [bool]:     Perform joint genotyping of samples based on value in first column of samplesheet
    --hpo [path]:           Path to file with hpo terms (only relevant for trios / family analyses)

NOTE:
If any of the parameters --skipVariants, --skipSV or --skipSTR is set, phasing of the data is disabled. 
In the current version of the script, phasing requires the output of deepVariant, Sawfish and TRGT to phase the data properly. This may be changed in future versions to allow e.g. phasing based solely on deepvariant.

## Execution plan
The script can be run on a single compute node (local), or using KG Vejles SLURM cluster
The script is run locally by default, but can use the SLURM cluster by adding "-profile slurm" to the commandline. Note that the "-profile" is a built in function of Nextflow, i.e. it is set using a single "-" (-profile instead of --profile)

## Usage examples

#### Default: Analyze all samples in samplesheet. Use all unmapped bam files available (across multiple SMRTcells) for each sample. Run all default analysis steps:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt

#### Default: Same as above, but use SLURM to execute the script:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt -profile slurm

#### Analyze all samples in samplesheet. Group output for all samples. Run joint genotyping for DeepVariant and Sawfish, use SLURM to execute the script:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt --jointCall --groupedOutput -profile slurm

#### Trio analysis: Run exomiser in addition to all default analysis steps (requires the use of --hpo, and that the samplesheet contains the trio samples with gender, relation and affection status set):
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/TrioSamplesheet.txt --jointCall --groupedOutput --hpo /path/to/TrioHPOfile.txt


# Output

Based on the options shown above, and the samplesheet used, output is either stored per sample (individual tools as subfolders for each sample), or grouped by tools and analysis (e.g. all deepVariant data for all samples stored in a single outputfolder).
See infonet for further details.


