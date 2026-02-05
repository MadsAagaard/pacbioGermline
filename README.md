# KG Vejle Germline PacBio LRS pipeline

## General info:
This pipeline is used for PacBio LRS germline WGS at Clinical Genetics, Vejle


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
- JointGenotyping is disabled by default for single genome analyses, but can be added with --jointCall or when analyzing family or other multisample analyses using --jointSS (see parameter section)
- Joint genotyping uses GLNexus for small variants and Sawfish for structural variants
- Tools and modules can be disabled using e.g. --skipQC, --skipVariants, --skipSV, --skipSTR (see parameter section)

# Usage

The tools used and output generated depends on how the pipeline is run. See below for instructions.
The script requires a samplesheet as input:

## Default samplesheet format (based on LabWare metadata output)
The default samplesheet is derived from the labWare metadata output generated per run.
The samplesheet is structured as single column, 1 sample per row. Each entry consists of 7 fields separated by underscore ("_"):
1. rekvisition
2. NPN / sampleID
3. Bio. material
4. Analysis testlist
5. Gender (M=male, K=female)
6. Proband (T=true, F=false)
7. Internal reference (if absent, i.e. single sample: "noInfo")

Example:
0000012345_123456789012_11_SL-LWG-CNV_K_F_noInfo


## Custom samplesheet format, unrelated samples
A custom samplesheet can be envoked with the --intSS parameter (see parameter section)

A custom samplesheet consists of at least 3 tab separated columns, in this particular order:

    CASE_GROUP NPN GENDER

Where CASE_GROUP can be either the NPN for unrelated samples, or e.g. contain a groupID for samples that should be analyzed together, e.g. "WCS_CNV", "TRIO_NAME" "PROJECT_A" etc.

Example: Unrelated samples, but collect sampleoutput per group based on values in CASE_GROUP

    WGS_CNV      123456789012    female
    WGS_CNV      234567890123    male
    Pseudogene   345678901234    male
    Pseudogene   456789012345    female
    ManualKey    567890123456    male
    ManualKey    678901234567    female

When using the above samplesheet with the --jointSS option, the output will be separated into WGS_CNV, Pseudogene and ManualKey.

## Custom samplesheet format, trios
A custom samplesheet for family analysis /trios etc., requires the same 3 columns described above, and an additional 2 columns describing the relation and affection status, e.g.
:
    CASEID  NPN  GENDER  RELATION  AFFECTED_STATUS

Example:

    trioID	113648565123	female	mater	normal
    trioID	345678965123	female	index	affected
    trioID	123456789123	male	pater	normal

For trios, if --hpo is used, the script will generate a pedigree file (.ped) and run exomiser for the trio, using the information in the samplesheet. Make sure to have each field set correctly!

Note: GENDER should be either male/female, RELATION should be either mater/index/pater and AFFECTED_STATUS should be one of normal/affected/unknown

## Options and parameters:
    --help                  Show this help menu with available options
    
    --samplesheet   [path]: Path to samplesheet to use. Required
    
    --jointSS       [bool]: Use jointGenotyping, and group output for all samples (use for trio and family analysis)
                                Default: Not set
    
    --intSS         [bool]: Use custom samplesheet format (tab separated, 3 mandatory columns)
                                Default: Not set

    --input         [path]: Path to data to use as input. 
                                Default: Not set. Instead, Search KG Vejle archive for input unmapped bams (search across all previous PacBio runs)
    
    --allReads      [bool]: Use allreads, i.e. HiFi reads and failed reads as input.
                                Default: Uses a combination of AllReads (for STR analysis) and HiFi reads for everything else
    
    --hifiReads     [bool]: Use HiFi reads only, for all analysis.
                                Default: Uses a combination of AllReads (for STR analysis) and HiFi reads for everything else

    --singleOnly    [bool]: Only analyze single genomes in samplesheet (i.e. "noInfo" in int ref)
                                Default: Not set - analyze all samples in samplesheet

    --skipQC        [bool]: Do not run QC module
                                Default: Not set

    --skipVariants  [bool]: Do not call small variants (i.e. skip DeepVariant)
                                Default: Not set

    --skipSV        [bool]: Do not call structural variants (i.e. skip Sawfish)
                                Default: Not set

    --skipSTR       [bool]: Do not call repeat expansions (i.e. skip TRGT and Kivvi)
                                Default: Not set

    --jointCall     [bool]: Perform joint genotyping of samples based on value in first column of samplesheet
                                Default: Not set

    --hpo           [path]: Path to file with hpo terms (only relevant for trios / family analyses)
                                Default: Not set

    --minGB         [int]:  Minimum size (in gigabytes) of all input unmapped bam files pr. sample.
                                Default: 36 GB for allreads (HiFi + failed), 30GB for HiFi reads

    ### Slurm Execution parameters:
    -profile slurm:         Run pipeline using KGVejle SLURM cluster
                                Default: Run pipeline on local server (where scrit is started)
    --slurmA        [bool]: Use secondary fast tmp storage (nfs_fast_a)
                                Default: Use primary fast tmp storage location at KGVejle                                        



NOTE:
If any of the parameters --skipVariants, --skipSV or --skipSTR is set, phasing of the data is disabled. 
In the current version of the script, HiPhase requires the output of DeepVariant, Sawfish and TRGT to phase the data properly. This may be changed in future versions to allow e.g. phasing based solely on DeepVariant.

## Script executor - local or SLURM
The script can be run on a single compute node (local), or using KG Vejles SLURM cluster
The script is run locally by default, but can use the SLURM cluster by adding "-profile slurm" to the commandline. Note that the "-profile" is a built in function of Nextflow, i.e. it is set using a single "-" (-profile instead of --profile)

## Usage examples

#### Default: Analyze all samples in default samplesheet. Use all unmapped bam files available (across multiple SMRTcells) for each sample. Run all default analysis steps:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt

#### Default: Same as above, but use SLURM to execute the script:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt -profile slurm

#### Same as above, but skip QC and Structural variantcalling:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt  -profile slurm --skipQC --skipSV

#### Default: Trio analysis, run exomiser using hpo.txt:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt -profile slurm --hpo /path/to/hpo.txt --jointSS

#### Analyze all samples in custom samplesheet. Run joint genotyping for DeepVariant and Sawfish, use SLURM to execute the script:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt --jointCall --profile slurm

#### Analyze samples in default samplesheet. Use only unmapped bam files available in subfolders under /input: for each sample. Run all default analysis steps:
    nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt --input /path/to/selected/rawData/




# Output

Based on the options shown above, and the samplesheet used, output is either stored per sample (individual tools as subfolders for each sample), or grouped by tools and analysis (e.g. all DeepVariant data for all samples stored in a single outputfolder).

## Default output structure

Per sample, 


See KG Vejle infonet for further details.


