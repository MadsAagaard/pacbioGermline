# KG Vejle Germline PacBio LRS pipeline

## General info

PacBio LRS germline WGS pipeline used at Clinical Genetics, Vejle.

## Read handling

HiFi reads and fail reads are aligned **separately**:

```
uBAM ─ create_fofn ─┬─ pbmm2_align_hifi       ──► HiFi BAM ──► all analyses
                    └─ pbmm2_align_failedOnly ──► fail BAM ──► TRGT only
```

TRGT receives the fail reads through `--fail-reads`; every other tool uses the
HiFi alignment. HiPhase phases both and publishes both:

- `alignments/<sample>.<genome>.HifiReads.hiphase.bam`
- `alignments/failedReads/<sample>.<genome>.failedReads.hiphase.bam`

There is no merged all-reads alignment and no `extracthifi` step.

## Default analysis steps and tools used

- Alignment (pbmm2)
- Small variants (DeepVariant)
- Structural variants (Sawfish)
- Complex SV visualisation (SVTopo)
- In-house allele frequency annotation of SVs (SVDB)
- Repeat expansions (TRGT v4 + v5; disease loci and genome-wide adotto catalog)
- Repeat contraction (Kivvi)
- Phasing (HiPhase)
- Pseudogenes (Paraphase)
- Pharmacogenomics (Starphase + PharmCAT)
- Mitochondrial variants (Mitorsaw)
- Methylation profiles (pb-CpG-tools and MethBat)
- QC (NanoStat, mosdepth, cramino, WhatsHap, MultiQC)

## Optional analysis steps (trios and related samples)

- Joint genotyping of small variants (GLNexus)
- Joint genotyping of structural variants (Sawfish)
- HPO-based gene and variant prioritisation (Exomiser / Genomiser)

---

# Usage

## Default samplesheet format (LabWare metadata output)

Single column, one sample per row, 7 underscore-separated fields:

1. rekvisition
2. NPN / sampleID
3. Bio. material
4. Analysis testlist
5. Gender (M=male, K=female)
6. Proband (T=true, F=false)
7. Internal reference (absent / single sample: `noInfo`)

```
0000012345_123456789012_11_SL-LWG-CNV_K_F_noInfo
```

**Gender must be M or K.** Anything else aborts the run at samplesheet parse
time rather than silently genotyping the sample as XX. Sawfish `--expected-cn`
and TRGT `--karyotype` depend on it.

## Custom samplesheet, unrelated samples

At least 3 tab-separated columns, in this order:

```
CASE_GROUP    NPN             GENDER
WGS_CNV       123456789012    female
WGS_CNV       234567890123    male
Pseudogene    345678901234    male
```

`CASE_GROUP` is the NPN for unrelated samples, or a group ID for samples that
should be analysed together. With `--jointSS`, output is grouped by this column.

## Custom samplesheet, trios / families

Two additional columns:

```
CASEID    NPN             GENDER    RELATION    AFFECTED_STATUS
trioID    113648565123    female    mater       normal
trioID    345678965123    female    index       affected
trioID    123456789123    male      pater       normal
```

`GENDER` must be `male`/`female`, `RELATION` one of `mater`/`index`/`pater`, and
`AFFECTED_STATUS` one of `normal`/`affected`/`unknown`. With `--hpo`, a pedigree
is generated and Exomiser runs for the family.

---

## Options and parameters

```
--help                  Show this help menu

--samplesheet   [path]: Path to samplesheet. Required (unless --input is used
                        without a samplesheet, or --familySS for family re-runs)

--input         [path]: Folder to use as input.
                            Default: search the KG Vejle archive across all runs

--jointSS       [bool]: Joint genotyping, output grouped per case
                            (use for trio and family analysis)

--customSS      [bool]: Custom samplesheet format (tab separated)

--jointCall     [bool]: Joint genotyping based on the first samplesheet column

--hpo           [path]: File with HPO terms (trios / family analyses)


### Read set selection
Default (no flag): HiFi + fail reads, aligned separately. TRGT uses both.

--hifiReads     [bool]: HiFi only. No fail alignment; TRGT runs without
                        --fail-reads.
--failedReads   [bool]: Fail reads only (debug / QC use).
--allReads      [bool]: DEPRECATED. Ignored with a warning — the default is now
                        what this used to mean.


### Sample filtering
--singleOnly    [bool]: Only samples with "noInfo" internal reference
--intrefOnly    [bool]: Only samples WITH an internal reference

--minHifiGB     [int]:  Minimum combined HiFi uBAM size per sample, in GB.
                            Default: 30
--minFailGB     [int]:  Minimum combined fail uBAM size per sample, in GB.
                            Default: 4
--minGB         [int]:  Legacy alias; overrides --minHifiGB if given.


### Skips
--skipQC        [bool]: Skip the QC module
--skipVariants  [bool]: Skip DeepVariant
--skipSV        [bool]: Skip Sawfish
--skipSTR       [bool]: Skip TRGT and Kivvi


### Safety / error handling
--errorMode     [str]:  'strict' (default) — retry transient failures, then
                            finish the run and fail loudly.
                        'cohort' — retry, then ignore, so one bad sample cannot
                            abort a large retrospective build. Failures are
                            recorded in runInfo/<date>_<ss>/FAILED_TASKS.txt.


### SLURM execution
-profile slurm:         Run on the KG Vejle SLURM cluster
                            Default: run locally on the launching node
--slurmA        [bool]: Use nfs_fast_a for the work directory instead of
                            nfs_fast_b
```

**Note:** if any of `--skipVariants`, `--skipSV` or `--skipSTR` is set, phasing
is disabled — HiPhase needs DeepVariant, Sawfish and TRGT output. A warning is
logged when this happens.

---

## Usage examples

Analyse all samples in the default samplesheet, all default steps:

```
nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt
```

Same, on SLURM:

```
nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt -profile slurm
```

Skip QC and SV calling:

```
nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt -profile slurm --skipQC --skipSV
```

Trio analysis with Exomiser:

```
nextflow run MadsAagaard/pacbioGermline -r main --samplesheet /path/to/samplesheet.txt -profile slurm --hpo /path/to/hpo.txt --jointSS
```

---

# Output

Depending on the options and samplesheet used, output is stored per sample
(tools as subfolders per sample) or grouped by case. See KG Vejle infonet for
details.
