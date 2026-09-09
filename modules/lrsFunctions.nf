#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


def sexFromGender(gender, sampleId) {
    def g = "${gender}".trim().toUpperCase()
    if (g == 'M') return 'male'
    if (g == 'K') return 'female'
    throw new IllegalArgumentException(
        "Sample ${sampleId}: gender field is '${gender}', expected M or K. " +
        "Fix the samplesheet — the pipeline will not guess a karyotype.")
}


/**
 * Parse a LabWare metadata line into a meta map.
 *
 * Accepts both layouts, so no preprocessing step is needed:
 *   underscore-joined (as emitted):  0000103212_113624121888_78_SL-NGC-SJAELDNE_K_T_113782522715
 *   tab-exploded (tr '_' '\t'):      same fields, one per column
 *   [0] rekv   [1] npn     [2] material  [3] testlist
 *   [4] gender [5] proband [6] intRef    [7] expectedCount (optional)
 */
def parseMetaLine(line) {

    def cols = line.toString().trim().split('\t') as List
    def f    = (cols.size() == 1) ? cols[0].tokenize('_') : cols

    if (f.size() < 7) {
        throw new IllegalArgumentException(
            "Metadata line '${line}' has ${f.size()} fields, expected 7 or 8 " +
            "(rekv_npn_material_testlist_gender_proband_intRef[_expectedCount]).")
    }

    return [
        rekv          : f[0].toString().trim(),
        id            : f[1].toString().trim(),
        material      : f[2].toString().trim(),
        testlist      : f[3].toString().trim(),
        gender        : f[4].toString().trim(),
        sex           : sexFromGender(f[4], f[1]),
        proband       : f[5].toString().trim(),
        intRef        : f[6].toString().trim(),
        expectedCount : (f.size() > 7) ? f[7].toString().trim() : null,
    ]
}

def parseUbamNpn(bamPath) {
    def base = bamPath.toString().tokenize('/').last().replaceFirst(/\.bam$/, '')
    def dots = base.tokenize('.')

    if (dots.size() < 3) return null

    def npn = dots[0].toString().tokenize('_')[0]

    return npn ? npn.toString() : null
}




def karyotypeFlag(meta) {
    if (meta.sex == 'male')   return '--karyotype XY'
    if (meta.sex == 'female') return '--karyotype XX'
    throw new IllegalStateException(
        "Sample '${meta.id}': meta.sex='${meta.sex}'. Sex should have been validated " +
        "at samplesheet parse time — refusing to guess a karyotype.")
}


/** Sawfish --expected-cn flag. Fail-closed, same contract as karyotypeFlag. */
def expectedCnFlag(meta) {
    if (meta.sex == 'male')   return "--expected-cn ${params.sawfishExpectedCnXY}"
    if (meta.sex == 'female') return "--expected-cn ${params.sawfishExpectedCnXX}"
    throw new IllegalStateException(
        "Sample '${meta.id}': meta.sex='${meta.sex}' — refusing to guess expected copy number.")
}


/*
 * Filename prefixes. Every published file is named
 *     <sampleOrCase>.<genomeVersion>.<tag>.<...>
 * Pick the variant by which read set produced the file:
 *   prefixHifi      HiFi-derived output              (params.tagHifi)
 *   prefixFail      fail-read-derived output         (params.tagFail)
 *   prefixStr       TRGT output, both read sets      (params.strTag)
 *   prefixBase      no read-set tag (mosdepth ROI, hiphase stats)
 *   prefixCaseHifi  as prefixHifi but keyed on meta.caseID, for joint calls
 *
 *   Example usage:
 *     output: path("${prefixHifi(meta)}.hiphase.paraphase/*")
 *     script: def prefix = prefixHifi(meta)
 */
def prefixBase(meta)     { "${meta.id}.${params.genomeVersion}" }
def prefixHifi(meta)     { "${meta.id}.${params.genomeVersion}.${params.tagHifi}" }
def prefixFail(meta)     { "${meta.id}.${params.genomeVersion}.${params.tagFail}" }
def prefixStr(meta)      { "${meta.id}.${params.genomeVersion}.${params.strTag}" }
def prefixCaseHifi(meta) { "${meta.caseID}.${params.genomeVersion}.${params.tagHifi}" }