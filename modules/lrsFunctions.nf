#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


def sexFromGender(gender, sampleId) {
    switch ((gender ?: '').toString().trim().toLowerCase()) {
        case ['m', 'male']:
            return 'male'
        case ['k', 'f', 'female']:
            return 'female'
        default:
            throw new IllegalArgumentException(
                "Sample ${sampleId}: gender field is '${gender}', expected M/K " +
                "(or male/female). Fix the samplesheet — the pipeline will not " +
                "guess a karyotype.")
    }
}


def parseMetaLine(line) {

    def cols = line.toString().trim().split('\t') as List
    def f    = (cols.size() == 1) ? cols[0].tokenize('_') : cols

    if (f.size() < 6) {
        throw new IllegalArgumentException(
            "Metadata line '${line}' has only ${f.size()} fields; expected at least " +
            "rekv_npn_material_testlist_gender_proband.")
    }

    def g = -1
    for (int i = 3; i < f.size() - 1; i++) {
        if (f[i].toString().trim().toUpperCase() in ['M', 'K'] &&
            f[i + 1].toString().trim().toUpperCase() in ['T', 'F']) {
            g = i
            break
        }
    }

    if (g < 0) {
        throw new IllegalArgumentException(
            "Metadata line '${line}': could not locate the gender/proband pair " +
            "(expected M or K followed by T or F). Either the line is malformed or " +
            "the encoding has changed — refusing to guess which field is the gender.")
    }

    def npn = f[1].toString().trim()
    if (!npn) throw new IllegalArgumentException("Metadata line '${line}' has an empty NPN.")

    return [
        rekv          : f[0].toString().trim(),
        id            : npn,
        material      : f[2].toString().trim(),
        testlist      : f[3..(g - 1)].collect { p -> p.toString().trim() }.join('_'),
        gender        : f[g].toString().trim(),
        sex           : sexFromGender(f[g], npn),
        proband       : f[g + 1].toString().trim(),
        intRef        : (f.size() > g + 2) ? f[g + 2].toString().trim() : null,
        expectedCount : (f.size() > g + 3) ? f[g + 3].toString().trim() : null,
    ]
}

def parseUbamNpn(bamPath) {
    def base = bamPath.toString().tokenize('/').last().replaceFirst(/\.bam$/, '')
    def dots = base.tokenize('.')

    if (dots.size() < 3) return null

    def npn = dots[0].toString().tokenize('_')[0]

    return npn ? npn.toString() : null
}
