#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * lrsFunctions.nf — shared helper functions
 * =============================================================================
 *
 * Included by name:
 *
 *     include { sexFromGender; parseMetaLine; parseUbamNpn } from './modules/lrsFunctions.nf'
 *
 * -------------------------------------------------------------------------
 * UNDERSCORES IN testlist BREAK POSITIONAL PARSING
 * -------------------------------------------------------------------------
 * The uBAM name blob and the LabWare metadata line are both underscore-joined,
 * and testlist is not guaranteed to be underscore-free. Older Revio output
 * contains e.g.
 *
 *     113720750803_78_NGC_NEUROGENETIK_M   ->  [0]=npn [1]=78 [2]=NGC
 *                                              [3]=NEUROGENETIK [4]=M
 *     113738820048_78_NGC-ARVELIG-KRFT_K   ->  [0]=npn [1]=78 [2]=testlist [3]=K
 *
 * Every index after material shifts. Counting from the left past index 2 is
 * therefore unsafe, and that is what these functions are built around:
 *
 *   parseUbamNpn  returns ONLY field 0 (npn), which cannot shift.
 *   parseMetaLine anchors on the gender/proband pair instead of counting.
 *
 * A NOTE ON CONCURRENCY
 * These are called from inside map closures. They are pure: no state is held
 * between calls, nothing outside is mutated. Do not add a script-level variable
 * to this file — closures on the dataflow thread pool would write to it
 * concurrently.
 *
 * A NOTE ON STRINGS
 * Everything returned goes through .toString(). Values from splitCsv and
 * GStrings are not always java.lang.String, and a GString has a different
 * hashCode than the equal String, so one used as a map key silently fails
 * lookup later.
 * =============================================================================
 */


/**
 * Normalise a raw gender field to 'male' or 'female'.
 *
 * Fail-closed: anything unrecognised throws, and there is no default. A missing
 * or malformed gender is a data-entry error, not a state the pipeline should
 * have a policy for — TRGT --karyotype and Sawfish --expected-cn both depend on
 * it being right.
 */
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


/**
 * Parse a LabWare metadata line into a meta map.
 *
 * Accepts BOTH layouts, so no preprocessing is needed:
 *
 *   underscore-joined (as emitted):
 *     0000103212_113624121888_78_SL-NGC-SJAELDNE_K_T_113782522715_3
 *
 *   tab-exploded (tr '_' '\t', for manual inspection):
 *     0000103212<TAB>113624121888<TAB>78<TAB>SL-NGC-SJAELDNE<TAB>K<TAB>...
 *
 * Detection is unambiguous: the joined blob contains no tabs and splits to
 * exactly one column, the exploded form never does.
 *
 *   rekv | npn | material | testlist | gender | proband | intRef | [count]
 *
 * PARSED BY ANCHOR, NOT BY POSITION.
 * The first three fields are read from the left (they cannot shift). The rest
 * are located by finding the gender/proband pair — gender in {M,K}, immediately
 * followed by proband in {T,F}. Anything between material and gender is the
 * testlist, however many underscores it contains.
 *
 * Counting positions instead would break on a testlist like NGC_NEUROGENETIK,
 * and would break silently in the 8-field case: an extra underscore turns a
 * 7-field line into an 8-field one, which parses "successfully" with every
 * field after material holding the wrong value.
 *
 * The trailing count field is optional (added late; older exports lack it).
 */
def parseMetaLine(line) {

    def cols = line.toString().trim().split('\t') as List
    def f    = (cols.size() == 1) ? cols[0].tokenize('_') : cols

    if (f.size() < 6) {
        throw new IllegalArgumentException(
            "Metadata line '${line}' has only ${f.size()} fields; expected at least " +
            "rekv_npn_material_testlist_gender_proband.")
    }

    // Locate gender: in {M,K}, immediately followed by proband in {T,F}.
    // Search from index 3 (the earliest a testlist could end).
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


/**
 * Extract the NPN from a PacBio uBAM filename.
 *
 *   113720750803_78_NGC_NEUROGENETIK_M.m84328_251121_123048_s1.hifi_reads.bc2017.bam
 *   \__________/
 *       npn
 *
 * ONLY the NPN is returned, and that is deliberate. It is field 0 of the name
 * blob, so it is the one value an underscore in testlist cannot displace.
 * material, testlist and gender are NOT extracted: their positions shift with
 * the testlist, so reading them would produce plausible-looking wrong values
 * rather than an error.
 *
 * The samplesheet is the sole source of sex, material and testlist. There is no
 * cross-check against the filename — a check that fires on correct samples is
 * worse than no check.
 *
 * Returns null for anything unparseable (UNASSIGNED, ad-hoc test data) so the
 * caller can filter rather than the run dying on one stray file.
 */
def parseUbamNpn(bamPath) {
    def base = bamPath.toString().tokenize('/').last().replaceFirst(/\.bam$/, '')
    def dots = base.tokenize('.')

    if (dots.size() < 3) return null

    def npn = dots[0].toString().tokenize('_')[0]

    return npn ? npn.toString() : null
}
