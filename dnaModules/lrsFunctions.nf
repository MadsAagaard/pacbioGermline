#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * lrsFunctions.nf — shared helper functions
 * =============================================================================
 *
 * Plain functions, included by name:
 *
 *     include { sexFromGender; parseUbamName } from './modules/lrsFunctions.nf'
 *
 * Only two things live here, and both earn it by being used in more than one
 * place. Anything used by a single script stays in that script — a shared file
 * that collects one-off helpers is harder to reason about than the duplication
 * it replaces.
 *
 * WHAT BELONGS HERE
 *   - pure functions parsing or validating input from outside the pipeline
 *   - logic that MUST behave identically everywhere (sex validation above all)
 *
 * WHAT DOES NOT
 *   - anything creating or operating on a Channel
 *   - publishDir / output-path logic (that is params.outBase, in the config)
 *   - process definitions
 *
 * A NOTE ON CONCURRENCY
 * These are called from inside map/each closures. They are pure: no state is
 * held between calls and nothing outside is mutated. Do not add a script-level
 * variable to this file — closures running on the dataflow thread pool would
 * write to it concurrently, which is the bug that clobbered date2 in the old
 * main.nf.
 *
 * A NOTE ON STRINGS
 * Everything returned goes through .toString(). Values arriving from splitCsv
 * and GStrings are not always java.lang.String, and a GString has a different
 * hashCode than the equal String — so a GString used as a map key silently
 * fails lookup later. Plain Strings only.
 * =============================================================================
 */


/**
 * Normalise a raw gender field to 'male' or 'female'.
 *
 * Fail-closed by design: anything unrecognised throws, and there is no default.
 * A missing or malformed gender is a data-entry error, not a state the pipeline
 * should have a policy for — TRGT --karyotype and Sawfish --expected-cn both
 * depend on it being right.
 *
 * This is the reason this file exists. Duplicating it is cheap right up until
 * someone relaxes one copy to accept a blank field for a test run, and the
 * other copy keeps everyone honest while that one silently genotypes XX.
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
 * Accepts BOTH layouts, so no preprocessing step is needed:
 *
 *   underscore-joined (as emitted):
 *     0000103212_113624121888_78_SL-NGC-SJAELDNE_K_T_113782522715_3
 *
 *   tab-exploded (tr '_' '\t', for manual inspection):
 *     0000103212<TAB>113624121888<TAB>78<TAB>SL-NGC-SJAELDNE<TAB>K<TAB>...
 *
 * Detection is unambiguous: the joined blob contains no tabs and therefore
 * splits to exactly one column, while the exploded form never yields one.
 *
 *   [0] rekv   [1] npn      [2] material  [3] testlist
 *   [4] gender [5] proband  [6] intRef    [7] expectedCount (recent addition)
 *
 * 7 and 8 fields are both accepted — expectedCount was added late and older
 * concatenated sheets will not have it. Fewer than 7 throws: a short line means
 * a malformed record, and guessing which field is missing is how the wrong
 * value ends up in `gender`.
 */
def parseMetaLine(line) {

    def cols = line.toString().trim().split('\t') as List
    def f    = (cols.size() == 1) ? (cols[0].tokenize('_')) : cols

    if (f.size() < 7 || f.size() > 8) {
        throw new IllegalArgumentException(
            "Metadata line '${line}' has ${f.size()} fields, expected 7 or 8 " +
            "(rekv_npn_material_testlist_gender_proband_intRef[_expectedCount]).")
    }

    def npn = f[1].toString().trim()
    if (!npn) throw new IllegalArgumentException("Metadata line '${line}' has an empty NPN (field 2).")

    return [
        rekv          : f[0].toString().trim(),
        id            : npn,
        material      : f[2].toString().trim(),
        testlist      : f[3].toString().trim(),
        gender        : f[4].toString().trim(),
        sex           : sexFromGender(f[4], npn),
        proband       : f[5].toString().trim(),
        intRef        : f[6].toString().trim(),
        expectedCount : (f.size() == 8) ? f[7].toString().trim() : null,
    ]
}


/**
 * Parse a PacBio uBAM filename into its components.
 *
 *   113738820048_78_NGC-ARVELIG-KRFT_K.m84313_260317_105642_s1.hifi_reads.bc2075.bam
 *   \__________/ \/ \_____________/ \/  \_______________/     \________/ \____/
 *      npn      mat     testlist   gender    instrument+run     readset  barcode
 *
 * NOTE the offsets differ from the samplesheet: the uBAM name blob has NO rekv
 * field, so npn is [0] and gender is [3], against [1] and [4] in the sheet.
 *
 * ONLY npn AND gender MAY EVER BE CROSS-CHECKED.
 * Not testlist, not intRef, not rekv. Those are properties of a *referral*,
 * not of a sample: the same NPN can be re-referred later under a different
 * testlist, join a different family, or arrive on a new requisition, and the
 * uBAM filename is frozen at sequencing time while the samplesheet is not. A
 * check on any of them would fire on samples that are entirely correct.
 * (The 'SL-' prefix present in sheet testlists and absent from uBAM testlists
 * is a second, independent reason — but the mutability is the real one.)
 *
 * npn and gender are safe precisely because they are immutable properties of
 * the individual.
 *
 * Returns null rather than throwing, so callers can branch on unparseable
 * filenames (UNASSIGNED, ad-hoc test data) instead of the run dying on one
 * stray file. Callers MUST handle null.
 *
 * This parser exists twice in main.nf today — in the samplesheet and
 * no-samplesheet branches — with subtly different field handling, and a third
 * time in backgroundPerSample.nf. One copy is enough.
 */
def parseUbamName(bamPath) {
    def base = bamPath.toString().tokenize('/').last().replaceFirst(/\.bam$/, '')
    def dots = base.tokenize('.')

    if (dots.size() < 3) return null

    def nameField = dots[0].toString()
    def pacbioID  = dots[1].toString()
    def readSet   = dots[2].toString()
    def barcode   = (dots.size() > 3) ? dots[3].toString() : null

    def nameParts = nameField.tokenize('_')
    if (!nameParts) return null

    def runParts = pacbioID.tokenize('_')

    return [
        id         : nameParts[0].toString(),
        material   : (nameParts.size() > 1) ? nameParts[1].toString() : null,
        testlist   : (nameParts.size() > 2) ? nameParts[2].toString() : null,
        genderFile : (nameParts.size() > 3) ? nameParts[3].toString() : null,
        instrument : (runParts.size()  > 0) ? runParts[0].toString()  : null,
        rundate    : (runParts.size()  > 1) ? runParts[1].toString()  : null,
        readSet    : readSet,
        barcode    : barcode,
    ]
}
