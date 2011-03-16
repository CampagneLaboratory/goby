/*
 * Helper structs to assist with read / writing the Goby compact formats using the C API.
 */

#ifndef C_ALIGNMENTS_H_
#define C_ALIGNMENTS_H_

#include "C_CompactHelpers.h"

#ifdef __cplusplus

extern "C" {
#endif
    void gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(char *basename, CAlignmentsWriterHelper **writerHelperpp);
    void gobyAlignments_openAlignmentsWriter(char *basename, unsigned int number_of_entries_per_chunk, CAlignmentsWriterHelper **writerHelperpp);

    /**
     * Provides a way to write data to a temporary file (or memory stream) for programs that,
     * for instance, write SAM data to a file. Use this to capture the SAM data for processing
     * by gobyAlignments_processSAM(...).
     */
    void gobyAlignments_openIntermediateOutputFiles(CAlignmentsWriterHelper *writerHelper, int openIgnoredOutputFile);
    FILE *gobyAlignments_intermediateOutputFileHandle(CAlignmentsWriterHelper *writerHelper);
    FILE *gobyAlignments_intermediateIgnoredOutputFileHandle(CAlignmentsWriterHelper *writerHelper);
    void gobyAlignments_intermediateOutputStartNew(CAlignmentsWriterHelper *writerHelper);
    void gobyAlignments_intermediateOutputFlush(CAlignmentsWriterHelper *writerHelper);
    char *gobyAlignments_intermediateOutputData(CAlignmentsWriterHelper *writerHelper);
    char *gobyAlignments_intermediateOutputIgnoredData(CAlignmentsWriterHelper *writerHelper);
    void gobyAlignments_closeIntermediateOutputFiles(CAlignmentsWriterHelper *writerHelper);

    void gobyAlignments_setAlignerName(CAlignmentsWriterHelper *writerHelper, char *value);
    void gobyAlignments_setAlignerVersion(CAlignmentsWriterHelper *writerHelper, char *value);
    int  gobyAlignments_getQualityAdjustment(CAlignmentsWriterHelper *writerHelper);
    void gobyAlignments_setQualityAdjustment(CAlignmentsWriterHelper *writerHelper, int value);
    void gobyAlignments_setSorted(CAlignmentsWriterHelper *writerHelper, int sorted /* bool */);
    void gobyAlignments_setIndexed(CAlignmentsWriterHelper *writerHelper, int indexed /* bool */);
    void gobyAlignments_setTargetLengths(CAlignmentsWriterHelper *writerHelper, const unsigned int* target_lengths);
    void gobyAlignments_addStatisticStr(CAlignmentsWriterHelper *writerHelper, const char *description, const char *value);
    void gobyAlignments_addStatisticInt(CAlignmentsWriterHelper *writerHelper, const char *description, const int value);
    void gobyAlignments_addStatisticDouble(CAlignmentsWriterHelper *writerHelper, const char *description, const double value);
    void gobyAlignments_addTarget(CAlignmentsWriterHelper *writerHelper, const unsigned int targetIndex, const char *targetName, unsigned int targetLength);
    unsigned gobyAlignments_addQueryIdentifier(CAlignmentsWriterHelper *writerHelper, const char *queryIdentifier);
    void gobyAlignments_addQueryIdentifierWithIndex(CAlignmentsWriterHelper *writerHelper, const char *queryIdentifier, unsigned int newQueryIndex);

    /** Create entries from from line/lines of SAM. */
    void gobyAlignments_processSAM(CAlignmentsWriterHelper *writerHelper, char *samLine, int npaths, int maxpaths);

    // get an empty alignment entry to populate
    void gobyAlignments_appendEntry(CAlignmentsWriterHelper *writerHelper);
    void gobyAlignments_debugSequences(CAlignmentsWriterHelper *writerHelper, int hitType, char *refSequence, char *readSequence, int startPos);
    void gobyAlEntry_setMultiplicity(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setMatchingReverseStrand(CAlignmentsWriterHelper *writerHelper, int value /* bool */);
    void gobyAlEntry_setQueryPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setScoreInt(CAlignmentsWriterHelper *writerHelper, int value);
    void gobyAlEntry_setNumberOfMismatches(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setNumberOfIndels(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setQueryAlignedLength(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setTargetAlignedLength(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setQueryLength(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setMappingQuality(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setFragmentIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setInsertSize(CAlignmentsWriterHelper *writerHelper, unsigned int value);

    // These are only used when dealing with a Query Pair
    void gobyAlEntry_setPairFlags(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setPairTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setPairPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setPairFragmentIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);

    // These are only used when dealing with a Splice
    void gobyAlEntry_setSplicedFlags(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setSplicedTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setSplicedFragmentIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value);
    void gobyAlEntry_setSplicedPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value);

    void gobyAlEntry_appendTooManyHits(CAlignmentsWriterHelper *writerHelper, unsigned int queryIndex, unsigned int alignedLength, int numberOfHits);
    void gobyAlEntry_addSequenceVariation(CAlignmentsWriterHelper *writerHelper, int readIndex, char refChar, char readChar, int hasQualCharInt /* bool */, char readQualChar);

    /**
     * Methods to assist reconstructing query and reference for SAM data.
     * Reconstructed query and ref strings will have '-' in the appropriate
     * sequence to denote insertions and deletions.
     */
    CSamHelper *samHelper_getResetSamHelper(CAlignmentsWriterHelper *writerHelper);
    void samHelper_addCigarItem(CSamHelper *samHelper, int length, char op);
    void samHelper_setCigar(CSamHelper *samHelper, const char *cigar);
    const char *samHelper_getCigarStr(CSamHelper *samHelper);
    void samHelper_setMd(CSamHelper *samHelper, const char *md);
    void samHelper_setQueryTranslate(CSamHelper *samHelper, char *reads, char *qual, unsigned int length, unsigned int reverseStrand);
    void samHelper_setQuery(CSamHelper *samHelper, const char *reads, const char *qual, unsigned int length, bool reverseStrand);
    void samHelper_constructRefAndQuery(CSamHelper *samHelper);
    const char *samHelper_sourceQuery(CSamHelper *samHelper);
    const char *samHelper_constructedQuery(CSamHelper *samHelper);
    const char *samHelper_sourceQual(CSamHelper *samHelper);
    const char *samHelper_constructedQual(CSamHelper *samHelper);
    const char *samHelper_constructedRef(CSamHelper *samHelper);

	void gobyAlignments_finished(CAlignmentsWriterHelper *alWriterHelper, unsigned int numberOfReads);
#ifdef __cplusplus
}
#endif


#endif /* C_ALIGNMENTS_H_ */
