/*
 * GsnapReads.h
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */

#ifndef GSNAPALIGNMENTS_H_
#define GSNAPALIGNMENTS_H_

#include "CompactFormatsHelpers.h"

#ifdef __cplusplus
#include "GsnapStructs.h"
bool intToBool(int value);   // This function is NOT extern'd to C
extern "C" {
#endif
	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(char *basename);
	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename, unsigned int number_of_entries_per_chunk);

    void gobyAlignments_setNumberOfQueries(CAlignmentsWriterHelper *writerHelper, unsigned number_of_queries);
    void gobyAlignments_setNumberOfTargets(CAlignmentsWriterHelper *writerHelper, unsigned number_of_targets);
    void gobyAlignments_setNumberOfAlignedReads(CAlignmentsWriterHelper *writerHelper, unsigned number_of_aligned_reads);
    void gobyAlignments_setConstantQuerylength(CAlignmentsWriterHelper *writerHelper, unsigned constant_query_length);
    void gobyAlignments_setSmallestSplitQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned smallest_split_query_index);
    void gobyAlignments_setLargestSplitQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned largest_split_query_index);
    void gobyAlignments_setSorted(CAlignmentsWriterHelper *writerHelper, int sorted /* bool */);
    void gobyAlignments_setIndexed(CAlignmentsWriterHelper *writerHelper, int indexed /* bool */);
    void gobyAlignments_setTargetLengths(CAlignmentsWriterHelper *writerHelper, const unsigned int* target_lengths);
    void gobyAlignments_addStatisticStr(CAlignmentsWriterHelper *writerHelper, const char *description, const char *value);
    void gobyAlignments_addStatisticInt(CAlignmentsWriterHelper *writerHelper, const char *description, const int value);
    void gobyAlignments_addStatisticDouble(CAlignmentsWriterHelper *writerHelper, const char *description, const double value);
    void gobyAlignments_addTargetIdentifier(CAlignmentsWriterHelper *writerHelper, const UINT4 targetIndex, const char *targetName);

    // get an empty alignment entry to populate
    void gobyAlignments_appendEntry(CAlignmentsWriterHelper *writerHelper); // ::goby::AlignmentEntry*
    void gobyAlEntry_setMultiplicity(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setQueryIndex(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setTargetIndex(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setPosition(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setMatchingReverseStrand(CAlignmentsWriterHelper *writerHelper, int value /* bool */);
    void gobyAlEntry_setQueryPosition(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setScoreInt(CAlignmentsWriterHelper *writerHelper, int value);
    void gobyAlEntry_setNumberOfMismatches(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setNumberOfIndels(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setQueryAlignedLength(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setTargetAlignedLength(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_setQueryLength(CAlignmentsWriterHelper *writerHelper, UINT4 value);
    void gobyAlEntry_addSequenceVariation(CAlignmentsWriterHelper *writerHelper, int readIndex, char refChar, char readChar, int hasQualCharInt /* bool */, char readQualChar);

	void gobyAlignments_finished(CAlignmentsWriterHelper *alWriterHelper);
#ifdef __cplusplus
}
#endif


#endif /* GSNAPALIGNMENTS_H_ */
