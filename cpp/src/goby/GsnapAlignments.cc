#include <string>
#include <iostream>

#include "Reads.h"
#include "GsnapAlignments.h"
#include "GsnapStructs.h"
#include "MessageChunks.h"

using namespace std;

#define DEBUG

/**
 * This class is a C interface so Gsnap can write Goby compact-alignments.
 */
extern "C" {

	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(char *basename) {
	    return gobyAlignments_openAlignmentsWriter(basename, GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK);
	}

	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename, unsigned number_of_entries_per_chunk) {
#ifdef DEBUG
        fprintf(stderr,"Writing alignment to basename, entries per chunk=%d\n", basename, number_of_entries_per_chunk);
#endif
		CAlignmentsWriterHelper *writerHelper = new CAlignmentsWriterHelper;
        string basenameStr(basename);

        writerHelper->alignmentWriter = new goby::AlignmentWriter(basenameStr, number_of_entries_per_chunk);
	    writerHelper->alignmentEntry = NULL;
	    writerHelper->sequenceVariation = NULL;
		writerHelper->numWritten = 0;

        return writerHelper;
	}

    void gobyAlignments_setNumberOfQueries(CAlignmentsWriterHelper *writerHelper, unsigned number_of_queries) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setNumberOfQueries=%d\n", number_of_queries);
#endif
        writerHelper->alignmentWriter->setNumberOfQueries(number_of_queries);
    }
    void gobyAlignments_setNumberOfTargets(CAlignmentsWriterHelper *writerHelper, unsigned number_of_targets) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setNumberOfTargets=%d\n", number_of_targets);
#endif
        writerHelper->alignmentWriter->setNumberOfTargets(number_of_targets);
    }
    void gobyAlignments_setNumberOfAlignedReads(CAlignmentsWriterHelper *writerHelper, unsigned number_of_aligned_reads) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setNumberOfAlignedReads=%d\n", number_of_aligned_reads);
#endif
        writerHelper->alignmentWriter->setNumberOfAlignedReads(number_of_aligned_reads);
    }
    void gobyAlignments_setConstantQuerylength(CAlignmentsWriterHelper *writerHelper, unsigned constant_query_length) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setConstantQuerylength=%d\n", constant_query_length);
#endif
        writerHelper->alignmentWriter->setConstantQuerylength(constant_query_length);
    }
    void gobyAlignments_setSmallestSplitQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned smallest_split_query_index) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setSmallestSplitQueryIndex=%d\n", smallest_split_query_index);
#endif
        writerHelper->alignmentWriter->setSmallestSplitQueryIndex(smallest_split_query_index);
    }
    void gobyAlignments_setLargestSplitQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned largest_split_query_index) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setLargestSplitQueryIndex=%d\n", largest_split_query_index);
#endif
        writerHelper->alignmentWriter->setLargestSplitQueryIndex(largest_split_query_index);
    }
    void gobyAlignments_setSorted(CAlignmentsWriterHelper *writerHelper, int sorted /* bool */) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setSorted=%d\n", sorted);
#endif
        writerHelper->alignmentWriter->setSorted(intToBool(sorted));
    }
    void gobyAlignments_setIndexed(CAlignmentsWriterHelper *writerHelper, int indexed /* bool */) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setIndexed=%d\n", indexed);
#endif
        writerHelper->alignmentWriter->setIndexed(intToBool(indexed));
    }
    void gobyAlignments_setTargetLengths(CAlignmentsWriterHelper *writerHelper, const unsigned int* target_lengths) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_setTargetLengths=...\n");
#endif
        writerHelper->alignmentWriter->setTargetLengths(target_lengths);
    }
    void gobyAlignments_addStatisticStr(CAlignmentsWriterHelper *writerHelper, const char *description, const char *value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_addStatisticStr %s=%s\n", description, value);
#endif
        string descriptionStr(description);
        string valueStr(value);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, valueStr);
    }
    void gobyAlignments_addStatisticInt(CAlignmentsWriterHelper *writerHelper, const char *description, const int value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_addStatisticInt %s=%d\n", description, value);
#endif
        string descriptionStr(description);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, value);
    }
    void gobyAlignments_addStatisticDouble(CAlignmentsWriterHelper *writerHelper, const char *description, const double value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_addStatisticDouble %s=%f\n", description, value);
#endif
        string descriptionStr(description);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, value);
    }
    void gobyAlignments_addTargetIdentifier(CAlignmentsWriterHelper *writerHelper, const UINT4 targetIndex, const char *targetName) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_addTargetIdentifier %s=%d\n", targetName, targetIndex);
#endif
        writerHelper->alignmentWriter->addTargetIdentifier(targetName, targetIndex);
    }

    // get an empty alignment entry to populate
    void gobyAlignments_appendEntry(CAlignmentsWriterHelper *writerHelper) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_appendEntry\n");
#endif
        writerHelper->numWritten++;
        writerHelper->alignmentEntry = writerHelper->alignmentWriter->appendEntry();
    }
    void gobyAlEntry_setMultiplicity(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setMultiplicity=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_multiplicity(value);
    }
    void gobyAlEntry_setQueryIndex(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setQueryIndex=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_query_index(value);
    }
    void gobyAlEntry_setTargetIndex(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setTargetIndex=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_target_index(value);
    }
    void gobyAlEntry_setPosition(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setPosition=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_position(value);
    }
    void gobyAlEntry_setMatchingReverseStrand(CAlignmentsWriterHelper *writerHelper, int value /* bool */) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setMatchingReverseStrand=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_matching_reverse_strand(intToBool(value));
    }
    void gobyAlEntry_setQueryPosition(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setQueryPosition=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_query_position(value);
    }
    void gobyAlEntry_setScoreInt(CAlignmentsWriterHelper *writerHelper, int value) {
        float fValue = (float) value;
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setScore=%d,%f\n", value, fValue);
#endif
        writerHelper->alignmentEntry->set_score(fValue);
    }
    void gobyAlEntry_setNumberOfMismatches(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setNumberOfMismatches=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_number_of_mismatches(value);
    }
    void gobyAlEntry_setNumberOfIndels(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setNumberOfIndels=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_number_of_indels(value);
    }
    void gobyAlEntry_setQueryAlignedLength(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setQueryAlignedLength=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_query_aligned_length(value);
    }
    void gobyAlEntry_setTargetAlignedLength(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setTargetAlignedLength=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_target_aligned_length(value);
    }
    void gobyAlEntry_setQueryLength(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_setQueryLength=%d\n", value);
#endif
        writerHelper->alignmentEntry->set_query_length(value);
    }

    void gobyAlEntry_addSequenceVariations(CAlignmentsWriterHelper *writerHelper) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_addSequenceVariations\n");
#endif
        writerHelper->sequenceVariation = writerHelper->alignmentEntry->add_sequence_variations();
    }
    void gobyAlSeqVar_setFrom(CAlignmentsWriterHelper *writerHelper, const char* value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlSeqVar_setFrom=%s\n", value);
#endif
        writerHelper->sequenceVariation->set_from(value);
    }
    void gobyAlSeqVar_setTo(CAlignmentsWriterHelper *writerHelper, const char* value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlSeqVar_setTo=%s\n", value);
#endif
        writerHelper->sequenceVariation->set_to(value);
    }
    void gobyAlSeqVar_setPosition(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlSeqVar_setPosition=%d\n", value);
#endif
        writerHelper->sequenceVariation->set_position(value);
    }
    void gobyAlSeqVar_setReadIndex(CAlignmentsWriterHelper *writerHelper, UINT4 value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlSeqVar_setReadIndex=%d\n", value);
#endif
        writerHelper->sequenceVariation->set_read_index(value);
    }
    void gobyAlSeqVar_setToQuality(CAlignmentsWriterHelper *writerHelper, const char* value) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlSeqVar_setToQuality=%s\n", value);
#endif
        writerHelper->sequenceVariation->set_to_quality(value);
    }

	void gobyAlignments_finished(CAlignmentsWriterHelper *writerHelper) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_finished\n");
#endif
        if (writerHelper != NULL) {
            // TODO: Write the stats ... # of entries, etc.

            writerHelper->alignmentWriter->close();
            delete writerHelper->alignmentWriter;
            delete writerHelper;
        }
	}
}

/** Convert int to bool. NOT extern'd to C. **/
bool intToBool(int value) {
    if (value == 0) {
        return false;
    } else {
        return true;
    }
}
