#include <string>
#include <iostream>

#include "Reads.h"
#include "GsnapAlignments.h"
#include "GsnapStructs.h"
#include "MessageChunks.h"

using namespace std;

/**
 * This class is a C interface so Gsnap can write Goby compact-alignments.
 */
extern "C" {

	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(char *basename) {
	    return gobyAlignments_openAlignmentsWriter(basename, GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK);
	}

	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename, unsigned number_of_entries_per_chunk) {
		CAlignmentsWriterHelper *writerHelper = new CAlignmentsWriterHelper;
        string basenameStr(basename);

        writerHelper->alignmentWriter = new goby::AlignmentWriter(basenameStr, number_of_entries_per_chunk);
	    writerHelper->alignmentEntry = NULL;
	    writerHelper->sequenceVariation = NULL;
		writerHelper->numRead = 0;

        return writerHelper;
	}

    void gobyAlignments_setNumberOfQueries(CAlignmentsWriterHelper *writerHelper, unsigned number_of_queries) {
        writerHelper->alignmentWriter->setNumberOfQueries(number_of_queries);
    }
    void gobyAlignments_setNumberOfTargets(CAlignmentsWriterHelper *writerHelper, unsigned number_of_targets) {
        writerHelper->alignmentWriter->setNumberOfTargets(number_of_targets);
    }
    void gobyAlignments_setNumberOfAlignedReads(CAlignmentsWriterHelper *writerHelper, unsigned number_of_aligned_reads) {
        writerHelper->alignmentWriter->setNumberOfAlignedReads(number_of_aligned_reads);
    }
    void gobyAlignments_setConstantQuerylength(CAlignmentsWriterHelper *writerHelper, unsigned constant_query_length) {
        writerHelper->alignmentWriter->setConstantQuerylength(constant_query_length);
    }
    void gobyAlignments_setSmallestSplitQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned smallest_split_query_index) {
        writerHelper->alignmentWriter->setSmallestSplitQueryIndex(smallest_split_query_index);
    }
    void gobyAlignments_setLargestSplitQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned largest_split_query_index) {
        writerHelper->alignmentWriter->setLargestSplitQueryIndex(largest_split_query_index);
    }
    void gobyAlignments_setSorted(CAlignmentsWriterHelper *writerHelper, int sorted /* bool */) {
        writerHelper->alignmentWriter->setSorted(intToBool(sorted));
    }
    void gobyAlignments_setIndexed(CAlignmentsWriterHelper *writerHelper, int indexed /* bool */) {
        writerHelper->alignmentWriter->setIndexed(intToBool(indexed));
    }
    void gobyAlignments_setTargetLengths(CAlignmentsWriterHelper *writerHelper, const unsigned int* target_lengths) {
        writerHelper->alignmentWriter->setTargetLengths(target_lengths);
    }
    void gobyAlignments_addStatisticStr(CAlignmentsWriterHelper *writerHelper, const char *description, const char *value) {
        string descriptionStr(description);
        string valueStr(value);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, valueStr);
    }
    void gobyAlignments_addStatisticInt(CAlignmentsWriterHelper *writerHelper, const char *description, const int value) {
        string descriptionStr(description);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, value);
    }
    void gobyAlignments_addStatisticDouble(CAlignmentsWriterHelper *writerHelper, const char *description, const double value) {
        string descriptionStr(description);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, value);
    }

    // get an empty alignment entry to populate
    void gobyAlignments_appendEntry(CAlignmentsWriterHelper *writerHelper) {
        writerHelper->alignmentEntry = writerHelper->alignmentWriter->appendEntry();
    }
    void gobyAlEntry_setMultiplicity(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_multiplicity(value);
    }
    void gobyAlEntry_setQueryIndex(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_query_index(value);
    }
    void gobyAlEntry_setTargetIndex(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_target_index(value);
    }
    void gobyAlEntry_setPosition(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_position(value);
    }
    void gobyAlEntry_setMatchingReverseStrand(CAlignmentsWriterHelper *writerHelper, int value /* bool */) {
        writerHelper->alignmentEntry->set_matching_reverse_strand(intToBool(value));
    }
    void gobyAlEntry_setQueryPosition(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_query_position(value);
    }
    void gobyAlEntry_setScore(CAlignmentsWriterHelper *writerHelper, float value) {
        writerHelper->alignmentEntry->set_score(value);
    }
    void gobyAlEntry_setNumberOfMismatches(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_number_of_mismatches(value);
    }
    void gobyAlEntry_setNumberOfIndels(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_number_of_indels(value);
    }
    void gobyAlEntry_setQueryAlignedLength(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_query_aligned_length(value);
    }
    void gobyAlEntry_setTargetAlignedLength(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_target_aligned_length(value);
    }
    void gobyAlEntry_setQueryLength(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->alignmentEntry->set_query_length(value);
    }

    void gobyAlEntry_addSequenceVariations(CAlignmentsWriterHelper *writerHelper) {
        writerHelper->sequenceVariation = writerHelper->alignmentEntry->add_sequence_variations();
    }
    void gobyAlSeqVar_setFrom(CAlignmentsWriterHelper *writerHelper, const char* value) {
        writerHelper->sequenceVariation->set_from(value);
    }
    void gobyAlSeqVar_setTo(CAlignmentsWriterHelper *writerHelper, const char* value) {
        writerHelper->sequenceVariation->set_to(value);
    }
    void gobyAlSeqVar_setPosition(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->sequenceVariation->set_position(value);
    }
    void gobyAlSeqVar_setReadIndex(CAlignmentsWriterHelper *writerHelper, uint32_t value) {
        writerHelper->sequenceVariation->set_read_index(value);
    }
    void gobyAlSeqVar_setToQuality(CAlignmentsWriterHelper *writerHelper, const char* value) {
        writerHelper->sequenceVariation->set_to_quality(value);
    }

	void gobyAlignments_finished(CAlignmentsWriterHelper *writerHelper) {
        if (writerHelper != NULL) {
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
