#include <string>
#include <iostream>

#include "Reads.h"
#include "GsnapAlignments.h"
#include "GsnapStructs.h"
#include "MessageChunks.h"

using namespace std;

#undef DEBUG

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
	    writerHelper->tmhWriter = new goby::TooManyHitsWriter(basenameStr, 1);
	    writerHelper->alignmentEntry = NULL;
	    writerHelper->sequenceVariation = NULL;
	    writerHelper->lastSeqVarReadIndex = -1;
	    writerHelper->smallestQueryIndex = -1;
	    writerHelper->largestQueryIndex = -1;

        return writerHelper;
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
        writerHelper->alignmentEntry = writerHelper->alignmentWriter->appendEntry();
        writerHelper->sequenceVariation = NULL;
	    writerHelper->lastSeqVarReadIndex = -1;
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

        if (writerHelper->smallestQueryIndex == -1) {
            writerHelper->smallestQueryIndex = value;
            writerHelper->largestQueryIndex = value;
        } else {
            writerHelper->smallestQueryIndex = min(value, writerHelper->smallestQueryIndex);
            writerHelper->largestQueryIndex = max(value, writerHelper->smallestQueryIndex);
        }
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

    void startNewSequenceVariation(CAlignmentsWriterHelper *writerHelper, int readIndex) {
#ifdef DEBUG
        fprintf(stderr,"... startNewSequenceVariation\n");
#endif
        writerHelper->sequenceVariation = writerHelper->alignmentEntry->add_sequence_variations();
        writerHelper->sequenceVariation->set_read_index(readIndex);
        writerHelper->sequenceVariation->set_position(readIndex + 1);
        if (!writerHelper->alignmentEntry->matching_reverse_strand()) {
            // This will be correct for NOT-reverse... for reverse we'll set this later
            writerHelper->sequenceVariation->set_read_index(readIndex + 1);
        }
    }

    /**
     * This method will accumulate sequence variations. If > 1 are registered at readIndex positions next to
     * each other, they will be put into the SAME sequence variation. If a gap of great than one in readIndex
     * is found an additional SeqVar will be created.
     * This will also take "reverse" and "position" into account, which should already have been defined for
     * writerHelper->alignmentEntry.
     * @param readIndex coming in, readIndex is >0< based. When string position and readIndex into the SeqVar
     *        position and readIndex are >1< based.
     * For contiguous sequence variations, this assumes readIndex will increment by one each time.
     */
    void gobyAlEntry_addSequenceVariation(CAlignmentsWriterHelper *writerHelper, int readIndex, char refChar, char readChar, int hasQualCharInt /* bool */, char readQualChar) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_addSequenceVariation readIndex=%d ref=%c read=%c hasQualChar=%d\n", readIndex, refChar, readChar, hasQualCharInt);
#endif
        bool hasQualChar = intToBool(hasQualCharInt);
        if (writerHelper->sequenceVariation == NULL || writerHelper->lastSeqVarReadIndex == -1) {
            // New sequence variation
            startNewSequenceVariation(writerHelper, readIndex);
        } else if (writerHelper->lastSeqVarReadIndex + 1 == readIndex) {
            // Append to prev SequenceVar entry
#ifdef DEBUG
            fprintf(stderr,"... appending to previous seqVar\n");
#endif
        } else {
            // Not contiguous to previous SeqVar
            startNewSequenceVariation(writerHelper, readIndex);
        }
        if (writerHelper->alignmentEntry->matching_reverse_strand()) {
            // For reverse, update read_index as we accumulate characters for this SeqVar
            google::protobuf::uint32 readLength = writerHelper->alignmentEntry->query_length();
            writerHelper->sequenceVariation->set_read_index(readLength - readIndex);
        }
        writerHelper->lastSeqVarReadIndex = readIndex;
        string *from = writerHelper->sequenceVariation->mutable_from();
        string *to = writerHelper->sequenceVariation->mutable_to();
        (*from) += refChar;
        (*to) += readChar;

        if (hasQualChar) {
            string *toQuality = writerHelper->sequenceVariation->mutable_to_quality();
            (*toQuality) += readQualChar;
        }
    }

    /**
     * This method should be called once AFTER you call
     *    gobyAlEntry_setQueryIndex(...)
     *    gobyAlEntry_setQueryAlignedLength(...)
     */
    void gobyAlEntry_appendTooManyHits(CAlignmentsWriterHelper *writerHelper, int numberOfHits) {
        google::protobuf::uint32 queryIndex = writerHelper->alignmentEntry->query_index();
        google::protobuf::uint32 alignedLength = writerHelper->alignmentEntry->query_aligned_length();
#ifdef DEBUG
        fprintf(stderr,"gobyAlEntry_appendTooManyHits queryIndex=%d numberOfHits=%d alignedLength=%d\n", queryIndex, numberOfHits, alignedLength);
#endif
        writerHelper->tmhWriter->append(queryIndex, numberOfHits, alignedLength);
    }

	void gobyAlignments_finished(CAlignmentsWriterHelper *writerHelper, unsigned int numberOfReads) {
#ifdef DEBUG
        fprintf(stderr,"gobyAlignments_finished\n");
#endif
        if (writerHelper != NULL) {
            // Stats for the alignment
            writerHelper->alignmentWriter->setNumberOfQueries(numberOfReads);
            writerHelper->alignmentWriter->setSmallestSplitQueryIndex(writerHelper->smallestQueryIndex);
            writerHelper->alignmentWriter->setLargestSplitQueryIndex(writerHelper->largestQueryIndex);
            writerHelper->alignmentWriter->setNumberOfAlignedReads(writerHelper->numberOfAlignedReads);

            // Close and delete
            writerHelper->alignmentWriter->close();
            delete writerHelper->alignmentWriter;
            writerHelper->tmhWriter->write();
            delete writerHelper->tmhWriter;
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
