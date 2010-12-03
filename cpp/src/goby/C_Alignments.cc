#include <string>
#include <iostream>
#include <pcrecpp.h>

#include "Reads.h"
#include "C_Alignments.h"
#include "MessageChunks.h"

using namespace std;

#ifdef C_WRITE_API_WRITE_ALIGNMENT_DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/**
 * C API to enable writing Goby compact-alignments in C.
 */
extern "C" {

	void gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(char *basename, CAlignmentsWriterHelper **writerHelperpp) {
	    gobyAlignments_openAlignmentsWriter(basename, GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK, writerHelperpp);
	}

	void gobyAlignments_openAlignmentsWriter(char *basename, unsigned number_of_entries_per_chunk, CAlignmentsWriterHelper **writerHelperpp) {
        debug(fprintf(stderr,"Writing alignment to basename, entries per chunk=%d\n", basename, number_of_entries_per_chunk));
        *writerHelperpp = new CAlignmentsWriterHelper;
        CAlignmentsWriterHelper *writerHelper = *writerHelperpp;
        string basenameStr(basename);
        
        writerHelper->alignmentWriter = new goby::AlignmentWriter(basenameStr, number_of_entries_per_chunk);
        writerHelper->tmhWriter = new goby::TooManyHitsWriter(basenameStr, 1);
        writerHelper->alignmentEntry = NULL;
        writerHelper->sequenceVariation = NULL;
        writerHelper->lastSeqVarReadIndex = -1;
        writerHelper->smallestQueryIndex = -1;
        writerHelper->largestQueryIndex = -1;
        writerHelper->numberOfAlignedReads = 0;
        writerHelper->qualityAdjustment = 0;
        writerHelper->samHelper = NULL;
	}

    int gobyAlignments_getQualityAdjustment(CAlignmentsWriterHelper *writerHelper) {
        return writerHelper->qualityAdjustment;
    }
    void gobyAlignments_setQualityAdjustment(CAlignmentsWriterHelper *writerHelper, int value) {
        writerHelper->qualityAdjustment = value;
    }
    void gobyAlignments_setAlignerName(CAlignmentsWriterHelper *writerHelper, char *value) {
        string valueStr(value);
        writerHelper->alignmentWriter->setAlignerName(valueStr);
    }
    void gobyAlignments_setAlignerVersion(CAlignmentsWriterHelper *writerHelper, char *value) {
        string valueStr(value);
        writerHelper->alignmentWriter->setAlignerVersion(valueStr);
    }
    void gobyAlignments_setSorted(CAlignmentsWriterHelper *writerHelper, int sorted /* bool */) {
        debug(fprintf(stderr,"gobyAlignments_setSorted=%d\n", sorted));
        writerHelper->alignmentWriter->setSorted(sorted == 0 ? false : true);
    }
    void gobyAlignments_setIndexed(CAlignmentsWriterHelper *writerHelper, int indexed /* bool */) {
        debug(fprintf(stderr,"gobyAlignments_setIndexed=%d\n", indexed));
        writerHelper->alignmentWriter->setIndexed(indexed == 0 ? false : true);
    }
    void gobyAlignments_addStatisticStr(CAlignmentsWriterHelper *writerHelper, const char *description, const char *value) {
        debug(fprintf(stderr,"gobyAlignments_addStatisticStr %s=%s\n", description, value));
        string descriptionStr(description);
        string valueStr(value);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, valueStr);
    }
    void gobyAlignments_addStatisticInt(CAlignmentsWriterHelper *writerHelper, const char *description, const int value) {
        debug(fprintf(stderr,"gobyAlignments_addStatisticInt %s=%d\n", description, value));
        string descriptionStr(description);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, value);
    }
    void gobyAlignments_addStatisticDouble(CAlignmentsWriterHelper *writerHelper, const char *description, const double value) {
        debug(fprintf(stderr,"gobyAlignments_addStatisticDouble %s=%f\n", description, value));
        string descriptionStr(description);
        writerHelper->alignmentWriter->addStatistic(descriptionStr, value);
    }

    /**
     * This should be called once for each target, targetIndex should start at 0
     * on the first call and increment by one on each following call.
     */
    void gobyAlignments_addTarget(CAlignmentsWriterHelper *writerHelper, const unsigned int targetIndex, const char *targetName, const unsigned int targetLength) {
        debug(fprintf(stderr,"gobyAlignments_addTargetIdentifier %s=%d\n", targetName, targetIndex));
        string targetNameStr(targetName);
        writerHelper->alignmentWriter->addTargetIdentifier(targetNameStr, targetIndex);
        writerHelper->alignmentWriter->addTargetLength(targetLength);
    }
    
    // get an empty alignment entry to populate
    void gobyAlignments_appendEntry(CAlignmentsWriterHelper *writerHelper) {
        debug(fprintf(stderr,"---------------------------------------\n"));
        debug(fprintf(stderr,"gobyAlignments_appendEntry\n"));
        writerHelper->alignmentEntry = writerHelper->alignmentWriter->appendEntry();
        writerHelper->sequenceVariation = NULL;
	    writerHelper->lastSeqVarReadIndex = -1;
    }
    void gobyAlEntry_setMultiplicity(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setMultiplicity=%d\n", value));
        writerHelper->alignmentEntry->set_multiplicity(value);
    }
    void gobyAlEntry_setQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryIndex=%d\n", value));
        if (writerHelper->smallestQueryIndex == -1) {
            writerHelper->smallestQueryIndex = value;
            writerHelper->largestQueryIndex = value;
        } else {
            writerHelper->smallestQueryIndex = min(value, writerHelper->smallestQueryIndex);
            writerHelper->largestQueryIndex = max(value, writerHelper->smallestQueryIndex);
        }
        writerHelper->alignmentEntry->set_query_index(value);
    }
    /**
     * Target index should be 0 based.
     */
    void gobyAlEntry_setTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setTargetIndex=%d\n", value));
        writerHelper->alignmentEntry->set_target_index(value);
    }
    void gobyAlEntry_setPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setPosition=%d\n", value));
        writerHelper->alignmentEntry->set_position(value);
    }
    void gobyAlEntry_setMatchingReverseStrand(CAlignmentsWriterHelper *writerHelper, int value /* bool */) {
        debug(fprintf(stderr,"gobyAlEntry_setMatchingReverseStrand=%d\n", value));
        writerHelper->alignmentEntry->set_matching_reverse_strand(value == 0 ? false : true);
    }
    void gobyAlEntry_setQueryPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryPosition=%d\n", value));
        writerHelper->alignmentEntry->set_query_position(value);
    }
    void gobyAlEntry_setScoreInt(CAlignmentsWriterHelper *writerHelper, int value) {
        float fValue = (float) value;
        debug(fprintf(stderr,"gobyAlEntry_setScore=%f\n", fValue));
        writerHelper->alignmentEntry->set_score(fValue);
    }
    void gobyAlEntry_setNumberOfMismatches(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setNumberOfMismatches=%d\n", value));
        writerHelper->alignmentEntry->set_number_of_mismatches(value);
    }
    void gobyAlEntry_setNumberOfIndels(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setNumberOfIndels=%d\n", value));
        writerHelper->alignmentEntry->set_number_of_indels(value);
    }
    void gobyAlEntry_setQueryAlignedLength(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryAlignedLength=%d\n", value));
        writerHelper->alignmentEntry->set_query_aligned_length(value);
    }
    void gobyAlEntry_setTargetAlignedLength(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setTargetAlignedLength=%d\n", value));
        writerHelper->alignmentEntry->set_target_aligned_length(value);
    }
    void gobyAlEntry_setQueryLength(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryLength=%d\n", value));
        writerHelper->alignmentEntry->set_query_length(value);
    }
    void gobyAlEntry_setMappingQuality(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->set_mapping_quality(value);
    }

    void startNewSequenceVariation(CAlignmentsWriterHelper *writerHelper, int readIndex) {
        debug(fprintf(stderr,"... startNewSequenceVariation\n"));
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
        bool hasQualChar = hasQualCharInt == 0 ? false : true;
        debug(fprintf(stderr,"gobyAlEntry_addSequenceVariation readIndex=%d ref=%c read=%c hasQualChar=%s\n", readIndex, refChar, readChar, hasQualChar ? "true" : "false"));
        if (writerHelper->sequenceVariation == NULL || writerHelper->lastSeqVarReadIndex == -1) {
            // New sequence variation
            startNewSequenceVariation(writerHelper, readIndex);
        } else if (writerHelper->lastSeqVarReadIndex + 1 == readIndex) {
            // Append to prev SequenceVar entry
            debug(fprintf(stderr,"... appending to previous seqVar\n"));
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
        debug(fprintf(stderr,"... sv->read_index=%d sv->position=%d from=%s to=%s\n",
            writerHelper->sequenceVariation->read_index(), writerHelper->sequenceVariation->position(),
            from->c_str(), to->c_str()));
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
    void gobyAlEntry_appendTooManyHits(CAlignmentsWriterHelper *writerHelper, unsigned int queryIndex, unsigned int alignedLength, int numberOfHits) {
        debug(fprintf(stderr,"gobyAlEntry_appendTooManyHits queryIndex=%d numberOfHits=%d alignedLength=%d\n", queryIndex, numberOfHits, alignedLength));
        writerHelper->tmhWriter->append(queryIndex, numberOfHits, alignedLength);
    }

    CSamHelper *samHelper_getResetSamHelper(CAlignmentsWriterHelper *writerHelper) {
        CSamHelper *samHelper = writerHelper->samHelper;
        if (samHelper == NULL) {
            writerHelper->samHelper = new CSamHelper;
            samHelper = writerHelper->samHelper;           
            samHelper->cpp_cigar = new string();
            samHelper->cpp_md = new string();
            samHelper->cpp_sourceQuery = new string();
            samHelper->cpp_sourceQual = new string();
            samHelper->cpp_query = new string();
            samHelper->cpp_qual = new string();
            samHelper->cpp_ref = new string();
            samHelper->minQualValue = writerHelper->qualityAdjustment;
        } else {
            samHelper->cpp_cigar->clear();
            samHelper->cpp_md->clear();
            samHelper->cpp_sourceQuery->clear();
            samHelper->cpp_sourceQual->clear();
            samHelper->cpp_query->clear();
            samHelper->cpp_qual->clear();
            samHelper->cpp_ref->clear();
        }
        samHelper->alignedLength = 0;
        samHelper->startPosition = 0;
        samHelper->numIndels = 0;
        samHelper->numMisMatches = 0;
        samHelper->score = 0;
        return samHelper;
    }

    void samHelper_addCigarItem(CSamHelper *samHelper, int length, char op) {
        std:stringstream lengthStream;
        lengthStream << length;
        (*samHelper->cpp_cigar) += lengthStream.str();
        (*samHelper->cpp_cigar) += op;
    }

    void samHelper_setMd(CSamHelper *samHelper, char *md) {
        if (md) {
            (*samHelper->cpp_md) += md;
        }
    }

    void samHelper_setQueryTranslate(CSamHelper *samHelper, char *reads, char *qual, int length, int reverseStrand) {
        if (!reverseStrand && qual) {
            (*samHelper->cpp_sourceQual) += qual;
        }
        int i;
        for (i = 0; i < length; i++) {
            if (reverseStrand) {
                (*samHelper->cpp_sourceQuery) += "TGCAN"[(int)(reads[length - 1 - i])];
                if (qual) {
                    (*samHelper->cpp_sourceQual) += (qual[length - 1 - i]);
                }
            } else {
                (*samHelper->cpp_sourceQuery) += "ACGTN"[(int)(reads[i])];
            }
        }
    }

    const char *samHelper_getCigarStr(CSamHelper *samHelper) {
        return samHelper->cpp_cigar->c_str();
    }

    void applyCigar(CSamHelper *samHelper) {
        pcrecpp::RE re("([0-9]+)([MID])");
        pcrecpp::StringPiece input(samHelper->cpp_cigar->c_str());
        int length;
        char op;
        int posInReads = 0;
        int i;
        while (re.Consume(&input, &length, &op)) {
            switch(op) {
                case 'M':
                    // Any mis-matches will be fixed in applyMd()
                    (*samHelper->cpp_query) += samHelper->cpp_sourceQuery->substr(posInReads, length);
                    if (samHelper->cpp_sourceQual->size() != 0) {
                        (*samHelper->cpp_qual) += samHelper->cpp_sourceQual->substr(posInReads, length);
                    }
                    (*samHelper->cpp_ref) += samHelper->cpp_sourceQuery->substr(posInReads, length);
                    posInReads += length;
                    break;
                case 'I':
                    (*samHelper->cpp_query) += samHelper->cpp_sourceQuery->substr(posInReads, length);
                    if (samHelper->cpp_sourceQual->size() != 0) {
                        (*samHelper->cpp_qual) += samHelper->cpp_sourceQual->substr(posInReads, length);
                    }
                    for (i = 0; i < length; i++) {
                        (*samHelper->cpp_ref) += '-';
                    }
                    posInReads += length;
                    break;
                case 'D':
                    for (i = 0; i < length; i++) {
                        (*samHelper->cpp_query) += '-';
                        if (samHelper->cpp_sourceQual->size() != 0) {
                            // Minimum qual
                            (*samHelper->cpp_qual) += ((char) samHelper->minQualValue); // min quality
                        }
                        (*samHelper->cpp_ref) += '?';  // provided in applyMd()
                    }
                    break;
            }
        }
    }

    void applyMd(CSamHelper *samHelper) {
        // My RE is simplified from the SAM spec but performs the same task
        // the major difference being mine would allow 5ACG where the current
        // spec would require 5A0C0G0 (which mine will still work with just fine).
        pcrecpp::RE re("([0-9]+|[ACGTN]|\\^[ACGTN]+)");
        pcrecpp::StringPiece input(samHelper->cpp_md->c_str());
        string mdPart;
        int position = 0;
        int i;
        while (re.Consume(&input, &mdPart)) {
            if (isdigit(mdPart[0])) {
                int length = atoi(mdPart.c_str());
                position += length;
            } else if (mdPart[0] == '^') {
                // Adjust the ref with these characters, ignoring the ^ character so start at 1
                for(i = 1; i < mdPart.size(); i++) {
                    (*samHelper->cpp_ref)[position++] = mdPart[i];
                }
            } else {
                // The regex should only allow a single character here, but we'll accept multiple
                for(i = 0; i < mdPart.size(); i++) {
                    (*samHelper->cpp_ref)[position++] = mdPart[i];
                }
            }
        }
    }

    void samHelper_constructRefAndQuery(CSamHelper *samHelper) {
        samHelper->cpp_query->clear();
        samHelper->cpp_qual->clear();
        samHelper->cpp_ref->clear();
        samHelper->alignedLength = 0;
        samHelper->startPosition = 0;
        samHelper->numIndels = 0;
        samHelper->numMisMatches = 0;
        samHelper->score = 0;

        applyCigar(samHelper);
        applyMd(samHelper);
        samHelper->alignedLength = samHelper->cpp_query->size();
        
        // Figure out start of alignment and alignment length, by observing mismatches at head and tail
        if (samHelper->cpp_query->size() != samHelper->cpp_ref->size()) {
            printf("ERROR! reconstructed reads and refs don't match in size!!\n");
        }
    }

    /**
     * Depends on samHelper_setQueryTranslate() being called.
     */
    const char *samHelper_sourceQuery(CSamHelper *samHelper) {
        return samHelper->cpp_sourceQuery->c_str();
    }
    /**
     * Depends on samHelper_setQueryTranslate() being called.
     */
    const char *samHelper_sourceQual(CSamHelper *samHelper) {
        if (samHelper->cpp_sourceQual->size() == 0) {
            return NULL;
        } else {
            return samHelper->cpp_sourceQual->c_str();
        }
    }
    /**
     * Depends on samHelper_constructRefAndQuery() being called.
     */
    const char *samHelper_constructedQuery(CSamHelper *samHelper) {
        return samHelper->cpp_query->c_str();
    }
    /**
     * Depends on samHelper_constructRefAndQuery() being called.
     */
    const char *samHelper_constructedQual(CSamHelper *samHelper) {
        if (samHelper->cpp_qual->size() == 0) {
            return NULL;
        } else {
            return samHelper->cpp_qual->c_str();
        }
    }

    /**
     * Depends on samHelper_constructRefAndQuery() being called.
     */
    const char *samHelper_constructedRef(CSamHelper *samHelper) {
        return samHelper->cpp_ref->c_str();
    }

	void gobyAlignments_finished(CAlignmentsWriterHelper *writerHelper, unsigned int numberOfReads) {
        debug(fprintf(stderr,"gobyAlignments_finished\n"));
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
            if (writerHelper->samHelper != NULL) {
                delete writerHelper->samHelper->cpp_cigar;
                delete writerHelper->samHelper->cpp_md;
                delete writerHelper->samHelper->cpp_sourceQuery;
                delete writerHelper->samHelper->cpp_sourceQual;
                delete writerHelper->samHelper->cpp_query;
                delete writerHelper->samHelper->cpp_qual;
                delete writerHelper->samHelper->cpp_ref;
                delete writerHelper->samHelper;
            }
            delete writerHelper;
        }
	}

    char* hitTypes[] = { "EXACT", "SUB", "INS", "DEL", "SPLICE", "TERMINAL" };

    void gobyAlignments_debugSequences(CAlignmentsWriterHelper *writerHelper, int hitType, char *refSequence, char *readSequence, int startPos) {
        debug(
            fprintf(stderr,":: type=%s\n", hitTypes[hitType]);
            string prefix;
            prefix.reserve(startPos);
            for (int i = 0; i < startPos; i++) {
                prefix += "_";
            }

            int suffixLength = writerHelper->alignmentEntry->query_length() - strlen(refSequence) - startPos;
            string suffix;
            suffix.reserve(suffixLength);
            for (int i = 0; i < suffixLength; i++) {
                suffix += "_";
            }

            fprintf(stderr,":: ref =%s%s%s\n", prefix.c_str(), refSequence, suffix.c_str());
            fprintf(stderr,":: read=%s%s%s\n", prefix.c_str(), readSequence, suffix.c_str());
        )
    }
}

