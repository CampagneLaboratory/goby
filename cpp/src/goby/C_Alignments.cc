#include <string>
#include <iostream>
#include <pcrecpp.h>
#include <stdio.h>

#include "Reads.h"
#include "C_Alignments.h"
#include "MessageChunks.h"
#include "SamFlags.h"

using namespace std;

#undef C_WRITE_API_WRITE_ALIGNMENT_DEBUG
#ifdef C_WRITE_API_WRITE_ALIGNMENT_DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef __GNUC__
#define INTERMEDIATE_OUTPUT_VIA_OPEN_MEMSTREAM
#endif


// TODO: When reading from COMPACT-READS, one should call
// TODO: addQueryIdentifierWithInt() but ONLY if there is a text
// TODO: identifier for the read in question, otherwise just
// TODO: use the actual queryIndex that comes from CompactReads
// TODO: as the query index and don't bother with addQueryIdentifier*
// TODO: the addQueryIdentifier() should only be used when reading
// TODO: from FA/FQ and what you have is a real string identifier.

/**
 * C API to enable writing Goby compact-alignments in C.
 */
extern "C" {

	void gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(char *basename, CAlignmentsWriterHelper **writerHelperpp) {
	    gobyAlignments_openAlignmentsWriter(basename, GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK, writerHelperpp);
	}

	void gobyAlignments_openAlignmentsWriter(char *basename, unsigned int number_of_entries_per_chunk, CAlignmentsWriterHelper **writerHelperpp) {
        debug(fprintf(stderr,"Writing alignment to basename, entries per chunk=%u\n", basename, number_of_entries_per_chunk));
        *writerHelperpp = new CAlignmentsWriterHelper;
        CAlignmentsWriterHelper *writerHelper = *writerHelperpp;
        string basenameStr(basename);
        
        writerHelper->alignmentWriter = new goby::AlignmentWriter(basenameStr, number_of_entries_per_chunk);
        writerHelper->tmhWriter = new goby::TooManyHitsWriter(basenameStr, 1);
        writerHelper->alignmentEntry = NULL;
        writerHelper->sequenceVariation = NULL;
        writerHelper->lastSeqVarReadIndex = -1;
        writerHelper->lastSeqVarReadChar = '\0';
        writerHelper->lastSeqVarRefChar = '\0';
        writerHelper->smallestQueryIndex = -1;
        writerHelper->largestQueryIndex = -1;
        writerHelper->numberOfAlignedReads = 0;
        writerHelper->qualityAdjustment = 0;
        writerHelper->samHelper = NULL;
        writerHelper->intermediateOutputFile = NULL;
        writerHelper->intermediateOutputBuffer = NULL;
        writerHelper->intermediateOutputBufferSize = 0;
        writerHelper->intermediateIgnoredOutputFile = NULL;
        writerHelper->intermediateIgnoredOutputBuffer = NULL;
        writerHelper->intermediateIgnoredOutputBufferSize = 0;
        writerHelper->alignerToGobyTargetIndexMap = NULL;
	}

    void gobyAlignments_openIntermediateOutputFiles(CAlignmentsWriterHelper *writerHelper, int openIgnoredOutputFile) {
        gobyAlignments_closeIntermediateOutputFiles(writerHelper);
#ifdef INTERMEDIATE_OUTPUT_VIA_OPEN_MEMSTREAM
        fprintf(stderr, "Opening intermediate output via open_memstream()\n");
#else
        fprintf(stderr, "Opening intermediate output via temporary files, slower but non-gcc compatible.\n");
#endif
#ifdef INTERMEDIATE_OUTPUT_VIA_OPEN_MEMSTREAM
        writerHelper->intermediateOutputFile = open_memstream(&writerHelper->intermediateOutputBuffer, &writerHelper->intermediateOutputBufferSize);
        if (openIgnoredOutputFile) {
            writerHelper->intermediateIgnoredOutputFile = open_memstream(&writerHelper->intermediateIgnoredOutputBuffer, &writerHelper->intermediateIgnoredOutputBufferSize);
        }
#else
        writerHelper->intermediateOutputFile = tmpfile();
        if (openIgnoredOutputFile) {
            writerHelper->intermediateIgnoredOutputFile = tmpfile();
        }
#endif
    }

    FILE *gobyAlignments_intermediateOutputFileHandle(CAlignmentsWriterHelper *writerHelper) {
        return writerHelper->intermediateOutputFile;
    }
    FILE *gobyAlignments_intermediateIgnoredOutputFileHandle(CAlignmentsWriterHelper *writerHelper) {
        return writerHelper->intermediateIgnoredOutputFile;
    }

    /**
     * Start a new section of output.
     */
	void gobyAlignments_intermediateOutputStartNew(CAlignmentsWriterHelper *writerHelper) {
        if (writerHelper->intermediateOutputFile) {
            rewind(writerHelper->intermediateOutputFile);
        }
        if (writerHelper->intermediateIgnoredOutputFile) {
            rewind(writerHelper->intermediateIgnoredOutputFile);
        }
	}

    char *gobyAlignments_intermediateOutputData(CAlignmentsWriterHelper *writerHelper) {
        return writerHelper->intermediateOutputBuffer;
    }

    char *gobyAlignments_intermediateOutputIgnoredData(CAlignmentsWriterHelper *writerHelper) {
        return writerHelper->intermediateIgnoredOutputBuffer;
    }

    /**
     * A section of output is complete. Populate what was written to the files into
     * writerHelper->intermediateOutputBuffer / writerHelper->intermediateIgnoredOutputBuffer.
     * The output should be processed. After you're done processing the output,
     * when you're ready to start writing new output call gobyAlignments_intermediateOutputStartNew.
     */
	void gobyAlignments_intermediateOutputFlush(CAlignmentsWriterHelper *writerHelper) {
	    if (writerHelper->intermediateOutputFile) {
            fflush(writerHelper->intermediateOutputFile);
        }
        if (writerHelper->intermediateIgnoredOutputFile) {
            fflush(writerHelper->intermediateIgnoredOutputFile);
        }
#ifndef INTERMEDIATE_OUTPUT_VIA_OPEN_MEMSTREAM
        size_t fileSize;
        if (writerHelper->intermediateOutputFile) {
            fileSize = ftell(writerHelper->intermediateOutputFile);
            if (fileSize + 1 > writerHelper->intermediateOutputBufferSize) {
                writerHelper->intermediateOutputBufferSize = fileSize + 1;
                writerHelper->intermediateOutputBuffer = (char *) realloc(writerHelper->intermediateOutputBuffer, writerHelper->intermediateOutputBufferSize);
            }
            rewind(writerHelper->intermediateOutputFile);
            fread(writerHelper->intermediateOutputBuffer, 1, fileSize, writerHelper->intermediateOutputFile);
            writerHelper->intermediateOutputBuffer[fileSize] = '\0';
        }
        if (writerHelper->intermediateIgnoredOutputFile) {
            fileSize = ftell(writerHelper->intermediateIgnoredOutputFile);
            if (fileSize + 1 > writerHelper->intermediateIgnoredOutputBufferSize) {
                writerHelper->intermediateIgnoredOutputBufferSize = fileSize + 1;
                writerHelper->intermediateIgnoredOutputBuffer = (char *) realloc(writerHelper->intermediateIgnoredOutputBuffer, writerHelper->intermediateIgnoredOutputBufferSize);
            }
            rewind(writerHelper->intermediateIgnoredOutputFile);
            fread(writerHelper->intermediateIgnoredOutputBuffer, 1, fileSize, writerHelper->intermediateIgnoredOutputFile);
            writerHelper->intermediateIgnoredOutputBuffer[fileSize] = '\0';
        }
#endif
    }

	void gobyAlignments_closeIntermediateOutputFiles(CAlignmentsWriterHelper *writerHelper) {
	    if (writerHelper->intermediateOutputFile || writerHelper->intermediateIgnoredOutputFile) {
#ifdef INTERMEDIATE_OUTPUT_VIA_OPEN_MEMSTREAM
            fprintf(stderr, "Closing intermediate output via open_memstream()\n");
#else
            fprintf(stderr, "Closing intermediate output via temporary files\n");
#endif
        }
        if (writerHelper->intermediateOutputFile) {
            fclose(writerHelper->intermediateOutputFile);
            writerHelper->intermediateOutputFile = NULL;
            free(writerHelper->intermediateOutputBuffer);
            writerHelper->intermediateOutputBuffer = NULL;
            writerHelper->intermediateOutputBufferSize = 0;
        }
        if (writerHelper->intermediateIgnoredOutputFile) {
            fclose(writerHelper->intermediateIgnoredOutputFile);
            writerHelper->intermediateIgnoredOutputFile = NULL;
            free(writerHelper->intermediateIgnoredOutputBuffer);
            writerHelper->intermediateIgnoredOutputBuffer = NULL;
            writerHelper->intermediateIgnoredOutputBufferSize = 0;
        }
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
     * This should be called once for each target, gobyTargetIndex should start at 0
     * but the alignerTargetIndex can be any value. If you use this version to register
     * targets, when you call gobyAlEntry_setTargetIndex you can provide the
     * alignerTargetIndex and it will be converted to the gobyTargetIndex.
     */
    void gobyAlignments_addTargetWithTranslation(CAlignmentsWriterHelper *writerHelper, const unsigned int gobyTargetIndex, const unsigned int alignerTargetIndex, const char *targetName, const unsigned int targetLength) {
        debug(fprintf(stderr,"gobyAlignments_addTargetIdentifier '%s' = [gi=%u, ai=%u] length=%u\n", targetName, gobyTargetIndex, alignerTargetIndex, targetLength));
        string targetNameStr(targetName);
        writerHelper->alignmentWriter->addTargetIdentifier(targetNameStr, gobyTargetIndex);
        writerHelper->alignmentWriter->addTargetLength(targetLength);
        if (writerHelper->alignerToGobyTargetIndexMap == NULL) {
            writerHelper->alignerToGobyTargetIndexMap = new map<unsigned int, unsigned int>;
        }
        (*(writerHelper->alignerToGobyTargetIndexMap))[alignerTargetIndex] = gobyTargetIndex;
    }

    /**
     * This should be called once for each target, targetIndex should start at 0
     * on the first call and increment by one on each following call.
     */
    void gobyAlignments_addTarget(CAlignmentsWriterHelper *writerHelper, const unsigned int targetIndex, const char *targetName, const unsigned int targetLength) {
        debug(fprintf(stderr,"gobyAlignments_addTargetIdentifier '%s' = [ai=%u] length=%u\n", targetName, targetIndex, targetLength));
        string targetNameStr(targetName);
        writerHelper->alignmentWriter->addTargetIdentifier(targetNameStr, targetIndex);
        writerHelper->alignmentWriter->addTargetLength(targetLength);
    }

    /**
     * If you wish to provide the identifier and have the queryIndex generated
     * automatically. If the identifier has already been registered, you'll get
     * back the same queryIndex as before. This should be called if your INPUT
     * is Fasta or Fastq and your query identifier is a real string and you need
     * to gernate a queryIndex automatically.
     */
    unsigned gobyAlignments_addQueryIdentifier(CAlignmentsWriterHelper *writerHelper, const char *queryIdentifier) {
        string queryIdentifierStr(queryIdentifier);
        unsigned queryIndex = writerHelper->alignmentWriter->addQueryIdentifier(queryIdentifierStr);
        debug(fprintf(stderr,"gobyAlignments_addQueryIdentifier %s=%u\n", queryIdentifier, queryIndex));
        return queryIndex;
    }

    /**
     * If you wish to provide the identifier AND the queryIndex. If the 
     * identifier has already been registered, nothing new will happen.
     * This should be called if your input is compact-reads and you HAVE
     * query-identifier strings for your queries. If you don't have query
     * identifiers, don't bother to call this, you know the queryIndex and
     * that's good enough.
     */
    void gobyAlignments_addQueryIdentifierWithIndex(CAlignmentsWriterHelper *writerHelper, const char *queryIdentifier, unsigned int newQueryIndex) {
        debug(fprintf(stderr,"gobyAlignments_addQueryIdentifierWithIndex %s=%u\n", queryIdentifier, newQueryIndex));
        string queryIdentifierStr(queryIdentifier);
        writerHelper->alignmentWriter->addQueryIdentifierWithIndex(queryIdentifierStr, newQueryIndex);
    }
    
    // get an empty alignment entry to populate
    void gobyAlignments_appendEntry(CAlignmentsWriterHelper *writerHelper) {
        debug(fprintf(stderr,"---------------------------------------\n"));
        debug(fprintf(stderr,"gobyAlignments_appendEntry\n"));
        writerHelper->alignmentEntry = writerHelper->alignmentWriter->appendEntry();
        writerHelper->sequenceVariation = NULL;
	    writerHelper->lastSeqVarReadIndex = -1;
        writerHelper->lastSeqVarReadChar = '\0';
        writerHelper->lastSeqVarRefChar = '\0';
    }
    void gobyAlEntry_setMultiplicity(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setMultiplicity=%u\n", value));
        writerHelper->alignmentEntry->set_multiplicity(value);
    }
    void gobyAlEntry_setQueryIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryIndex=%u\n", value));
        if (writerHelper->smallestQueryIndex == -1) {
            writerHelper->smallestQueryIndex = value;
            writerHelper->largestQueryIndex = value;
        } else {
            writerHelper->smallestQueryIndex = min(value, writerHelper->smallestQueryIndex);
            writerHelper->largestQueryIndex = max(value, writerHelper->smallestQueryIndex);
        }
        writerHelper->alignmentEntry->set_query_index(value);
    }
    unsigned int gobyAlEntry_getQueryIndex(CAlignmentsWriterHelper *writerHelper) {
        return writerHelper->alignmentEntry->query_index();
    }

    /**
     * The target index. If you registered targets with gobyAlignments_addTargetWithTranslation()
     * value should in the range of the native aligner (can be 1-based or even any-based).
     * If you registered targets with gobyAlignments_addTarget() you should be providing
     * targets that are 0-based.
     */
    void gobyAlEntry_setTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setTargetIndex=%u\n", value));
        unsigned int indexToUse;
        if (writerHelper->alignerToGobyTargetIndexMap == NULL) {
            indexToUse = value;
        } else {
            indexToUse = (*(writerHelper->alignerToGobyTargetIndexMap))[value];
        }
        writerHelper->alignmentEntry->set_target_index(indexToUse);
    }
    void gobyAlEntry_setPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setPosition=%u\n", value));
        writerHelper->alignmentEntry->set_position(value);
    }
    void gobyAlEntry_setMatchingReverseStrand(CAlignmentsWriterHelper *writerHelper, int value /* bool */) {
        debug(fprintf(stderr,"gobyAlEntry_setMatchingReverseStrand=%d\n", value));
        writerHelper->alignmentEntry->set_matching_reverse_strand(value == 0 ? false : true);
    }
    void gobyAlEntry_setQueryPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryPosition=%u\n", value));
        writerHelper->alignmentEntry->set_query_position(value);
    }
    void gobyAlEntry_setScoreInt(CAlignmentsWriterHelper *writerHelper, int value) {
        float fValue = (float) value;
        debug(fprintf(stderr,"gobyAlEntry_setScore=%f\n", fValue));
        writerHelper->alignmentEntry->set_score(fValue);
    }
    void gobyAlEntry_setNumberOfMismatches(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setNumberOfMismatches=%u\n", value));
        writerHelper->alignmentEntry->set_number_of_mismatches(value);
    }
    void gobyAlEntry_setNumberOfIndels(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setNumberOfIndels=%u\n", value));
        writerHelper->alignmentEntry->set_number_of_indels(value);
    }
    void gobyAlEntry_setQueryAlignedLength(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryAlignedLength=%u\n", value));
        writerHelper->alignmentEntry->set_query_aligned_length(value);
    }
    void gobyAlEntry_setTargetAlignedLength(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setTargetAlignedLength=%u\n", value));
        writerHelper->alignmentEntry->set_target_aligned_length(value);
    }
    void gobyAlEntry_setQueryLength(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        debug(fprintf(stderr,"gobyAlEntry_setQueryLength=%u\n", value));
        writerHelper->alignmentEntry->set_query_length(value);
    }
    void gobyAlEntry_setMappingQuality(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->set_mapping_quality(value);
    }

    void startNewSequenceVariation(CAlignmentsWriterHelper *writerHelper, unsigned int readIndex, unsigned int refPosition) {
        debug(fprintf(stderr,"... startNewSequenceVariation\n"));
        writerHelper->sequenceVariation = writerHelper->alignmentEntry->add_sequence_variations();
        writerHelper->sequenceVariation->set_read_index(readIndex);
        writerHelper->sequenceVariation->set_position(refPosition);
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
    void gobyAlEntry_addSequenceVariation(CAlignmentsWriterHelper *writerHelper, unsigned int readIndex, unsigned int refPosition, char refChar, char readChar, int hasQualCharInt /* bool */, char readQualChar) {
        bool hasQualChar = hasQualCharInt == 0 ? false : true;
        debug(fprintf(stderr,"gobyAlEntry_addSequenceVariation readIndex=%u refPosition=%u ref=%c read=%c hasQualChar=%s\n", readIndex, refPosition, refChar, readChar, hasQualChar ? "true" : "false"));
        if (writerHelper->sequenceVariation == NULL || writerHelper->lastSeqVarReadIndex == -1) {
            // New sequence variation
            startNewSequenceVariation(writerHelper, readIndex, refPosition);
        } else if (writerHelper->lastSeqVarReadIndex + 1 == readIndex || writerHelper->lastSeqVarReadIndex == readIndex) {
            // Append to prev SequenceVar entry
            debug(fprintf(stderr,"... appending to previous seqVar\n"));
        } else if  ((readChar == '-' && writerHelper->lastSeqVarReadChar != '-') ||
                    (refChar == '-' && writerHelper->lastSeqVarRefChar != '-') ||
                    (writerHelper->lastSeqVarReadChar == '-' && readChar != '-')  ||
                    (writerHelper->lastSeqVarRefChar == '-' && refChar != '-')) {
            // transitioning into or out of an insert or delete, make a new SeqVar
            startNewSequenceVariation(writerHelper, readIndex, refPosition);
        } else {
            // Not contiguous to previous SeqVar
            startNewSequenceVariation(writerHelper, readIndex, refPosition);
        }
        writerHelper->lastSeqVarReadIndex = readIndex;
        writerHelper->lastSeqVarReadChar = readChar;
        writerHelper->lastSeqVarRefChar = refChar;
        string *from = writerHelper->sequenceVariation->mutable_from();
        string *to = writerHelper->sequenceVariation->mutable_to();
        (*from) += refChar;
        (*to) += readChar;
        debug(fprintf(stderr,"... sv->read_index=%u sv->position=%u from=%s to=%s\n",
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
        debug(fprintf(stderr,"gobyAlEntry_appendTooManyHits queryIndex=%u alignedLength=%u numberOfHits=%u\n", queryIndex, alignedLength, numberOfHits));
        writerHelper->tmhWriter->append(queryIndex, numberOfHits, alignedLength);
    }

    void gobyAlEntry_setFragmentIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->set_fragment_index(value);
    }

    void gobyAlEntry_setInsertSize(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->set_insert_size(value);
    }

    // These are only used when dealing with a Query Pair
    void gobyAlEntry_setPairFlags(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->set_pair_flags(value);
    }
    void gobyAlEntry_setPairTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->mutable_pair_alignment_link()->set_target_index(value);
    }
    void gobyAlEntry_setPairFragmentIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->mutable_pair_alignment_link()->set_fragment_index(value);
    }
    void gobyAlEntry_setPairPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->mutable_pair_alignment_link()->set_position(value);
    }

    // These are only used when dealing with a Splice
    void gobyAlEntry_setSplicedFlags(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->set_spliced_flags(value);
    }
    void gobyAlEntry_setSplicedTargetIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->mutable_spliced_alignment_link()->set_target_index(value);
    }
    void gobyAlEntry_setSplicedFragmentIndex(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->mutable_spliced_alignment_link()->set_fragment_index(value);
    }
    void gobyAlEntry_setSplicedPosition(CAlignmentsWriterHelper *writerHelper, unsigned int value) {
        writerHelper->alignmentEntry->mutable_spliced_alignment_link()->set_position(value);
    }

    /**
     * Split a string into a vector.
     */
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;
        while(std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }

    /**
     * Input is SAM, ready to be converted to Compact-Alignment. The queryName [col=0] and
     * targetName [col=3] should already be zero based integers.
     * the samInput may contain multiple lines separated by '\n', if multiple lines of input
     * exist, the queryName (first column) should be the same for each of the lines.
     *
     * This is not fully tested or debugged.
     */
    void gobyAlignments_processSAM(CAlignmentsWriterHelper *writerHelper, char *samInput, int npaths, int maxpaths) {
        string samInputStr(samInput);

        std::vector<std::string> samLines;
        std::vector<std::string> samCols;
        split(samInputStr, '\n', samLines);
        
        int sameLineNo;
        for (sameLineNo = 0; sameLineNo < samLines.size(); sameLineNo++) {
            samCols.clear();
            split(samLines[sameLineNo], '\t', samCols);
            
            unsigned int queryIndex = strtoul(samCols[0].c_str(), NULL, 10);
            unsigned int pairFlags = strtoul(samCols[1].c_str(), NULL, 10);
            unsigned int fragmentIndex = 0;
            unsigned int mateFragmentIndex = 0;
            if (pairFlags & SAM_FLAGS_FIRST_READ_P) {
                fragmentIndex = 0;
                mateFragmentIndex = 1;
            } else if (pairFlags & SAM_FLAGS_SECOND_READ_P) {
                fragmentIndex = 1;
                mateFragmentIndex = 0;
            }
            bool matchingReverseStrand = ((pairFlags & SAM_FLAGS_QUERY_MINUSP) > 0);
            
            unsigned int targetIndex = strtoul(samCols[2].c_str(), NULL, 10);
            unsigned int position = strtoul(samCols[3].c_str(), NULL, 10);
            int score = strtol(samCols[4].c_str(), NULL, 10);
            string cigar = samCols[5];
            unsigned int mateTargetIndex;
            bool hasMate;
            if (samCols[6].length() == 0) {
                // No data in mateTargetIndex column
                mateTargetIndex = 0;
                hasMate = false;
            } else if (samCols[6].compare(0, 1, "=") == 0) {
                // Mate on same targetIndex
                mateTargetIndex = targetIndex;
                hasMate = true;
            } else if (samCols[6].compare(0, 1, "*") == 0) {
                // No mate, this value shouldn't be used.
                mateTargetIndex = 0;
                hasMate = false;
            } else {
                // Mate on different targetIndex
                mateTargetIndex = strtoul(samCols[6].c_str(), NULL, 10);
                hasMate = true;
            }
            unsigned int matePosition = strtoul(samCols[7].c_str(), NULL, 10);
            int templateLength = strtol(samCols[8].c_str(), NULL, 10);
            string query = samCols[9];
            string quality;
            if (samCols[10].length() > 0 && samCols[10].compare(0, 1, "*") != 0) {
                quality = samCols[10];
            }
            string md;
            int nm = -1;
            int sm = -1;
            for (int i = 11; i < samCols.size(); i++) {
                if (samCols[i].compare(0, 5, "MD:Z:") == 0) {
                    md = samCols[i].substr(5);
                }
                if (samCols[i].compare(0, 5, "NM:i:") == 0) {
                    nm = strtol(samCols[i].substr(5).c_str(), NULL, 10);
                }
                if (samCols[i].compare(0, 5, "SM:i:") == 0) {
                    sm = strtol(samCols[i].substr(5).c_str(), NULL, 10);
                }
            }
            debug(fprintf(stderr,"%u:%u:%u:%u:%d:%s:%u:%u:%d:%s:%s:%s:%d:%d   ",
                queryIndex,
                pairFlags,
                targetIndex,
                position,
                score,
                cigar.c_str(),
                mateTargetIndex,
                matePosition,
                templateLength,
                query.c_str(),
                quality.c_str(),
                md.c_str(),
                nm,
                sm);)
            debug(fprintf(stderr,"hasMate=%s, fragmentIndex=%u, mateFragmentIndex=%u, matchingReverseStrand=%s\n",
                hasMate?"yes":"no", fragmentIndex, mateFragmentIndex,
                matchingReverseStrand?"yes":"no");)
            CSamHelper *samHelper = samHelper_getResetSamHelper(writerHelper);
            samHelper_setCigar(samHelper, cigar.c_str());
            samHelper_setMd(samHelper, md.c_str());
            samHelper_setQuery(samHelper, query.c_str(), quality.c_str(), query.length(), matchingReverseStrand == true ? 1 : 0);
            debug(fprintf(stderr,"cpp_cigar=%s, cpp_md=%s, cpp_query=%s, cpp_qual=%s\n",
                samHelper->cpp_cigar->c_str(),
                samHelper->cpp_md->c_str(),
                samHelper->cpp_sourceQuery->c_str(),
                samHelper->cpp_sourceQual->c_str());)
            samHelper_constructRefAndQuery(samHelper);
            
            debug(fprintf(stderr,"  score=%d,   qual=%s\n  query=%s\n    ref=%s\n",
                samHelper->score,
                samHelper->cpp_qual->c_str(),
                samHelper->cpp_query->c_str(),
                samHelper->cpp_ref->c_str());)
            debug(fprintf(stderr, "----------------------------------------\n");)
            
        }
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
        samHelper->numIndels = 0;
        samHelper->numMisMatches = 0;
        samHelper->score = 0;
        samHelper->numLeftClipped = 0;
        return samHelper;
    }

    void samHelper_setCigar(CSamHelper *samHelper, const char *cigar) {
        if (cigar) {
            (*samHelper->cpp_cigar) += cigar;
        }
    }

    void samHelper_addCigarItem(CSamHelper *samHelper, int length, char op) {
        std:stringstream lengthStream;
        lengthStream << length;
        (*samHelper->cpp_cigar) += lengthStream.str();
        (*samHelper->cpp_cigar) += op;
    }

    void samHelper_setMd(CSamHelper *samHelper, const char *md) {
        if (md) {
            (*samHelper->cpp_md) += md;
        }
    }

    void samHelper_setQueryTranslate(CSamHelper *samHelper, char *reads, char *qual, 
            unsigned int length, unsigned int reverseStrand) {
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

    void samHelper_setQuery(CSamHelper *samHelper, const char *reads, const char *qual,
            unsigned int length, unsigned int reverseStrand) {
        // TODO: Do we need to reverse this if reverse query?
        if (reads && strlen(reads) > 0) {
            (*samHelper->cpp_sourceQuery) += reads;
        }
        if (qual && strlen(qual) > 0) {
            (*samHelper->cpp_sourceQual) += qual;
        }
    }

    const char *samHelper_getCigarStr(CSamHelper *samHelper) {
        return samHelper->cpp_cigar->c_str();
    }

    void applyCigar(CSamHelper *samHelper) {
        pcrecpp::RE re("([0-9]+)([SMID])");
        pcrecpp::StringPiece input(samHelper->cpp_cigar->c_str());
        debug (fprintf(stderr, ":: Applying cigar=%s\n", samHelper->cpp_cigar->c_str());)
        int length;
        char op;
        int posInReads = 0;
        int i;
        bool startOfCigar = true;
        samHelper->numIndels = 0;
        samHelper->numMisMatches = 0;
        while (re.Consume(&input, &length, &op)) {
            switch(op) {
                case 'S':
                    // Soft clipping
                    for (i = 0; i < length; i++) {
                        (*samHelper->cpp_ref) += '-';
                        (*samHelper->cpp_query) += '-';
                        (*samHelper->cpp_qual) += ((char) samHelper->minQualValue); // min quality
                    }
                    if (startOfCigar) {
                        samHelper->numLeftClipped = length;
                    }
                    posInReads += length;
                    break;
                case 'M':
                    // Account for matches AND mismatches. Any mis-matches will be fixed in applyMd()
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
                    samHelper->numIndels += length;
                    posInReads += length;
                    break;
                case 'D':
                    for (i = 0; i < length; i++) {
                        (*samHelper->cpp_query) += '-';
                        if (samHelper->cpp_sourceQual->size() != 0) {
                            // Minimum qual, placing min value here but it shouldn't get written to
                            // sequence variations
                            (*samHelper->cpp_qual) += ((char) samHelper->minQualValue); // min quality
                        }
                        (*samHelper->cpp_ref) += '?';  // provided in applyMd()
                    }
                    samHelper->numIndels += length;
                    break;
            }
            startOfCigar = false;
        }
    }

    void applyMd(CSamHelper *samHelper) {
        // My RE is simplified from the SAM spec but performs the same task
        // the major difference being mine would allow 5ACG where the current
        // spec would require 5A0C0G0 (which mine will still work with just fine).
        pcrecpp::RE re("([0-9]+|[ACGTN]|\\^[ACGTN]+)");
        pcrecpp::StringPiece input(samHelper->cpp_md->c_str());
        debug (fprintf(stderr, ":: Applying md=%s\n", samHelper->cpp_md->c_str());)
        string mdPart;
        int position = samHelper->numLeftClipped;
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
                    samHelper->numMisMatches++;
                }
            }
        }
    }

    void samHelper_constructRefAndQuery(CSamHelper *samHelper) {
        samHelper->cpp_query->clear();
        samHelper->cpp_qual->clear();
        samHelper->cpp_ref->clear();
        samHelper->alignedLength = 0;
        samHelper->numIndels = 0;
        samHelper->numMisMatches = 0;
        samHelper->score = 0;

        debug (
            fprintf(stderr, ":: Reference and query before construction\n");
            fprintf(stderr, ":: read=%s\n", samHelper->cpp_sourceQuery->c_str());
        )

        applyCigar(samHelper);
        applyMd(samHelper);
        samHelper->alignedLength = samHelper->cpp_query->size();
        samHelper->score = samHelper->alignedLength - samHelper->numIndels - samHelper->numMisMatches;
        
        debug (
            fprintf(stderr, ":: Reference and query constructed via SAM\n");
            fprintf(stderr, ":: ref =%s\n", samHelper->cpp_ref->c_str());
            fprintf(stderr, ":: read=%s\n", samHelper->cpp_query->c_str());
        )
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
            if (writerHelper->alignerToGobyTargetIndexMap != NULL) {
                delete writerHelper->alignerToGobyTargetIndexMap;
            }
            delete writerHelper;
        }
	}

    char* hitTypes[] = {
              "EXACT", "SUB", "INSERTION", "DELETION", "HALFSPLICE_DONOR",
              "HALFSPLICE_ACCEPTOR", "SPLICE", "ONE_THIRD_SHORTEXON",
              "TWO_THIRDS_SHORTEXON", "SHORTEXON", "TERMINAL" };

    void gobyAlignments_debugSequences(CAlignmentsWriterHelper *writerHelper, int hitType, char *refSequence, char *readSequence, int startPos, int hasBeenReversed) {
        debug(
            int prefixLength;
            int suffixLength;
            if (!hasBeenReversed) {
                prefixLength = startPos;
                suffixLength= writerHelper->alignmentEntry->query_length() - strlen(refSequence) - startPos;
            } else {
                prefixLength = writerHelper->alignmentEntry->query_length() - strlen(refSequence) - startPos;
                suffixLength= startPos;
            }

            fprintf(stderr,":: type=%s\n", hitTypes[hitType]);
            fprintf(stderr,":: ref =");
            for (int i = 0; i < prefixLength; i++) {
                fprintf(stderr, "_");
            }
            fprintf(stderr,"%s", refSequence);
            for (int i = 0; i < suffixLength; i++) {
                fprintf(stderr, "_");
            }
            fprintf(stderr,"\n");

            fprintf(stderr,":: read=");
            for (int i = 0; i < prefixLength; i++) {
                fprintf(stderr, "_");
            }
            fprintf(stderr,"%s", readSequence);
            for (int i = 0; i < suffixLength; i++) {
                fprintf(stderr, "_");
            }
            fprintf(stderr,"\n");
        )
    }
}

