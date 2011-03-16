/*
 * CompactFormatsHelpers.h
 *
 *  Created on: Sep 14, 2010
 *      Author: kdorff
 */

#ifndef C_COMPACTHELPERS_H_
#define C_COMPACTHELPERS_H_

/**************************************************************
 * Helper for reading Compact-Reads format.
 **************************************************************/
#ifdef __cplusplus
	#include <queue>
	#include <string>
	#include "Reads.h"
	// More complex structure for C++
	struct CReadsHelper {
		goby::ReadsReader *readsReader;
		goby::ReadEntryIterator *it;
		const goby::ReadEntryIterator *end;
		std::queue<std::string> *unopenedFiles;
		unsigned char circular;
		unsigned int numberOfReads;
		// We will reuse the same memory while reading CR data.
		// If the calling app needs to malloc() and free() every time, they can.
	    char *lastReadIdentifier; int lastReadIdentifier_m;
	    char *lastDescription; int lastDescription_m;
	    char *lastSequence; int lastSequence_m;
	    char *lastQuality; int lastQuality_m;
	    char *lastPairSequence; int lastPairSequence_m;
	    char *lastPairQuality; int lastPairQuality_m;
	    int qualityAdjustment;
	};
#else
	// Opaque structure for C
	typedef struct {
		void *readsReader;
		void *it;
		void *end;
		void *unopenedFiles;
		unsigned char circular;
		unsigned int numberOfReads;
		// We will reuse the same memory while reading CR data.
		// If the calling app needs to malloc() and free() every time, they can.
	    char *lastReadIdentifier; int lastReadIdentifier_m;
	    char *lastDescription; int lastDescription_m;
	    char *lastSequence; int lastSequence_m;
	    char *lastQuality; int lastQuality_m;
	    char *lastPairSequence; int lastPairSequence_m;
	    char *lastPairQuality; int lastPairQuality_m;
	    int qualityAdjustment;
	} CReadsHelper;
#endif

/**************************************************************
 * Helper for writing Compact-Alignments format.
 **************************************************************/

#ifdef __cplusplus
    #include "Alignments.h"
    #include "TooManyHits.h"
    // More complex structure for C++
    struct CSamHelper {
        int minQualValue;
        std::string *cpp_cigar;
        std::string *cpp_md;
        std::string *cpp_sourceQuery;
        std::string *cpp_sourceQual;
        std::string *cpp_query;
        std::string *cpp_qual;
        std::string *cpp_ref;
        int alignedLength;
        int numIndels;
        int numMisMatches;
        int score;
    };
    struct CAlignmentsWriterHelper {
        goby::AlignmentWriter *alignmentWriter;
        goby::TooManyHitsWriter *tmhWriter;
        goby::AlignmentEntry *alignmentEntry;
        goby::SequenceVariation *sequenceVariation;
        int lastSeqVarReadIndex;
        unsigned int smallestQueryIndex;
        unsigned int largestQueryIndex;
        unsigned int numberOfAlignedReads;
        CSamHelper *samHelper;
        int qualityAdjustment;
	    FILE *intermediateOutputFile;
	    char *intermediateOutputBuffer;
	    size_t  intermediateOutputBufferSize;
	    FILE *intermediateIgnoredOutputFile;
	    char *intermediateIgnoredOutputBuffer;
	    size_t intermediateIgnoredOutputBufferSize;
    };
#else
	// Opaque structure for C
	typedef struct {
        int minQualValue;
	    void *cpp_cigar;
	    void *cpp_md;
	    void *cpp_sourceQuery;
	    void *cpp_sourceQual;
	    void *cpp_query;
	    void *cpp_qual;
	    void *cpp_ref;
	    int alignedLength;
	    int numIndels;
	    int numMisMatches;
	    int score;
	} CSamHelper;
	typedef struct {
	    void *alignmentWriter;
	    void *tmhWriter;
	    void *alignmentEntry;
	    void *sequenceVariation;
	    int lastSeqVarReadIndex;
	    unsigned int smallestQueryIndex;
	    unsigned int largestQueryIndex;
	    unsigned int numberOfAlignedReads;
	    void *samHelper;
        int qualityAdjustment;
	    FILE *intermediateOutputFile;
	    char *intermediateOutputBuffer;
	    size_t  intermediateOutputBufferSize;
	    FILE *intermediateIgnoredOutputFile;
	    char *intermediateIgnoredOutputBuffer;
	    size_t intermediateIgnoredOutputBufferSize;
	} CAlignmentsWriterHelper;
#endif

#endif /* C_COMPACTHELPERS_H_ */
