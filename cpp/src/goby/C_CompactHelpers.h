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
	struct CAlignmentsWriterHelper {
	    goby::AlignmentWriter *alignmentWriter;
	    goby::TooManyHitsWriter *tmhWriter;
	    goby::AlignmentEntry *alignmentEntry;
	    goby::SequenceVariation *sequenceVariation;
	    int lastSeqVarReadIndex;
	    unsigned int smallestQueryIndex;
	    unsigned int largestQueryIndex;
	    unsigned int numberOfAlignedReads;
	    std::string *currentCigar;
	    std::string *currentMd;
	    std::string *currentSourceQuery;
	    std::string *currentSourceQual;
	    std::string *currentQuery;
	    std::string *currentQual;
	    std::string *currentRef;
	    int currentAlignedLength;
	    int currentStartPosition;
	    int currentNumIndels;
	    int currentMisMatches;
	    int currentScore;
	};
#else
	// Opaque structure for C
	typedef struct {
	    void *alignmentWriter;
	    void *tmhWriter;
	    void *alignmentEntry;
	    void *sequenceVariation;
	    int lastSeqVarReadIndex;
	    unsigned int smallestQueryIndex;
	    unsigned int largestQueryIndex;
	    unsigned int numberOfAlignedReads;
	    void *currentCigar;
	    void *currentMd;
	    void *currentSourceQuery;
	    void *currentSourceQual;
	    void *currentQuery;
	    void *currentQual;
	    void *currentRef;
	    int currentAlignedLength;
	    int currentStartPosition;
	    int currentNumIndels;
	    int currentMisMatches;
	    int currentScore;
	} CAlignmentsWriterHelper;
#endif

#endif /* C_COMPACTHELPERS_H_ */
