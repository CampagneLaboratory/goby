/*
 * CompactFormatsHelpers.h
 *
 *  Created on: Sep 14, 2010
 *      Author: kdorff
 */

#ifndef COMPACTFORMATSHELPERS_H_
#define COMPACTFORMATSHELPERS_H_

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
		unsigned int numRead;
	};
#else
	// Opaque structure for C
	typedef struct {
		void *readsReader;
		void *it;
		void *end;
		void *unopenedFiles;
		unsigned char circular;
		unsigned int numRead;
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
		unsigned int numWritten;
	};
#else
	// Opaque structure for C
	typedef struct {
	    void *alignmentWriter;
	    void *tmhWriter;
	    void *alignmentEntry;
	    void *sequenceVariation;
	    int lastSeqVarReadIndex;
		unsigned int numWritten;
	} CAlignmentsWriterHelper;
#endif

#endif /* COMPACTFORMATSHELPERS_H_ */
