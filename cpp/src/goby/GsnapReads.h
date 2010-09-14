/*
 * GsnapReads.h
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */

#ifndef GSNAPREADS_H_
#define GSNAPREADS_H_


#ifdef __cplusplus
extern "C" {
	#include "Reads.h"
	#include "GsnapSequence.h"
	goby::ReadsReader* openReadsReader(char *filename);
	goby::ReadEntryIterator* getReadsIterator(goby::ReadsReader *readerReader);
	int hasNext(goby::ReadsReader *readerReaderP, goby::ReadEntryIterator *itP, int numRead);
	// Sequence_T* next(goby::ReadEntryIterator *itP);
	void next(goby::ReadEntryIterator *itP);
	void finished();
	int testReturnsTen();
}
#else
	void* openReadsReader(char *filename);
	void* getReadsIterator(void *reader);
	int hasNext(void *reader, void *iterator, int numRead);
	// Sequence_T* next(void *iterator);
	void next(void *iterator);
	void finished();
	int testReturnsTen();
#endif

#endif /* GSNAPREADS_H_ */
