/*
 * GsnapReads.h
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */

#ifndef GSNAPREADS_H_
#define GSNAPREADS_H_

#include "CompactFormatsHelpers.h"

#ifdef __cplusplus
#include "GsnapStructs.h"
extern "C" {
#endif
	void gobyReads_openReadsReader(
			char **unopenedFiles, int numUnopenedFiles, unsigned char circular, CReadsHelper **readsHelperpp);
	void gobyReads_openReadsReaderSingleWindowed(
			char *filename,  unsigned long startOffset, unsigned long endOffset, CReadsHelper **readsHelperpp);
	void gobyReads_openReadsReaderWindowed(
			char **unopenedFiles, int numUnopenedFiles, unsigned char circular,
			unsigned long startOffset, unsigned long endOffset, CReadsHelper **readsHelperpp);

	int gobyReads_hasNext(CReadsHelper *readsHelper);
	void gobyReads_next(CReadsHelper *readsHelper, Sequence_T **queryseq1pp, Sequence_T **queryseq2pp);
	unsigned long gobyReads_nextSequence(
	    CReadsHelper *readsHelper,
	    char **readIdentifierpp, char **descriptionpp,
	    char **sequencepp, int *sequenceLength,
	    char **qualitypp, int *qualityLength);
	unsigned long gobyReads_nextSequencePair(
	    CReadsHelper *readsHelper,
	    char **readIdentifierpp, char **descriptionpp,
	    char **sequencepp, int *sequenceLength,
	    char **qualitypp, int *qualityLength,
	    char **pairSequencepp, int *pairSequenceLength,
	    char **pairQualitypp, int *pairQualityLength);
	void gobyReads_finished(CReadsHelper *readsHelper);
	void goby_shutdownProtobuf();
#ifdef __cplusplus
}
#endif

#endif /* GSNAPREADS_H_ */
