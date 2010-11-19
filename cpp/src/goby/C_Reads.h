/*
 * Definition of functions to enable reading Goby compact-reads in C.
 */

#ifndef C_READS_H_
#define C_READS_H_

#include "C_CompactHelpers.h"

#ifdef __cplusplus
extern "C" {
#endif
	void gobyReads_openReadsReader(
			char **unopenedFiles, int numUnopenedFiles, unsigned char circular, CReadsHelper **readsHelperpp);
	void gobyReads_openReadsReaderSingleWindowed(
			char *filename,  unsigned long startOffset, unsigned long endOffset, CReadsHelper **readsHelperpp);
	void gobyReads_openReadsReaderWindowed(
			char **unopenedFiles, int numUnopenedFiles, unsigned char circular,
			unsigned long startOffset, unsigned long endOffset, CReadsHelper **readsHelperpp);

	int gobyReads_getQualityAdjustment(CReadsHelper *readsHelper);
	void gobyReads_setQualityAdjustment(CReadsHelper *readsHelper, int value);
	int gobyReads_hasNext(CReadsHelper *readsHelper);
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

#endif /* C_READS_H_ */
