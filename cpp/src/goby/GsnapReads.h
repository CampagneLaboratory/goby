/*
 * GsnapReads.h
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */

#ifndef GSNAPREADS_H_
#define GSNAPREADS_H_

#include "CReadsHelper.h"

#ifdef __cplusplus
#include "GsnapSequence.h"
extern "C" {
#endif
	CReadsHelper *openReadsReader(char *filename);
	int hasNext(CReadsHelper *readsHelper);
	// Sequence_T* next(CReadsHelper *readsHelper);
	void next(CReadsHelper *readsHelper);
	void finished();
#ifdef __cplusplus
}
#endif

#endif /* GSNAPREADS_H_ */
