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
	CReadsHelper *gobyReads_openReadsReader(char **unopenedFiles, int numUnopenedFiles, unsigned char circular);
	int gobyReads_hasNext(CReadsHelper *readsHelper);
	Sequence_T *gobyReads_next(CReadsHelper *readsHelper);
	void gobyReads_finished(CReadsHelper *readsHelper);
	void goby_shutdownProtobuf();
#ifdef __cplusplus
}
#endif

#endif /* GSNAPREADS_H_ */
