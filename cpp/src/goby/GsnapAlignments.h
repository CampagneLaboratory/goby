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
	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename);
	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename, unsigned int number_of_entries_per_chunk);

	void gobyAlignments_finished(CAlignmentsWriterHelper *alWriterHelper);
#ifdef __cplusplus
}
#endif

#endif /* GSNAPREADS_H_ */
