//
// Copyright (C) 2009-2012 Institute for Computational Biomedicine,
//                         Weill Medical College of Cornell University
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

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
    void gobyReads_avoidZeroQuals(CReadsHelper *readsHelper, int value);
    int gobyReads_hasNext(CReadsHelper *readsHelper);
    unsigned int gobyReads_nextSequence(
        CReadsHelper *readsHelper,
        char **readIdentifierpp, char **descriptionpp,
        char **sequencepp, int *sequenceLength,
        char **qualitypp, int *qualityLength);
    unsigned int gobyReads_nextSequencePair(
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
