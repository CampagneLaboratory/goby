//
// Copyright (C) 2009-2012 Institute for Computational Biomedicine,
//                         Weill Medical College of Cornell University
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <stdio.h>
#include <goby/C_Reads.h>

int main(int argc, char** argv) {
    if (argc < 2) { 
        printf("Specify Goby compact-reads file to open.\n");
        return;
    }
    char *input = argv[1];
    CReadsHelper *targetReadsHelper;
    gobyReads_openReadsReaderSingleWindowed(input, 0, 0, &targetReadsHelper);
    while (gobyReads_hasNext(targetReadsHelper)) {
        char *readIdentifier;
        char *description;
        char *sequence;
        int sequenceLength;
        char *quality;
        int qualityLength;
        unsigned long readIndex = gobyReads_nextSequence(
            targetReadsHelper,
            &readIdentifier, &description,
            &sequence, &sequenceLength,
            &quality, &qualityLength);
            
        printf("read-index: %d read-id: %s sequence: %s\n",
            readIndex, readIdentifier, sequence);
    }
    gobyReads_finished(targetReadsHelper);
    goby_shutdownProtobuf();
}
