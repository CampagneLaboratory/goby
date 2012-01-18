//
// Copyright (C) 2010 Institute for Computational Biomedicine,
//                    Weill Medical College of Cornell University
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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include "goby/C_Gsnap.h"
#include "goby/C_Alignments.h"
#include "goby/C_Reads.h"

using namespace std;

char *read_file(char *filename) {
    FILE *pFile = fopen(filename, "r" );
    if (pFile == NULL) {
        fputs("File error", stderr);
        exit(1);
    }

    // obtain file size:
    fseek(pFile, 0, SEEK_END);
    long lSize = ftell(pFile);
    rewind(pFile);

    // allocate memory to contain the whole file + null terminator for the string
    char *buffer = (char*) malloc((sizeof(char) * lSize) + sizeof(char));
    if (buffer == NULL) {
        fputs("Memory error", stderr);
        exit(1);
    }

    // copy the file into the buffer:
    size_t result = fread(buffer, 1, lSize, pFile);
    if (result != lSize) {
        fputs("Reading error", stderr);
        exit(1);
    }
    // null terminate the buffer
    buffer[lSize] = '\0';

    fclose (pFile);
    return buffer;
}

void registerChromosomes(CAlignmentsWriterHelper *writerHelper) {
    char *test1 = read_file("test-data/zebrafish.chromosome.list.txt");
    gobyGsnap_test_registerTargets(writerHelper, test1);
    free(test1);
}

void targetIdentiferTest(CAlignmentsWriterHelper *writerHelper) {
    cout << "MT?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "MT") << " ";
    cout << "Zv9...?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Zv9_scaffold3564") << " ";
    cout << "25?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "25") << " ";
    cout << "z?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Z") << endl;

    unsigned t_mt = gobyAlignments_targetIndexForIdentifier(writerHelper, "MT");
    unsigned t_zv9 = gobyAlignments_targetIndexForIdentifier(writerHelper, "Zv9_scaffold3564");
    unsigned t_25 = gobyAlignments_targetIndexForIdentifier(writerHelper, "25");
    unsigned t_z = gobyAlignments_targetIndexForIdentifier(writerHelper, "Z");
    cout << "MT=" << t_mt << " " << "Zv9...=" << t_zv9 << " " << "25=" << t_25 << " " << "z=" << t_z << endl;

    // Verify that getting targetIndex didn't effect the value of the index
    cout << "MT?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "MT") << " ";
    cout << "Zv9...?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Zv9_scaffold3564") << " ";
    cout << "25?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "25") << " ";
    cout << "z?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Z") << endl;
}

void testPairedEnd(CAlignmentsWriterHelper *writerHelper) {
    char *test1 = read_file("test-data/gsnap-output-pair-test-1.gsnap");
    gobyGsnap_parse(writerHelper, test1);
    free(test1);
}


int main(int argc, const char *const argv[]) {
    CAlignmentsWriterHelper *writerHelper;
    gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk("deleteme", &writerHelper);
    registerChromosomes(writerHelper);
    targetIdentiferTest(writerHelper);
    testPairedEnd(writerHelper);
    gobyAlignments_finished(writerHelper, 1);
    goby_shutdownProtobuf();
    return 0;
}
