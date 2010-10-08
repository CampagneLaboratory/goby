#include <string>
#include <iostream>

#include "Reads.h"
#include "GsnapSequence.h"
#include "GsnapAlignments.h"
#include "MessageChunks.h"

using namespace std;

/**
 * This class is a C interface so Gsnap can write Goby compact-alignments.
 */
extern "C" {

	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename) {
	    return gobyAlignments_openAlignmentsWriter(basename, GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK);
	}

	CAlignmentsWriterHelper *gobyAlignments_openAlignmentsWriter(char *basename, unsigned number_of_entries_per_chunk) {
		CAlignmentsWriterHelper *writerHelper = new CAlignmentsWriterHelper;
        string basenameToWriteTo(basename);
        writerHelper->alignmentWriter = new goby::AlignmentWriter(filename, number_of_entries_per_chunk);

		writerHelper->numRead = 0;
        return writerHelper;
	}

	void gobyAlignments_finished(CAlignmentsWriterHelper *writerHelper) {
        if (writerHelper != NULL) {
            writerHelper->close();
            delete writerHelper->alignmentWriter;
            delete writerHelper;
        }
	}

}