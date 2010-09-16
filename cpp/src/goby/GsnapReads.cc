#include <string>
#include <iostream>
#include "Reads.h"
#include "GsnapSequence.h"
#include "GsnapReads.h"

using namespace std;
/**
 * This class is a C interface so GSnap can reads from a Goby ReadsReader.
 */
extern "C" {

	/**
	 * Open the .compact-reads file to read from.
	 */
	CReadsHelper* gobyReads_openReadsReader(char *filename) {
		string strFilename (filename);
		cout << "Opening file " << strFilename << endl;
		CReadsHelper *readsHelper = new CReadsHelper;
		readsHelper->readsReader = new goby::ReadsReader(strFilename);
		readsHelper->it = readsHelper->readsReader->beginPointer();
		readsHelper->numRead = 0;
		return readsHelper;
	}

	/**
	 * This should be called ONCE per read.
	 * Return values, 0 == false,  1 == true
	 */
	int gobyReads_hasNext(CReadsHelper *readsHelper) {
		if (*(readsHelper->it) != readsHelper->readsReader->end()) {
			return 1;
		} else {
			return 0;
		}
	}

	/**
	 * This should be called ONCE per read AFTER hasNext(...) has been called and returned TRUE.
	 */
	Sequence_T gobyReads_next(CReadsHelper *readsHelper) {
		readsHelper->numRead++;
	    goby::ReadEntry entry = *(*(readsHelper->it));
	    Sequence_T query;

	    if (entry.has_sequence()) {
		    query.contents = new char[entry.sequence().size()+1];
		    strcpy (query.contents, entry.sequence().c_str());
	    } else {
		    query.contents = (char *) NULL;
	    }

	    // Increment to the next ReadsEntry
		(*(readsHelper->it))++;
	    return query;
	}

	/**
	 * Call after you are _completely_ done reading Goby Reads.
	 */
	void gobyReads_finished(CReadsHelper *readsHelper) {
		  google::protobuf::ShutdownProtobufLibrary();
	}
}
