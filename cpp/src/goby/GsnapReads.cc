#include <string>
#include <iostream>
#include "Reads.h"
#include "GSnapSequence.h"
#include "GSnapReads.h"

using namespace std;
/**
 * This class is a C interface so GSnap can reads from a Goby ReadsReader.
 */
extern "C" {

	/**
	 * Open the .compact-reads file to read from.
	 */
	CReadsHelper* openReadsReader(char *filename) {
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
	int hasNext(CReadsHelper *readsHelper) {
		if (*(readsHelper->it) != readsHelper->readsReader->end()) {
			if (readsHelper->numRead != 0) {
				(*(readsHelper->it))++;
			}
			return 1;
		} else {
			return 0;
		}
	}

	/**
	 * This should be called ONCE per read AFTER hasNext(...) has been called and returned TRUE.
	 */
//	Sequence_T* next(goby::ReadEntryIterator it) {
	void next(CReadsHelper *readsHelper) {
	    goby::ReadEntry entry = *(*(readsHelper->it));
		readsHelper->numRead++;
	}

	/**
	 * Call after you are _completely_ done reading Goby Reads.
	 */
	void finished() {
		  google::protobuf::ShutdownProtobufLibrary();
	}
}
