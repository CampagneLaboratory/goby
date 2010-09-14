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
	goby::ReadsReader* openReadsReader(char *filename) {
		string strFilename (filename);
		cout << "Opening file " << strFilename << endl;
		return new goby::ReadsReader(strFilename);
	}

	/**
	 * Obtain an iterator to read entries from the already open .compact-reads file.
	 */
	goby::ReadEntryIterator* getReadsIterator(goby::ReadsReader *readerReader) {
		return readerReader->beginPointer();
	}

	/**
	 * This should be called ONCE per read.
	 * Return values, 0 == false,  1 == true
	 */
	int hasNext(goby::ReadsReader *readsReader, goby::ReadEntryIterator *it, int numRead) {
		if ((*it) != readsReader->end()) {
			if (numRead != 0) {
				(*it)++;
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
	void next(goby::ReadEntryIterator *it) {
	    goby::ReadEntry entry = *(*it);
	}

	/**
	 * Call after you are _completely_ done reading Goby Reads.
	 */
	void finished() {
		  google::protobuf::ShutdownProtobufLibrary();
	}

	int testReturnsTen() {
		return 10;
	}
}
