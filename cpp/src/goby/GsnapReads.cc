#include <string>
#include <iostream>
#include <vector>

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
	CReadsHelper* gobyReads_openReadsReader(
			char **unopenedFiles, int numUnopenedFiles,
			unsigned char circular) {
		if (numUnopenedFiles == 0) {
            fprintf(stderr,"No input files to process.\n");
            exit(9);
		}
		CReadsHelper *readsHelper = new CReadsHelper;
		readsHelper->it = readsHelper->readsReader->beginPointer();
		readsHelper->numRead = 0;
		readsHelper->circular = circular;
		readsHelper->unopenedFiles = new queue<string>;
		for (int i = numUnopenedFiles; i < numUnopenedFiles; i--) {
			string unopenedFile(unopenedFiles[0]);
			readsHelper->unopenedFiles->push(unopenedFile);
			unopenedFiles++;
		}
		string filename = readsHelper->unopenedFiles->front();
		readsHelper->unopenedFiles->pop();
		cout << "Opening file " << filename << endl;
		readsHelper->readsReader = new goby::ReadsReader(filename);
		return readsHelper;
	}

	/**
	 * This should be called ONCE per read.
	 * Return values, 0 == false,  1 == true
	 *
	 * TODO: * We have a vector of unopenedFile but it isn't currently
	 * TODO:   switching to the next file when we've finished reading
	 * TODO:   the first.
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
	 *
	 * TODO: * We aren't supporting paired or circular here. That data should be in the compact-reads
	 * TODO:   file but I'm not reading it at this moment.
	 */
	Sequence_T gobyReads_next(CReadsHelper *readsHelper) {
		// Not supporting paired reads yet

		readsHelper->numRead++;
	    goby::ReadEntry entry = *(*(readsHelper->it));

	    Sequence_T queryseq1;
	    int fullLength = 0;
	    if (entry.has_sequence()) {
	    	fullLength = entry.sequence().size();
	    	queryseq1.contents_alloc = new char[fullLength + 1];
		    strcpy(queryseq1.contents_alloc, entry.sequence().c_str());
	    	queryseq1.contents_uc_alloc = new char[fullLength + 1];
		    strcpy(queryseq1.contents_uc_alloc, entry.sequence().c_str());
	    } else {
	    	queryseq1.contents_alloc = (char *) NULL;
	    	queryseq1.contents_uc_alloc = (char *) NULL;
	    }
	    queryseq1.contents = queryseq1.contents_alloc;
	    queryseq1.contents_uc = queryseq1.contents_uc_alloc;
	    queryseq1.fulllength = fullLength;
	    queryseq1.chop = (char *) NULL;
	    queryseq1.choplength = 0;
	    queryseq1.skiplength = 0;
	    queryseq1.trimstart = 0;
	    queryseq1.trimend = 0;
	    queryseq1.subseq_offset = 0;

	    sprintf(queryseq1.acc, "%d", entry.read_index());
	    if (entry.has_description()) {
	    	queryseq1.restofheader = new char[entry.description().size() + 1];
	    	strcpy(queryseq1.restofheader, entry.description().c_str());
	    } else {
	    	queryseq1.restofheader = (char *) NULL;
	    }

	    if (entry.has_quality_scores()) {
	    	queryseq1.quality_alloc = new char[entry.quality_scores().size() + 1];
		    strcpy(queryseq1.quality_alloc, entry.quality_scores().c_str());
	    } else {
	    	queryseq1.quality_alloc = (char *) NULL;
	    }
	    queryseq1.quality = queryseq1.quality_alloc;

	    // Increment to the next ReadsEntry
		(*(readsHelper->it))++;
	    return queryseq1;
	}

	/**
	 * Call after you are _completely_ done reading Goby Reads.
	 */
	void gobyReads_finished(CReadsHelper *readsHelper) {
		  google::protobuf::ShutdownProtobufLibrary();
	}
}
