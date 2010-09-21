#include <string>
#include <iostream>
#include <queue>

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

		printf("Adding %d files to the queue\n", numUnopenedFiles);
		if (numUnopenedFiles == 0) {
			fprintf(stderr,"No input files to process.\n");
			exit(9);
		}
		printf("First item should be %s\n", unopenedFiles[0]);
		printf("Making CReadsHelper\n");
		CReadsHelper *readsHelper = new CReadsHelper;
		printf("... done\n");
		readsHelper->numRead = 0;
		readsHelper->circular = circular;
		printf("Making queue\n");
		readsHelper->unopenedFiles = new queue<string>;
		printf("... done\n");
		for (int i = 0; i < numUnopenedFiles; i++) {
			printf("Grabbing queue item # %d\n", i);
			printf("Grabbing for queue the filename %s\n", unopenedFiles[0]);
			string unopenedFile(unopenedFiles[0]);
			cout << "... adding to the queue: " << unopenedFile << endl;
			readsHelper->unopenedFiles->push(unopenedFile);
			unopenedFiles++;
		}
		string filename = readsHelper->unopenedFiles->front();
		readsHelper->unopenedFiles->pop();
		cout << "Opening file " << filename << endl;
		readsHelper->readsReader = new goby::ReadsReader(filename);
		printf("Making begin\n");
		readsHelper->it = readsHelper->readsReader->beginPointer();
		printf("... done\n");
		return readsHelper;
	}

	/**
	 * This should be called ONCE per read.
	 * Return values, 0 == false,  1 == true
	 *
	 * TODO: * We have a queue of unopenedFile but it isn't currently
	 * TODO:   switching to the next file when we've finished reading
	 * TODO:   the first.
	 */
	int gobyReads_hasNext(CReadsHelper *readsHelper) {
		printf("Has next? ");
		goby::ReadEntryIterator end = (*readsHelper).readsReader->end();
        if (*((*readsHelper).it) != end) {
			printf(" Yes 33\n");
			return 33;
		} else {
			printf(" No\n");
			return 0;
		}
	}

	/**
	 * This should be called ONCE per read AFTER hasNext(...) has been called and returned TRUE.
	 *
	 * TODO: * We aren't supporting paired or circular here. That data should be in the compact-reads
	 * TODO:   file but I'm not reading it at this moment.
	 */
	Sequence_T *gobyReads_next(CReadsHelper *readsHelper) {
		// Not supporting paired reads yet
        printf("Starting next\n");
		printf("Next... last # read was %d\n", (*readsHelper).numRead);
		printf("Getting entry\n");
	    goby::ReadEntry entry = *(*(*readsHelper).it);
	    printf("Got entry\n");
		(*readsHelper).numRead++;
		printf("Incremented newval=%d\n", (*readsHelper).numRead);

	    Sequence_T* queryseq1 = new Sequence_T;
	    // queryseq1 = (Sequence_T) malloc(sizeof(Sequence_T));
	    int fullLength = 0;
	    if (entry.has_sequence()) {
	    	printf("Populating contents\n");
	    	fullLength = entry.sequence().size();
	    	queryseq1->contents_alloc = new char[fullLength + 1];
		    strcpy(queryseq1->contents_alloc, entry.sequence().c_str());
	    	queryseq1->contents_uc_alloc = new char[fullLength + 1];
		    strcpy(queryseq1->contents_uc_alloc, entry.sequence().c_str());
		    printf("Done\n");
	    } else {
	    	printf("No contents\n");
	    	queryseq1->contents_alloc = (char *) NULL;
	    	queryseq1->contents_uc_alloc = (char *) NULL;
	    }
	    queryseq1->contents = queryseq1->contents_alloc;
	    queryseq1->contents_uc = queryseq1->contents_uc_alloc;
	    queryseq1->fulllength = fullLength;
	    queryseq1->chop = (char *) NULL;
	    queryseq1->choplength = 0;
	    queryseq1->skiplength = 0;
	    queryseq1->trimstart = 0;
	    queryseq1->trimend = 0;
	    queryseq1->subseq_offset = 0;

	    printf("populating acc from read index\n");
	    queryseq1->acc = new char[50];
	    sprintf(queryseq1->acc, "%d", entry.read_index());
    	printf("... done\n");
	    if (entry.has_description()) {
	    	printf("Populating description\n");
	    	queryseq1->restofheader = new char[entry.description().size() + 1];
	    	strcpy(queryseq1->restofheader, entry.description().c_str());
	    } else {
	    	printf("No description\n");
	    	queryseq1->restofheader = (char *) NULL;
	    }

	    /**
	     * TODO: Goby is storing Quality Scores modified. I need a decode function.
	     */
	    if (false && entry.has_quality_scores()) {
	    	printf("Populating quality\n");
	    	queryseq1->quality_alloc = new char[entry.quality_scores().size() + 1];
		    strcpy(queryseq1->quality_alloc, entry.quality_scores().c_str());
	    } else {
	    	printf("No quality\n");
	    	queryseq1->quality_alloc = (char *) NULL;
	    }
    	queryseq1->quality = queryseq1->quality_alloc;

	    // Increment to the next ReadsEntry
		printf("Incrementing to next entry\n");
		(*(*readsHelper).it)++;
		printf("Done with next\n");
	    return queryseq1;
	}

	/**
	 * Call after you are _completely_ done reading Goby Reads.
	 */
	void gobyReads_finished(CReadsHelper *readsHelper) {
		  google::protobuf::ShutdownProtobufLibrary();
	}
}
