#include <string>
#include <iostream>
#include <queue>

#include "Reads.h"
#include "GsnapStructs.h"
#include "GsnapReads.h"

using namespace std;

#ifdef GSNAP_READ_COMPACT_DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/**
 * This class is a C interface so Gsnap can read Goby compact-reads.
 */
extern "C" {

	void gobyReads_openReadsReader(
			char **unopenedFiles, int numUnopenedFiles, unsigned char circular, CReadsHelper **readsHelperpp) {
	    gobyReads_openReadsReaderWindowed(unopenedFiles, numUnopenedFiles, circular, 0, 0, readsHelperpp);
    }

	void gobyReads_openReadsReaderSingleWindowed(
			char *filename,  unsigned long startOffset, unsigned long endOffset, CReadsHelper **readsHelperpp) {
	    gobyReads_openReadsReaderWindowed(&filename, 1, 0, startOffset, endOffset, readsHelperpp);
    }

	/**
	 * Open the .compact-reads file to read from.
	 */
	void gobyReads_openReadsReaderWindowed(
			char **unopenedFiles, int numUnopenedFiles, unsigned char circular,
			unsigned long startOffset, unsigned long endOffset, CReadsHelper **readsHelperpp) {

		if (numUnopenedFiles == 0) {
			fprintf(stderr,"No input files to process.\n");
			exit(9);
		}
        *readsHelperpp = new CReadsHelper;
		CReadsHelper *readsHelper = *readsHelperpp;
		readsHelper->numberOfReads = 0;
		readsHelper->circular = circular;
		readsHelper->unopenedFiles = new queue<string>;
		for (int i = 0; i < numUnopenedFiles; i++) {
			string unopenedFile(unopenedFiles[0]);
			readsHelper->unopenedFiles->push(unopenedFile);
			unopenedFiles++;
		}
		string filename = readsHelper->unopenedFiles->front();
		readsHelper->unopenedFiles->pop();
		readsHelper->readsReader = new goby::ReadsReader(filename);
		readsHelper->it = readsHelper->readsReader->beginPointer(startOffset, endOffset);
		readsHelper->end = readsHelper->readsReader->endPointer();
		readsHelper->numberOfReads = 0;
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
        if (*((*readsHelper).it) != *(readsHelper->end)) {
			return 1;
		} else {
			return 0;
		}
	}

	/**
     * This should be called ONCE per read AFTER hasNext(...) has been called and returned TRUE.
     * This is specicially for Gsnap. It is important that the _alloc portions of queryseq1, queryseq2
     * are free'd by the calling code.
     * See the more generic gobyReads_nextSequence(), gobyReads_nextSequencePair()
     * for non-GSnap implementations.
	 *
	 * TODO: * We aren't supporting paired or circular here. That data should be in the compact-reads
	 * TODO:   file but I'm not reading it at this moment.
	 */
	void gobyReads_next(CReadsHelper *readsHelper, Sequence_T **queryseq1pp, Sequence_T **queryseq2pp) {
        goby::ReadEntry entry = *(*(*readsHelper).it);
		(*readsHelper).numberOfReads++;

	    *queryseq1pp = (Sequence_T*) malloc(sizeof(Sequence_T));
        Sequence_T *queryseq1 = *queryseq1pp;

	    int fullLength = 0;
	    if (entry.has_sequence()) {
	    	fullLength = entry.sequence().size();
	    	queryseq1->contents_alloc = (char *) calloc(fullLength + 1, sizeof(char));
		    strcpy(queryseq1->contents_alloc, entry.sequence().c_str());
	    	queryseq1->contents_uc_alloc = (char *) calloc(fullLength + 1, sizeof(char));
		    strcpy(queryseq1->contents_uc_alloc, entry.sequence().c_str());
	    } else {
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
        //TODO introduce a read index and preserve read id in acc?
	    queryseq1->acc = (char *) malloc(50);
	    sprintf(queryseq1->acc, "%d", entry.read_index());
	    if (entry.has_description()) {
	    	queryseq1->restofheader = (char *) calloc(entry.description().size() + 1, sizeof(char));
	    	strcpy(queryseq1->restofheader, entry.description().c_str());
	    } else {
	    	queryseq1->restofheader = (char *) NULL;
	    }

	    /**
	     * TODO: Goby is storing Quality Scores in Phred units. What encoding does GSnap require?
	     * We need to convert quality score appropriately here.
	     */
	    if (entry.has_quality_scores()) {
	        int qualSize = entry.quality_scores().size();
	    	queryseq1->quality_alloc = (char *) malloc(qualSize + 1);
		    memcpy(queryseq1->quality_alloc, entry.quality_scores().c_str(), qualSize);
		    queryseq1->quality_alloc[qualSize] = '\0';
	    } else {
	    	queryseq1->quality_alloc = (char *) NULL;
	    }
    	queryseq1->quality = queryseq1->quality_alloc;

        *queryseq2pp = NULL;
	    if (entry.has_sequence_pair()) {
            // Populate the paired sequence into queryseq2
            fullLength = entry.sequence_pair().size();
            if (fullLength > 0) {
                *queryseq2pp = (Sequence_T*) malloc(sizeof(Sequence_T));
                Sequence_T *queryseq2 = *queryseq2pp;

                queryseq2->contents_alloc = (char *) calloc(fullLength + 1, sizeof(char));
                strcpy(queryseq2->contents_alloc, entry.sequence_pair().c_str());
                queryseq2->contents_uc_alloc = (char *) calloc(fullLength + 1, sizeof(char));
                strcpy(queryseq2->contents_uc_alloc, entry.sequence_pair().c_str());

                queryseq2->contents = queryseq2->contents_alloc;
                queryseq2->contents_uc = queryseq2->contents_uc_alloc;
                queryseq2->fulllength = fullLength;
                queryseq2->chop = (char *) NULL;
                queryseq2->choplength = 0;
                queryseq2->skiplength = 0;
                queryseq2->trimstart = 0;
                queryseq2->trimend = 0;
                queryseq2->subseq_offset = 0;
                //TODO introduce a read index and preserve read id in acc?
                queryseq2->acc = (char *) malloc(strlen(queryseq1->acc) + 1);
                strcpy(queryseq2->acc, queryseq1->acc);
                if (queryseq1->restofheader != NULL) {
                    queryseq2->restofheader = (char *) malloc(strlen(queryseq1->restofheader) + 1);
                    strcpy(queryseq2->restofheader, queryseq1->restofheader);
                } else {
                    queryseq2->restofheader = (char *) NULL;
                }

                /**
                 * TODO: Goby is storing Quality Scores in Phred units. What encoding does GSnap require?
                 * We need to convert quality score appropriately here.
                 */
                if (entry.has_quality_scores_pair()) {
                    int qualSize = entry.quality_scores_pair().size();
                    queryseq2->quality_alloc = (char *) malloc(qualSize + 1);
                    memcpy(queryseq2->quality_alloc, entry.quality_scores_pair().c_str(), qualSize);
                    queryseq2->quality_alloc[qualSize] = '\0';
                } else {
                    queryseq2->quality_alloc = (char *) NULL;
                }
                queryseq2->quality = queryseq2->quality_alloc;
            }
	    }

	    // Increment to the next ReadsEntry
		(*(*readsHelper).it)++;
	}

	void transferString(string src, char **dest) {
        int size = src.size();
        *dest = (char *) malloc(size + 1);
        memcpy(*dest, src.c_str(), size);
        (*dest)[size] = '\0';
	}

    /**
     * For implementations OTHER than GSnap. Read the sequence but ignore the pair even if it exists.
     * It is important that each of thes char ** parameters are free()'d in the calling code.
     * @return the Goby read index
     */
	unsigned long gobyReads_nextSequence(
	    CReadsHelper *readsHelper,
	    char **readIdentifierpp, char **descriptionpp,
	    char **sequencepp, int *sequenceLength,
	    char **qualitypp, int *qualityLength) {

        /** Default is nothing populated. */
	    *readIdentifierpp = NULL;
	    *descriptionpp = NULL;
	    *sequencepp = NULL;
	    *sequenceLength = 0;
	    *qualitypp = NULL;
	    *qualityLength = 0;

        goby::ReadEntry entry = *(*(*readsHelper).it);
		(*readsHelper).numberOfReads++;

	    if (entry.has_read_identifier()) {
	        transferString(entry.read_identifier(), readIdentifierpp);
	    }

	    if (entry.has_description()) {
	        transferString(entry.description(), descriptionpp);
	    }

	    if (entry.has_sequence()) {
	        transferString(entry.sequence(), sequencepp);
	        *sequenceLength = entry.sequence().size();
	        if (entry.has_quality_scores()) {
	            transferString(entry.quality_scores(), qualitypp);
	            *qualityLength = entry.quality_scores().size();
            }
	    }

	    // Increment to the next ReadsEntry
		(*(*readsHelper).it)++;

	    return entry.read_index();
    }

    /**
     * For implementations OTHER than GSnap. Read the sequence WITH pair if applicable.
     * It is important that each of thes char ** parameters are free()'d in the calling code.
     * @return the Goby read index
     */
	unsigned long gobyReads_nextSequencePair(
	    CReadsHelper *readsHelper,
	    char **readIdentifierpp, char **descriptionpp,
	    char **sequencepp, int *sequenceLength,
	    char **qualitypp, int *qualityLength,
	    char **pairSequencepp, int *pairSequenceLength,
	    char **pairQualitypp, int *pairQualityLength) {

        /** Default is nothing populated. */
	    *readIdentifierpp = NULL;
	    *descriptionpp = NULL;
	    *sequencepp = NULL;
	    *qualitypp = NULL;
	    *qualityLength = 0;
	    *pairSequencepp = NULL;
	    *pairQualitypp = NULL;
	    *pairQualityLength = 0;

        goby::ReadEntry entry = *(*(*readsHelper).it);
		(*readsHelper).numberOfReads++;

	    if (entry.has_read_identifier()) {
	        transferString(entry.read_identifier(), readIdentifierpp);
	    }

	    if (entry.has_description()) {
	        transferString(entry.description(), descriptionpp);
	    }

	    if (entry.has_sequence()) {
	        transferString(entry.sequence(), sequencepp);
	        *sequenceLength = entry.sequence().size();
	        if (entry.has_quality_scores()) {
	            transferString(entry.quality_scores(), qualitypp);
	            *qualityLength = entry.quality_scores().size();
            }
	    }

	    if (entry.has_sequence_pair()) {
	        transferString(entry.sequence_pair(), pairSequencepp);
	        *pairSequenceLength = entry.sequence_pair().size();
	        if (entry.has_quality_scores_pair()) {
	            transferString(entry.quality_scores_pair(), pairQualitypp);
	            *pairQualityLength = entry.quality_scores_pair().size();
	        }
	    }

	    // Increment to the next ReadsEntry
		(*(*readsHelper).it)++;

	    return entry.read_index();
	}

	/**
	 * Call after you are _completely_ done reading Goby Reads.
	 */
	void gobyReads_finished(CReadsHelper *readsHelper) {
	    if (readsHelper != NULL) {
            while (!readsHelper->unopenedFiles->empty()) {
                string unopenedFile = readsHelper->unopenedFiles->front();
                readsHelper->unopenedFiles->pop();
            }
            delete readsHelper->unopenedFiles;
            delete readsHelper->readsReader;
            delete readsHelper->it;
            delete readsHelper;
        }
	}

	void goby_shutdownProtobuf() {
        google::protobuf::ShutdownProtobufLibrary();
    }
}
