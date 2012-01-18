#include <string>
#include <iostream>
#include <vector>

#include "C_Gsnap.h"
#include "C_Gsnap_structs.h"
#include "C_Alignments.h"

using namespace std;

#undef C_GSNAP_DEBUG
#ifdef C_GSNAP_DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/**
 * C API and functions to support GSNAP alignment / parsing.
 */
extern "C" {

    vector<char *> *gobyGsnap_split(char *str, const char *delim) {
        char* token = strtok(str, delim);
        vector<char *> *result = new vector<char *>();
        while(token != NULL) {
            result->push_back(token);
            token = strtok(NULL, delim);
        }
        return result;
    }

    void gobyGsnap_test_registerTargets(CAlignmentsWriterHelper *writerHelper, char *targetLines) {
        vector<char *> *targets = gobyGsnap_split(targetLines, "\n\r");
        int numTargets = targets->size();
        for (int i = 0; i < numTargets; i++) {
            gobyAlignments_addTarget(writerHelper, i, (*targets)[i], 2000000000);
        }
        delete targets;
    }

    void clearAlignmentSegments(vector<GsnapAlignmentSegment*> *alignmentSegments)  {
        int size = alignmentSegments->size();
        for (int i = 0; i < size; i++) {
            delete (*alignmentSegments)[i];
        }
    }

    void clearAlignmentEntries(vector<GsnapAlignmentEntry*> *alignmentEntries)  {
        int size = alignmentEntries->size();
        for (int i = 0; i < size; i++) {
            vector<GsnapAlignmentSegment*> *alignmentSegments = (*alignmentEntries)[i]->alignmentSegments;
            if (alignmentSegments != NULL) {
                // Delete the elements of the alignmentSegments vector
                clearAlignmentSegments(alignmentSegments);
                // Delete the the alignmentSegments vector
                delete alignmentSegments;
            }
        }
        // Keep the alignmentEntries vector, just clear it.
        alignmentEntries->clear();
    }

    void initializeAlignmentSegment(GsnapAlignmentSegment *alignmentSegment) {
        alignmentSegment->readSequence = NULL;
        alignmentSegment->readStart = 0;
        alignmentSegment->readEnd = 0;
        alignmentSegment->reverseStrand = 0;
        alignmentSegment->targetIdentifier = NULL;
        alignmentSegment->targetIndex = 0;
        alignmentSegment->startType = NULL;
        alignmentSegment->startClip = 0;
        alignmentSegment->endType = NULL;
        alignmentSegment->endClip = 0;
        alignmentSegment->matches = 0;
        alignmentSegment->subs = 0;
        alignmentSegment->alignScore = 0;
        alignmentSegment->mapq = 0;
        alignmentSegment->pairScore = 0;
        alignmentSegment->insertLength = 0;
    }

    void resetAlignment(GsnapAlignment *alignment) {
        alignment->pairedEnd = 0;
        alignment->pairType = NULL;
        alignment->referenceSequence = NULL;
        alignment->referenceQuality = NULL;
        alignment->queryIdentifer = NULL;
        alignment->queryIndex = 0;
        clearAlignmentEntries(alignment->alignmentEntries);
        clearAlignmentEntries(alignment->alignmentEntriesPair);
    }

    void initializeAlignment(GsnapAlignment *alignment) {
        alignment->alignmentEntries = new vector<GsnapAlignmentEntry*>;
        alignment->alignmentEntriesPair = new vector<GsnapAlignmentEntry*>;
        resetAlignment(alignment);
    }

    void gobyGsnap_destoryAlignment(CAlignmentsWriterHelper *writerHelper) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        if (alignment != NULL) {
            resetAlignment(alignment);
            delete alignment->alignmentEntries;
            delete alignment->alignmentEntriesPair;
            delete alignment;
            writerHelper->gsnapAlignment = NULL;
        }
    }

    /**
     * Parse a single gsnap alignment entry
     * @param writerHelper the alignment writer helper
     * @param alignment the gsnap alignment entry to parse. Supports single or paired end.
     * This should be ONE alignment entry (can be two stanzas if paired end)
     */
    void gobyGsnap_parse(CAlignmentsWriterHelper *writerHelper, char *alignmentStr) {
        vector<char *> *lines = gobyGsnap_split(alignmentStr, "\n\r");
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        if (alignment == NULL) {
            alignment = new GsnapAlignment;
            writerHelper->gsnapAlignment = alignment;
            initializeAlignment(alignment);
        } else {
            resetAlignment(alignment);
        }

        int lineNo = 0;

        delete lines;
    }
}
