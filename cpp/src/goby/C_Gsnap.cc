#include <string>
#include <iostream>
#include <vector>
#include <stdlib.h>

#include "C_Gsnap.h"
#include "C_Gsnap_structs.h"
#include "C_Alignments.h"
#include "hash.h"

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
            gobyAlignments_addTarget(writerHelper, i, targets->at(i), 2000000000);
        }
        delete targets;
    }

    void clearAlignmentSegments(vector<GsnapAlignmentSegment*> *alignmentSegments)  {
        int size = alignmentSegments->size();
        for (int i = 0; i < size; i++) {
            delete alignmentSegments->at(i);
        }
    }

    void clearAlignmentEntries(vector<GsnapAlignmentEntry*> *alignmentEntries)  {
        int size = alignmentEntries->size();
        for (int i = 0; i < size; i++) {
            GsnapAlignmentEntry *alignmentEntry = alignmentEntries->at(i);
            vector<GsnapAlignmentSegment*> *alignmentSegments = alignmentEntry->alignmentSegments;
            if (alignmentSegments != NULL) {
                // Delete the elements of the alignmentSegments vector
                clearAlignmentSegments(alignmentSegments);
                // Delete the the alignmentSegments vector
                delete alignmentSegments;
            }
            delete alignmentEntry;
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
    }

    void resetAlignment(GsnapAlignment *alignment) {
        alignment->lineNum = 0;
        alignment->pairedEnd = 0;
        alignment->pairType = NULL;
        alignment->referenceSequence = NULL;
        alignment->referenceQuality = NULL;
        alignment->queryIndex = 0;
        alignment->numAlignmentEntries = 0;
        clearAlignmentEntries(alignment->alignmentEntries);
        alignment->numAlignmentEntriesPair = 0;
        clearAlignmentEntries(alignment->alignmentEntriesPair);
        alignment->alignScore = 0;
        alignment->mapq = 0;
        alignment->pairScore = 0;
        alignment->insertLength = 0;
        alignment->pairSubType = NULL;
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
    int mapHasKey(LIBGOBY_HASH_MAP<string, char*> *whichMap, string toFind) {
        return (whichMap->find(toFind) != whichMap->end());
    }

    char *mapValueForKey(LIBGOBY_HASH_MAP<string, char*> *whichMap, char *toFind) {
        if (whichMap) {
            string keyToFind(toFind);
            if (mapHasKey(whichMap, keyToFind)) {
                return whichMap->at(keyToFind);
            }
        }
        return NULL;
    }

    int mapValueForKeyInt(LIBGOBY_HASH_MAP<string, char*> *whichMap, char *toFind, int defaultVal = 0) {
        if (whichMap) {
            string keyToFind(toFind);
            if (mapHasKey(whichMap, keyToFind)) {
                return atoi(whichMap->at(keyToFind));
            }
        }
        return defaultVal;
    }

    unsigned int mapValueForKeyUnsignedInt(LIBGOBY_HASH_MAP<string, char*> *whichMap, char *toFind, unsigned int defaultVal = 0) {
        if (whichMap) {
            string keyToFind(toFind);
            if (mapHasKey(whichMap, keyToFind)) {
                return (unsigned) strtoul(whichMap->at(keyToFind), NULL, 10);
            }
        }
        return defaultVal;
    }

    pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*> *mapOfValues(char *a) {
        LIBGOBY_HASH_MAP<string, char*> *result = new LIBGOBY_HASH_MAP<string, char*>();
        vector<char *> *keys = new vector<char *>();
        if (a != NULL) {
            vector<char *> *aParts = gobyGsnap_split(a, ",.");
            int size = aParts->size();
            for (int i = 0; i < size; i++) {
                vector<char *> *kvParts = gobyGsnap_split(aParts->at(i), ":");
                char *value = "1";
                if (kvParts->size() > 1) {
                    value = kvParts->at(1);
                }
                (*result)[kvParts->at(0)] = value;
                keys->push_back(kvParts->at(0));
                delete kvParts;
            }
            delete aParts;
        }
        return new pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*>(keys, result);
    }

    void parseEntryHeader(vector<char *> *lines, GsnapAlignment *alignment, int isPair) {
        char *line = lines->at(alignment->lineNum++);
        vector<char *> *parts = gobyGsnap_split(line, "\t");
        int numParts = parts->size();
        if (numParts < 3 || numParts > 4) {
            fprintf(stderr, "Unrecoverable error: Gsnap alignment header line contains more than <3 or >4 fields (%d).\n", numParts);
            exit(99);
        }
        char *refSeq = parts->at(0);
        refSeq++; // Remove the initial character
        vector<char *> *numEntriesParts = gobyGsnap_split(parts->at(1), " ");
        char *numEntriesStr = numEntriesParts->at(0);
        char *pairType = NULL;
        if (numEntriesParts->size() > 1) {
            pairType = numEntriesParts->at(1);
        }
        char *qual = NULL;
        char *queryIndexStr = NULL;
        if (numParts == 4) {
            qual = parts->at(2);
            queryIndexStr = parts->at(3);
        } else if (numParts == 3) {
            queryIndexStr = parts->at(2);
        }

        alignment->pairedEnd = pairType != NULL;
        alignment->pairType = pairType;
        alignment->referenceSequence = refSeq;
        alignment->referenceQuality = qual;
        alignment->queryIndex = (unsigned) strtoul(queryIndexStr, NULL, 10);
        if (isPair) {
            alignment->numAlignmentEntriesPair = atoi(numEntriesStr);
        } else {
            alignment->numAlignmentEntries = atoi(numEntriesStr);
        }
        delete numEntriesParts;
        delete parts;
    }

    void parseSegment(CAlignmentsWriterHelper *writerHelper,
            vector<char *> *lines, vector<char *> *currentParts,
            GsnapAlignment *alignment, GsnapAlignmentEntry *alignmentEntry) {
        vector<char *> *parts;
        if (currentParts != NULL) {
            parts = currentParts;
        } else {
            char *line = lines->at(alignment->lineNum++);
            parts = gobyGsnap_split(line, "\t");
        }

        char *querySeq = parts->at(0);
        querySeq++; // Skip first char of query sequence

        vector<char *> *readStartEndParts = gobyGsnap_split(parts->at(1), ".");
        char *readStartStr = readStartEndParts->at(0);
        char *readEndStr = readStartEndParts->at(1);
        delete readStartEndParts;
        
        vector<char *> *targetParts = gobyGsnap_split(parts->at(2), ":.");
        char *targetIdentifier = targetParts->at(0);
        int reverseStrand = 0;
        if (targetIdentifier[0] == '+') {
            reverseStrand = 0;
            targetIdentifier++;
        } else if (targetIdentifier[0] == '-') {
            reverseStrand = 1;
            targetIdentifier++;
        }
        char *targetStartStr = targetParts->at(1);
        char *targetEndStr = targetParts->at(2);
        delete targetParts;
        
        pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*> *settingsMapTuple = NULL;
        vector<char *>* settingsMapKeys = NULL;
        LIBGOBY_HASH_MAP<string, char*>* settingsMap = NULL;
        settingsMapTuple = mapOfValues(parts->at(3));
        settingsMapKeys = settingsMapTuple->first;
        settingsMap = settingsMapTuple->second;

        GsnapAlignmentSegment *newSegment = new GsnapAlignmentSegment;
        initializeAlignmentSegment(newSegment);
        alignmentEntry->alignmentSegments->push_back(newSegment);

        newSegment->readSequence = querySeq;
        newSegment->readStart = (unsigned) strtoul(readStartStr, NULL, 10);
        newSegment->readEnd = (unsigned) strtoul(readEndStr, NULL, 10);
        newSegment->reverseStrand = reverseStrand;
        newSegment->targetIdentifier = targetIdentifier;
        newSegment->targetIndex = gobyAlignments_targetIndexForIdentifier(writerHelper, targetIdentifier);
        newSegment->startType = settingsMapKeys->at(0);
        newSegment->startClip = mapValueForKeyInt(settingsMap, settingsMapKeys->at(0));
        newSegment->endType = settingsMapKeys->at(1);
        newSegment->endClip = mapValueForKeyInt(settingsMap, settingsMapKeys->at(1));
        newSegment->matches = mapValueForKeyInt(settingsMap, "matches");
        newSegment->subs = mapValueForKeyInt(settingsMap, "sub");
        
        delete settingsMapKeys;
        delete settingsMap;
        delete settingsMapTuple;
        if (currentParts == NULL) {
            delete parts;
        }
    }

    void parseEntries(CAlignmentsWriterHelper *writerHelper,
            vector<char *> *lines, GsnapAlignment *alignment, int isPair) {
        char *line = lines->at(alignment->lineNum++);
        vector<char *> *parts = gobyGsnap_split(line, "\t");
        int numParts = parts->size();
        if (numParts < 5) {
            fprintf(stderr, "Unrecoverable error: Gsnap alignment entry (first) line contains more than <5 (%d).\n", numParts);
            exit(99);
        }
        pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*> *settingsMapTuple = NULL;
        vector<char *>* settingsMapKeys = NULL;
        LIBGOBY_HASH_MAP<string, char*>* settingsMap = NULL;
        settingsMapTuple = mapOfValues(parts->at(4));
        settingsMapKeys = settingsMapTuple->first;
        settingsMap = settingsMapTuple->second;

        pair<vector<char *>*,LIBGOBY_HASH_MAP<string, char*>*> *pairSettingsMapTuple = NULL;
        vector<char *>* pairSettingsMapKeys = NULL;
        LIBGOBY_HASH_MAP<string, char*>* pairSettingsMap = NULL;
        if (numParts >= 6) {
            pairSettingsMapTuple = mapOfValues(parts->at(5));
            pairSettingsMapKeys = pairSettingsMapTuple->first;
            pairSettingsMap = pairSettingsMapTuple->second;
        }
        
        if (!isPair) {
            // Don't set this if we are parsing the second in a pair
            alignment->alignScore = mapValueForKeyInt(settingsMap, "align_score");
            alignment->mapq = mapValueForKeyInt(settingsMap, "mapq");
            alignment->pairScore = mapValueForKeyInt(pairSettingsMap, "pair_score");
            alignment->insertLength = mapValueForKeyInt(pairSettingsMap, "insert_length");
            alignment->pairSubType = mapValueForKey(pairSettingsMap, "pairtype");
        }
        
        int numSegs = mapValueForKeyInt(settingsMap, "segs");
        GsnapAlignmentEntry *newAlignmentEntry = new GsnapAlignmentEntry;
        if (isPair) {
            alignment->alignmentEntriesPair->push_back(newAlignmentEntry);
        } else {
            alignment->alignmentEntries->push_back(newAlignmentEntry);
        }
        for (int i = 0; i < numSegs; i++) {
            if (i == 0) {
                newAlignmentEntry->alignmentSegments = new vector<GsnapAlignmentSegment*>();
                parseSegment(writerHelper, lines, parts, alignment, newAlignmentEntry);
            } else {
                parseSegment(writerHelper, lines, NULL, alignment, newAlignmentEntry);
            }
        }
        delete settingsMapKeys;
        delete settingsMap;
        delete settingsMapTuple;
        if (pairSettingsMapTuple != NULL) {
            delete pairSettingsMapKeys;
            delete pairSettingsMap;
            delete pairSettingsMapTuple;
        }
        delete parts;
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

        parseEntryHeader(lines, alignment, 0);
        int numEntries;
        for (numEntries = 0; numEntries < alignment->numAlignmentEntries; numEntries++) {
            parseEntries(writerHelper, lines, alignment, 0);
        }
        if (alignment->pairedEnd) {
            // If a paired-end alignment
            parseEntryHeader(lines, alignment, 1);
            for (numEntries = 0; numEntries < alignment->numAlignmentEntriesPair; numEntries++) {
                parseEntries(writerHelper, lines, alignment, 1);
            }
        }
        delete lines;
    }
}
