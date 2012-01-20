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

    /**
     * Take a string and a set of one or more delimiters and split the string
     * into a vector of strings (char*'s). This method allocates no new
     * memory (except for the vector). As this method uses strtok,
     * IS DESTRUCTIVE TO THE STRING (it inserts '\0' in place of the)
     * delimieters. Repeated delimiters ARE IGNORED, such as if you
     * want to parse "a::b::c:d:e" delimited by ":" this will return
     * ["a","b","c","d","e"].
     * @param writerHelper the Goby writer helper
     * @param targetLines the lines of targets (one per line) to register.
     */
    vector<char *> *gobyGsnap_split(char *str, const char *delim) {
        char* token = strtok(str, delim);
        vector<char *> *result = new vector<char *>();
        while(token != NULL) {
            result->push_back(token);
            token = strtok(NULL, delim);
        }
        return result;
    }

    /**
     * A test method that will register a set of targets before a test
     * alignment.
     * @param writerHelper the Goby writer helper
     * @param targetLines the lines of targets (one per line) to register.
     */
    void gobyGsnap_test_registerTargets(CAlignmentsWriterHelper *writerHelper,
            char *targetLines) {
        vector<char *> *targets = gobyGsnap_split(targetLines, "\n\r");
        int numTargets = targets->size();
        for (int i = 0; i < numTargets; i++) {
            gobyAlignments_addTarget(writerHelper, i, targets->at(i),
                    2000000000);
        }
        delete targets;
    }

    /**
     * Delete existing alignment segments that might exist from a previous
     * alignment.
     * @param alignmentSegments the vector of alignment segments to clear out.
     */
    void clearAlignmentSegments(
            vector<GsnapAlignmentSegment*> *alignmentSegments)  {
        int size = alignmentSegments->size();
        for (int i = 0; i < size; i++) {
            delete alignmentSegments->at(i);
        }
    }

    /**
     * Delete existing alignment entries that might exist from a previous
     * alignment.
     * @param alignmentEntries the vector of alignment entries to clear out.
     */
    void clearAlignmentEntries(
            vector<GsnapAlignmentEntry*> *alignmentEntries)  {
        int size = alignmentEntries->size();
        for (int i = 0; i < size; i++) {
            GsnapAlignmentEntry *alignmentEntry = alignmentEntries->at(i);
            vector<GsnapAlignmentSegment*> *alignmentSegments =
                    alignmentEntry->alignmentSegments;
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

    /**
     * When creating a new alignment segment, this should be called to
     * initialize the values in the data structure.
     * @param alignmentSegment The Gsnap alignment segment data structure
     * to initialize.
     */
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

    /**
     * Before parsing a new Gsnap alignment entry, this needs to be called to
     * make sure all data from the previous alignment is cleared out.
     * @param alignment The Gsnap data structure used to parse Gsnap alignments
     */
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

    /**
     * Initialize the Gsnap data structures. This needs to be done once when we
     * first start parsing Gsnap data.
     * @param alignment The Gsnap data structure used to parse Gsnap alignments
     */
    void initializeAlignment(GsnapAlignment *alignment) {
        alignment->alignmentEntries = new vector<GsnapAlignmentEntry*>;
        alignment->alignmentEntriesPair = new vector<GsnapAlignmentEntry*>;
        resetAlignment(alignment);
    }

    /**
     * Destory the Gsnap parsing data structures if they exist. This should
     * be done after the alignment has been completely written, as the program
     * exits.
     * @param writerHelper the Goby writer helper
     */
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
     * Given a map created by createMap(), return 0 if the key isn't found
     * in the map or not zero (probably 1?) if the key is found.
     * @param whichMap The map to obtain the value for
     * @param toFind The key to lookup in the map
     */
    int mapHasKey(LIBGOBY_HASH_MAP<string, char*> *whichMap, string toFind) {
        return (whichMap->find(toFind) != whichMap->end());
    }

    /**
     * Given a map created by createMap(), return the value for the key
     * toFind as a string (char*). If the key isn't found, NULL is returned.
     * @param whichMap The map to obtain the value for
     * @param toFind The key to lookup in the map
     * @param defaultVal the default value to use if the key toFind isn't found.
     * @return the value for the key toFind (or defaultVal)
     */
    char *mapValueForKey(LIBGOBY_HASH_MAP<string, char*> *whichMap,
            char *toFind) {
        if (whichMap) {
            string keyToFind(toFind);
            if (mapHasKey(whichMap, keyToFind)) {
                return whichMap->at(keyToFind);
            }
        }
        return NULL;
    }

    /**
     * Given a map created by createMap(), return the value for the key
     * toFind as an int. If the key isn't found, defaultVal is returned.
     * @param whichMap The map to obtain the value for
     * @param toFind The key to lookup in the map
     * @param defaultVal the default value to use if the key toFind isn't found.
     * @return the value for the key toFind (or defaultVal)
     */
    int mapValueForKeyInt(LIBGOBY_HASH_MAP<string, char*> *whichMap,
            char *toFind, int defaultVal = 0) {
        if (whichMap) {
            string keyToFind(toFind);
            if (mapHasKey(whichMap, keyToFind)) {
                return atoi(whichMap->at(keyToFind));
            }
        }
        return defaultVal;
    }

    /**
     * Given a map created by createMap(), return the value for the key
     * toFind as an unsigned int. If the key isn't found, defaultVal is
     * returned.
     * @param whichMap The map to obtain the value for
     * @param toFind The key to lookup in the map
     * @param defaultVal the default value to use if the key toFind isn't found.
     * @return the value for the key toFind (or defaultVal)
     */
    unsigned int mapValueForKeyUnsignedInt(
            LIBGOBY_HASH_MAP<string, char*> *whichMap, char *toFind,
            unsigned int defaultVal = 0) {
        if (whichMap) {
            string keyToFind(toFind);
            if (mapHasKey(whichMap, keyToFind)) {
                return (unsigned) strtoul(whichMap->at(keyToFind), NULL, 10);
            }
        }
        return defaultVal;
    }

    /**
     * Parses key:value from the input string (a). key-values are
     * deparated by "," or ".". and ":" separates keys from values.
     * Returns a pair of values:
     * - The list of keys, in the order they are found in the input
     * - A map of key to value.
     * Calling createMap("ins:1..end:0,matches:23,sub:0") will return
     * [  ["ins", "end", "matches", "sub"],
     *    ["ins" : "1", "end" : "0", "matches" : "23", "sub" : "0"]  ]
     * @param a the input string that is being split into the map
     * @return the pair of (keys in parsed order, map of the key/valu pairs)
     */
    pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*> *createMap(
            char *a) {
        LIBGOBY_HASH_MAP<string, char*> *result =
            new LIBGOBY_HASH_MAP<string, char*>();
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
        return new pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*>(
                keys, result);
    }

    /**
     * Parse one of the 1 or more segments of an alignment entry.
     * @param writerHelper the Goby writer helper
     * @param parts A vector storing the parts of the alignment entry
     * segment already split by the tab character.
     * @param alignmentEntry the alignmentEntry that this new segment
     * belongs to.
     */
    void parseSegment(CAlignmentsWriterHelper *writerHelper,
            vector<char *> *parts, GsnapAlignmentEntry *alignmentEntry) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;

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
        
        pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*>
            *settingsMapTuple = NULL;
        vector<char *>* settingsMapKeys = NULL;
        LIBGOBY_HASH_MAP<string, char*>* settingsMap = NULL;
        settingsMapTuple = createMap(parts->at(3));
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
        newSegment->targetIndex = gobyAlignments_targetIndexForIdentifier(
                writerHelper, targetIdentifier);
        newSegment->startType = settingsMapKeys->at(0);
        newSegment->startClip = mapValueForKeyInt(
                settingsMap, settingsMapKeys->at(0));
        newSegment->endType = settingsMapKeys->at(1);
        newSegment->endClip = mapValueForKeyInt(
                settingsMap, settingsMapKeys->at(1));
        newSegment->matches = mapValueForKeyInt(settingsMap, "matches");
        newSegment->subs = mapValueForKeyInt(settingsMap, "sub");
        
        delete settingsMapKeys;
        delete settingsMap;
        delete settingsMapTuple;
    }

    /**
     * Parse the entries for an alignment (for a single side of the pair if
     * the alignment is paired end).
     * @param writerHelper the Goby writer helper
     * @param lines the vector which contains the lines of the alignment.
     * @param isPair true (1) if we are parsing entries for the paired end
     * portion of the alignment, otherwise 0.
     */
    void parseEntries(CAlignmentsWriterHelper *writerHelper,
            vector<char *> *lines, int isPair) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        char *line = lines->at(alignment->lineNum++);
        vector<char *> *parts = gobyGsnap_split(line, "\t");
        int numParts = parts->size();
        if (numParts < 5) {
            fprintf(stderr, "Unrecoverable error: Gsnap alignment entry (first) line contains more than <5 (%d).\n", numParts);
            exit(99);
        }
        pair<vector<char *>*, LIBGOBY_HASH_MAP<string, char*>*>
            *settingsMapTuple = NULL;
        vector<char *>* settingsMapKeys = NULL;
        LIBGOBY_HASH_MAP<string, char*>* settingsMap = NULL;
        settingsMapTuple = createMap(parts->at(4));
        settingsMapKeys = settingsMapTuple->first;
        settingsMap = settingsMapTuple->second;

        pair<vector<char *>*,LIBGOBY_HASH_MAP<string, char*>*>
            *pairSettingsMapTuple = NULL;
        vector<char *>* pairSettingsMapKeys = NULL;
        LIBGOBY_HASH_MAP<string, char*>* pairSettingsMap = NULL;
        if (numParts >= 6) {
            pairSettingsMapTuple = createMap(parts->at(5));
            pairSettingsMapKeys = pairSettingsMapTuple->first;
            pairSettingsMap = pairSettingsMapTuple->second;
        }
        
        if (!isPair) {
            // Don't set this if we are parsing the second in a pair
            alignment->alignScore = mapValueForKeyInt(
                    settingsMap, "align_score");
            alignment->mapq = mapValueForKeyInt(settingsMap, "mapq");
            alignment->pairScore = mapValueForKeyInt(
                    pairSettingsMap, "pair_score");
            alignment->insertLength = mapValueForKeyInt(
                    pairSettingsMap, "insert_length");
            alignment->pairSubType = mapValueForKey(
                    pairSettingsMap, "pairtype");
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
                newAlignmentEntry->alignmentSegments =
                        new vector<GsnapAlignmentSegment*>();
            } else {
                line = lines->at(alignment->lineNum++);
                parts = gobyGsnap_split(line, "\t");
            }
            parseSegment(writerHelper, parts, newAlignmentEntry);
            delete parts;
        }
        delete settingsMapKeys;
        delete settingsMap;
        delete settingsMapTuple;
        if (pairSettingsMapTuple != NULL) {
            delete pairSettingsMapKeys;
            delete pairSettingsMap;
            delete pairSettingsMapTuple;
        }
    }

    /**
     * Parse the header line for an alignment.
     * @param writerHelper the Goby writer helper
     * @param lines the vector which contains the lines of the alignment.
     * @param isPair true (1) if we are parsing a header entry for the
     * paired end portion of the alignment, otherwise 0.
     */
    void parseEntryHeader(CAlignmentsWriterHelper *writerHelper,
            vector<char *> *lines, int isPair) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
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

    /**
     * Parse a single Gsnap alignment into Gsnap data structures.
     * @param writerHelper the Goby writer helper
     * @param alignmentStr the lines that make up a single Gsnap alignment
     */
    void parseToGsnapDataStructures(CAlignmentsWriterHelper *writerHelper,
            char *alignmentStr) {
        vector<char *> *lines = gobyGsnap_split(alignmentStr, "\n\r");
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        if (alignment == NULL) {
            alignment = new GsnapAlignment;
            writerHelper->gsnapAlignment = alignment;
            initializeAlignment(alignment);
        } else {
            resetAlignment(alignment);
        }

        parseEntryHeader(writerHelper, lines, 0);
        int numEntries;
        for (numEntries = 0;
                numEntries < alignment->numAlignmentEntries;
                numEntries++) {
            parseEntries(writerHelper, lines, 0);
        }
        if (alignment->pairedEnd) {
            // If a paired-end alignment
            parseEntryHeader(writerHelper, lines, 1);
            for (numEntries = 0;
                    numEntries < alignment->numAlignmentEntriesPair;
                    numEntries++) {
                parseEntries(writerHelper, lines, 1);
            }
        }
        delete lines;
    }
    
    /**
     * Write the Gsnap alignment from the Gsnap data structures into Goby
     * compact-alignment format.
     * @param writerHelper the Goby writer helper
     * @param alignmentStr the lines that make up a single Gsnap alignment
     */
    void writeGobyAlignment(CAlignmentsWriterHelper *writerHelper) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
    }
    
    /**
     * Parse and save a single gsnap alignment entry. Supports single or
     * paired end. Supports spliced alighments.
     * @param writerHelper the alignment writer helper
     * @param alignment the gsnap alignment entry to parse. Supports single
     * or paired end.
     * This should be ONE alignment entry (can be two stanzas if paired end)
     * @param writerHelper the Goby writer helper
     * @param alignmentStr the lines that make up a single Gsnap alignment
     */
    void gobyGsnap_parse(CAlignmentsWriterHelper *writerHelper,
            char *alignmentStr) {
        parseToGsnapDataStructures(writerHelper, alignmentStr);
        writeGobyAlignment(writerHelper);
    }
}
