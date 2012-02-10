#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#include "C_Gsnap.h"
#include "C_Gsnap_structs.h"
#include "C_Alignments.h"
#include "hash.h"
#include "pcrecpp.h"

/*
 * TODO: Make sure the fields I moved from segment to alignment really
 *       DO belong in alignment.
 * TODO: Support paired writing. This might involve supporting more than
 *       one pair target in the proto.
 * TODO: Support spliced writing. Probably need to grab more gsnap fields.
 * TODO: Make a pool of GsnapAlignmentSegments. Check out from pool for
 *       each alignment, clean then, and then use. The return them to the
 *       pool. This should be a difference in parse speed. If none are in
 *       the pool, create a new one. This also means don't check for NULL
 *       strings, check for empty strings. 
 */

using namespace std;

using pcrecpp::StringPiece;
using pcrecpp::RE;
using pcrecpp::RE_Options;

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

    static char *complCode = "???????????????????????????????? ??#$%&')(*+,-./0123456789:;>=<??TVGHEFCDIJMLKNOPQYSAABWXRZ]?[^_`tvghefcdijmlknopqysaabwxrz}|{~?";

    /**
     * Take a string and a set of one or more delimiters and split the string
     * into a vector of strings (char*'s). This method allocates no new
     * memory (except for the vector). As this method uses strtok,
     * IS DESTRUCTIVE TO THE STRING -- it inserts '\0' in place of the
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

    void destroyAlignmentSegment(GsnapAlignmentSegment *seg) {
        if (seg->segmentSequence != NULL) {
            delete seg->segmentSequence;
        }
        if (seg->querySequence != NULL) {
            delete seg->querySequence;
        }
        if (seg->queryQuality != NULL) {
            delete seg->queryQuality;
        }
        if (seg->referenceSequence != NULL) {
            delete seg->referenceSequence;
        }
        if (seg->deletesSequence != NULL) {
            delete seg->deletesSequence;
        }
        if (seg->insertsSequence != NULL) {
            delete seg->insertsSequence;
        }
        delete seg;
    }

    /**
     * Delete existing alignment segments that might exist from a previous
     * alignment.
     * @param alignmentSegments the vector of alignment segments to clear out.
     */
    void destroyAlignmentSegments(
            vector<GsnapAlignmentSegment*> *segs)  {
        int size = segs->size();
        for (int i = 0; i < size; i++) {
            destroyAlignmentSegment(segs->at(i));
        }
        delete segs;
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
                destroyAlignmentSegments(alignmentSegments);
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
        alignmentSegment->segmentSequence = NULL;
        alignmentSegment->referenceSequence = NULL;
        alignmentSegment->querySequence = NULL;
        alignmentSegment->queryQuality = NULL;
        alignmentSegment->deletesSequence = NULL;
        alignmentSegment->insertsSequence = NULL;
        alignmentSegment->queryStart = 0;
        alignmentSegment->queryEnd = 0;
        alignmentSegment->reverseStrand = 0;
        alignmentSegment->targetIdentifier = NULL;
        alignmentSegment->targetIndex = 0;
        alignmentSegment->targetStart = 0;
        alignmentSegment->targetEnd = 0;
        alignmentSegment->startType = STARTENDTYPE_UNKNOWN;
        alignmentSegment->startClip = 0;
        alignmentSegment->startProb = 0.0;
        alignmentSegment->endType = STARTENDTYPE_UNKNOWN;
        alignmentSegment->endClip = 0;
        alignmentSegment->endProb = 0.0;
        alignmentSegment->matches = 0;
        alignmentSegment->subs = 0;
        alignmentSegment->inserts = 0;
        alignmentSegment->deletes = 0;
        alignmentSegment->spliceDir = SPLICEDIR_UNKNOWN;
        alignmentSegment->spliceType = SPLICETYPE_UNKNOWN;

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
        alignment->querySequence->resize(0);
        alignment->queryQuality->resize(0);
        alignment->queryIndex = 0;
        alignment->numAlignmentEntries = 0;
        clearAlignmentEntries(alignment->alignmentEntries);
        alignment->numAlignmentEntriesPair = 0;
        clearAlignmentEntries(alignment->alignmentEntriesPair);
        alignment->alignScore = 0;
        alignment->mapq = 0;
        alignment->pairScore = 0;
        alignment->insertLength = 0;
        if (alignment->pairSubType != NULL) {
            delete alignment->pairSubType;
            alignment->pairSubType = NULL;
        }
    }

    /**
     * Initialize the Gsnap data structures. This needs to be done once when we
     * first start parsing Gsnap data.
     * @param alignment The Gsnap data structure used to parse Gsnap alignments
     */
    void initializeAlignment(GsnapAlignment *alignment) {
        alignment->alignmentEntries = new vector<GsnapAlignmentEntry*>;
        alignment->alignmentEntriesPair = new vector<GsnapAlignmentEntry*>;
        alignment->querySequence = new string("");
        alignment->queryQuality = new string("");
        alignment->pairSubType = NULL;
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
            delete alignment->querySequence;
            delete alignment->queryQuality;
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
    bool mapHasKey(LIBGOBY_HASH_MAP<string, string> whichMap, string toFind) {
        return (whichMap.find(toFind) != whichMap.end());
    }

    /**
     * Given a map created by createMap(), return the value for the key
     * toFind as a string*. If the key isn't found, NULL is returned.
     * IF NOT NULL YOU MUST DELETE THIS VALUE WHEN YOU ARE DONE WITH IT.
     * @param whichMap The map to obtain the value for
     * @param toFind The key to lookup in the map
     * @param defaultVal the default value to use if the key toFind isn't found.
     * @return the value for the key toFind (or defaultVal)
     */
    string *mapValueForKey(
            LIBGOBY_HASH_MAP<string, string> whichMap, string toFind) {
        string keyToFind(toFind);
        if (mapHasKey(whichMap, keyToFind)) {
            return new string(whichMap[keyToFind].c_str());
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
    long mapValueForKeyLong(
            LIBGOBY_HASH_MAP<string, string> whichMap, string toFind) {
        string keyToFind(toFind);
        if (mapHasKey(whichMap, keyToFind)) {
            return atoi(whichMap[keyToFind].c_str());
        }
        return 0;
    }

        /**
     * Given a map created by createMap(), return the value for the key
     * toFind as an int. If the key isn't found, defaultVal is returned.
     * @param whichMap The map to obtain the value for
     * @param toFind The key to lookup in the map
     * @param defaultVal the default value to use if the key toFind isn't found.
     * @return the value for the key toFind (or defaultVal)
     */
    double mapValueForKeyDouble(
            LIBGOBY_HASH_MAP<string, string> whichMap, string toFind) {
        string keyToFind(toFind);
        double value;
        if (mapHasKey(whichMap, keyToFind)) {
            istringstream ss(whichMap[keyToFind]);
            ss >> value;
        } else {
            value = 0.0;
        }
        return value;
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
        LIBGOBY_HASH_MAP<string, string> whichMap, const char *toFind) {
        string keyToFind(toFind);
        if (mapHasKey(whichMap, keyToFind)) {
            return (unsigned) strtoul(whichMap.at(keyToFind).c_str(), NULL, 10);
        }
        return 0;
    }

    pair<vector<string>, LIBGOBY_HASH_MAP<string, string> > createMap(char *a) {
        LIBGOBY_HASH_MAP<string, string> result;
        vector<string> keys;

        if (a != NULL) {
        
            string temp1;
            string da;     string daClip;   string daProb;
            string temp2;
            string std;    string stdClip;

            pcrecpp::StringPiece input(a);  // Wrap in a StringPiece
            pcrecpp::RE re("(([A-Za-z_\\-]+):(\\d+)\\.(\\d+))|(([A-Za-z_\\-]+):(\\d+))");

            while (re.FindAndConsume(&input, &temp1,
                                     &da, &daClip, &daProb,
                                     &temp2,
                                     &std, &stdClip)) {
                if (da.length() > 0) {
                    if (da.length() > 12 && da.substr(0, 12) == "splice_dist_") {
                        // Remove the _1, _2 suffix
                        da.resize(11);
                    }
                    keys.push_back(da);
                    if (daProb.length() > 0) {
                        result[da] = "0";
                        result[da + "-prob"] = daClip + "." + daProb;
                    } else {
                        result[da] = daClip;
                    }
                } else {
                    keys.push_back(std);
                    result[std] = stdClip;
                }
            }
        }
        return pair<vector<string>,
                    LIBGOBY_HASH_MAP<string, string> >(keys, result);
    }

    SegmentStartEndType segmentStartEndTypeFromString(string type) {
        //START, END, INS, DEL, TERM, DONOR, ACCEPTOR, UNKNOWN
        if (type.length() == 0) {
            return STARTENDTYPE_UNKNOWN;
        } else if (type == "start") {
            return STARTENDTYPE_START;
        } else if (type == "end") {
            return STARTENDTYPE_END;
        } else if (type == "ins") {
            return STARTENDTYPE_INS;
        } else if (type == "del") {
            return STARTENDTYPE_DEL;
        } else if (type == "term") {
            return STARTENDTYPE_TERM;
        } else if (type == "donor") {
            return STARTENDTYPE_DONOR;
        } else if (type == "acceptor") {
            return STARTENDTYPE_ACCEPTOR;
        } else {
            return STARTENDTYPE_UNKNOWN;
        }
    }

    SegmentSpliceDir segmentSpliceDirFromString(string *type) {
        //"sense", "antisense"
        // NONE, SENSE, ANTISENSE, UNKNOWN
        if (type == NULL || type->length() == 0) {
            return SPLICEDIR_NONE;
        } else if ((*type) == "sense") {
            return SPLICEDIR_SENSE;
        } else if ((*type) == "antisense") {
            return SPLICEDIR_ANTISENSE;
        } else {
            return SPLICEDIR_UNKNOWN;
        }
    }

    SegmentSpliceType segmentSpliceTypeFromString(string *type) {
        // "consistent", "inversion", "scramble", "translocation"
        //  NONE, CONSISTENT, INVERSION, SCRAMBLE, TRANSLOCATION, UNKNOWN };
        if (type == NULL || type->length() == 0) {
            return SPLICETYPE_NONE;
        } else if ((*type) == "consistent") {
            return SPLICETYPE_CONSISTENT;
        } else if ((*type) == "inversion") {
            return SPLICETYPE_INVERSION;
        } else if ((*type) == "scramble") {
            return SPLICETYPE_SCRAMBLE;
        } else if ((*type) == "translocation") {
            return SPLICETYPE_TRANSLOCATION;
        } else {
            return SPLICETYPE_UNKNOWN;
        }
    }

    string *reverseCharArray(const char *toReverse, int size, bool complement) {
        if (toReverse == NULL) {
            return NULL;
        }
        string *result = new string("");
        for (int i = size - 1; i >= 0; i--) {
            const char current = toReverse[i];
            result->append(1, complement ? complCode[current] : current);
        }
        return result;
    }

    string *reverseString(string *toReverse, bool complement) {
        return reverseCharArray(toReverse->c_str(), toReverse->size(), complement);
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

        vector<char *> *queryStartEndParts = gobyGsnap_split(parts->at(1), ".");
        char *queryStartStr = queryStartEndParts->at(0);
        char *queryEndStr = queryStartEndParts->at(1);
        delete queryStartEndParts;
        
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

        char *segmentSeqBuf = parts->at(0);
        segmentSeqBuf++; // Skip first char of query sequence
        string *segmentSeq = NULL;
        string *querySequence = NULL;
        string *queryQuality = NULL;
        if (reverseStrand) {
            segmentSeq = reverseCharArray(segmentSeqBuf, strlen(segmentSeqBuf), true);
            querySequence = reverseString(alignment->querySequence, true);
            if (alignment->queryQuality->size() > 0) {
                queryQuality = reverseString(alignment->queryQuality, false);
            }
        } else {
            segmentSeq = new string(segmentSeqBuf);
            querySequence = new string(*(alignment->querySequence));
            if (alignment->queryQuality->size() > 0) {
                queryQuality = new string(*(alignment->queryQuality));
            }
        }
        int queryLength = segmentSeq->size();
        
        pair<vector<string>, LIBGOBY_HASH_MAP<string, string> >
                settingsMapTuple = createMap(parts->at(3));
        vector<string> settingsMapKeys = settingsMapTuple.first;
        LIBGOBY_HASH_MAP<string, string> settingsMap = settingsMapTuple.second;

        GsnapAlignmentSegment *newSegment = new GsnapAlignmentSegment;
        initializeAlignmentSegment(newSegment);
        alignmentEntry->alignmentSegments->push_back(newSegment);
        
        newSegment->segmentSequence = segmentSeq;
        newSegment->querySequence = querySequence;
        newSegment->queryQuality = queryQuality;      
        
        newSegment->reverseStrand = reverseStrand;
        newSegment->targetIdentifier = targetIdentifier;
        newSegment->targetIndex = gobyAlignments_targetIndexForIdentifier(
                writerHelper, targetIdentifier);

        // Make targetStart, targetEnd 0-based
        unsigned targetStart = (unsigned) strtoul(targetStartStr, NULL, 10) - 1;
        unsigned targetEnd = (unsigned) strtoul(targetEndStr, NULL, 10) - 1;
        if (reverseStrand) {
            newSegment->targetStart = targetEnd;
            newSegment->targetEnd = targetStart;
        } else {
            newSegment->targetStart = targetStart;
            newSegment->targetEnd = targetEnd;
        }

        if (reverseStrand) {
            newSegment->startType =  segmentStartEndTypeFromString(settingsMapKeys.at(1));
            newSegment->startClip = mapValueForKeyLong(settingsMap, settingsMapKeys.at(1));
            newSegment->startProb = mapValueForKeyDouble(settingsMap, settingsMapKeys.at(1) + "-prob");
            newSegment->endType = segmentStartEndTypeFromString(settingsMapKeys.at(0));
            newSegment->endClip = mapValueForKeyLong(settingsMap, settingsMapKeys.at(0));
            newSegment->endProb = mapValueForKeyDouble(settingsMap, settingsMapKeys.at(0) + "-prob");
        } else {
            newSegment->startType = segmentStartEndTypeFromString(settingsMapKeys.at(0));
            newSegment->startClip = mapValueForKeyLong(settingsMap, settingsMapKeys.at(0));
            newSegment->startProb = mapValueForKeyDouble(settingsMap, settingsMapKeys.at(0) + "-prob");
            newSegment->endType =  segmentStartEndTypeFromString(settingsMapKeys.at(1));
            newSegment->endClip = mapValueForKeyLong(settingsMap, settingsMapKeys.at(1));
            newSegment->endProb = mapValueForKeyDouble(settingsMap, settingsMapKeys.at(1) + "-prob");
        }

        newSegment->spliceDir = segmentSpliceDirFromString(
                mapValueForKey(settingsMap, "dir"));
        newSegment->spliceType = segmentSpliceTypeFromString(
                mapValueForKey(settingsMap, "splice_type"));
        newSegment->spliceDistance = mapValueForKeyUnsignedInt(
                settingsMap, "splice_dist");
        
        // Change to 0-based
        int queryStart = atoi(queryStartStr) - 1;
        int queryEnd = atoi(queryEndStr) - 1;
        if (reverseStrand) {
            newSegment->queryStart = queryLength - queryEnd - 1;
            newSegment->queryEnd = queryLength - queryStart - 1;
            if (newSegment->startType == STARTENDTYPE_DEL) {
                newSegment->deletesSequence = new string(
                        newSegment->segmentSequence->substr(
                            newSegment->queryStart - newSegment->startClip, newSegment->startClip));
                newSegment->segmentSequence->erase(
                        newSegment->queryStart - newSegment->startClip, newSegment->startClip);
                newSegment->queryStart -= newSegment->startClip;
                newSegment->queryEnd -= newSegment->startClip;
            } else if (newSegment->startType == STARTENDTYPE_INS) {
                newSegment->insertsSequence = new string(
                        newSegment->querySequence->substr(
                            newSegment->queryStart - newSegment->startClip, newSegment->startClip));
            }
        } else {
            newSegment->queryStart = queryStart;
            newSegment->queryEnd = queryEnd;
            if (newSegment->endType == STARTENDTYPE_DEL) {
                newSegment->deletesSequence = new string(
                        newSegment->segmentSequence->substr(
                            newSegment->queryEnd + 1, newSegment->endClip));
                newSegment->segmentSequence->erase(
                        newSegment->queryEnd + 1, newSegment->endClip);
            } else if (newSegment->endType == STARTENDTYPE_INS) {
                newSegment->insertsSequence = new string(
                        newSegment->querySequence->substr(
                            newSegment->queryEnd + 1, newSegment->endClip));
            }
        }

        newSegment->matches = mapValueForKeyLong(settingsMap, "matches");
        newSegment->subs = mapValueForKeyLong(settingsMap, "sub");
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
        debug(printf("%s\n", line););
        vector<char *> *parts = gobyGsnap_split(line, "\t");
        int numParts = parts->size();
        if (numParts < 5) {
            fprintf(stderr, "Unrecoverable error: Gsnap alignment entry (first) line contains < 5 parts (%d).\n", numParts);
            exit(99);
        }
        pair<vector<string>, LIBGOBY_HASH_MAP<string, string> > settingsMapTuple = createMap(parts->at(4));
        vector<string> settingsMapKey = settingsMapTuple.first;
        LIBGOBY_HASH_MAP<string, string> settingsMap = settingsMapTuple.second;

        pair<vector<string>,LIBGOBY_HASH_MAP<string, string> > pairSettingsMapTuple;
        vector<string> pairSettingsMapKeys;
        LIBGOBY_HASH_MAP<string, string> pairSettingsMap;
        if (numParts >= 6) {
            pairSettingsMapTuple = createMap(parts->at(5));
            pairSettingsMapKeys = pairSettingsMapTuple.first;
            pairSettingsMap = pairSettingsMapTuple.second;
        }
        
        if (!isPair) {
            // Don't set this if we are parsing the second in a pair
            alignment->alignScore = mapValueForKeyLong(settingsMap, "align_score");
            alignment->mapq = mapValueForKeyLong(settingsMap, "mapq");
            alignment->pairScore = mapValueForKeyLong(pairSettingsMap, "pair_score");
            alignment->insertLength = mapValueForKeyLong(pairSettingsMap, "insert_length");
            alignment->pairSubType = mapValueForKey(pairSettingsMap, "pairtype");
        }
        
        int numSegs = mapValueForKeyLong(settingsMap, "segs");
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
                debug(printf("%s\n", line););
                parts = gobyGsnap_split(line, "\t");
            }
            parseSegment(writerHelper, parts, newAlignmentEntry);
            delete parts;
        }

        if ((numSegs > 1) &&
                (newAlignmentEntry->alignmentSegments->at(0)->reverseStrand)) {
            reverse(newAlignmentEntry->alignmentSegments->begin(),
                    newAlignmentEntry->alignmentSegments->end());
            for (int segNum = 0; segNum < numSegs - 1; segNum++) {
                GsnapAlignmentSegment *cur =
                        newAlignmentEntry->alignmentSegments->at(segNum);
                GsnapAlignmentSegment *next =
                        newAlignmentEntry->alignmentSegments->at(segNum + 1);
                cur->deletesSequence = next->deletesSequence;
                next->deletesSequence = NULL;
                cur->insertsSequence = next->insertsSequence;
                next->insertsSequence = NULL;
                
            }
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
        debug(printf("%s\n", line););
        vector<char *> *parts = gobyGsnap_split(line, "\t");
        int numParts = parts->size();
        if (numParts < 3 || numParts > 4) {
            fprintf(stderr, "Unrecoverable error: Gsnap alignment header line contains more than <3 or >4 fields (%d).\n", numParts);
            exit(99);
        }
        char *querySeq = parts->at(0);
        querySeq++; // Remove the initial character
        vector<char *> *numEntriesParts = gobyGsnap_split(parts->at(1), " ");
        char *numEntriesStr = numEntriesParts->at(0);
        char *pairType = NULL;
        if (numEntriesParts->size() > 1) {
            if (numEntriesParts->at(1)[0] != '(') {
                pairType = numEntriesParts->at(1);
            }
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
        alignment->querySequence->resize(0);
        alignment->querySequence->append(querySeq);
        alignment->queryQuality->resize(0);
        if (qual != NULL) {
            int qualSize = strlen(qual);
            char convertedChar;
            for (int qualPos = 0; qualPos < qualSize; qualPos++) {
                // Convert from Gsnap's SANGER to Phred
                alignment->queryQuality->append(1, qual[qualPos]);
            }
        }
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
        // Throw these away so one isn't tempted to use them. The
        // versions you are looking for exist in Segment.
        delete lines;
    }

    bool contiguousSegment(GsnapAlignmentSegment *seg) {
        SegmentStartEndType startType = seg->startType;
        SegmentStartEndType endType = seg->endType;
        // Check of segments are contiguous
        return ((startType == STARTENDTYPE_START || 
                    startType == STARTENDTYPE_END ||
                    startType == STARTENDTYPE_INS ||
                    startType == STARTENDTYPE_DEL || 
                    startType == STARTENDTYPE_TERM)  &&
                (endType == STARTENDTYPE_START ||
                    endType == STARTENDTYPE_END ||
                    endType == STARTENDTYPE_INS || 
                    endType == STARTENDTYPE_DEL ||
                    endType == STARTENDTYPE_TERM));
    }

    bool splicedSegment(GsnapAlignmentSegment *seg) {
        SegmentStartEndType startType = seg->startType;
        SegmentStartEndType endType = seg->endType;
        // Check if segments are spliced
        return (startType == STARTENDTYPE_DONOR || 
                startType == STARTENDTYPE_ACCEPTOR ||
                endType == STARTENDTYPE_DONOR || 
                endType == STARTENDTYPE_ACCEPTOR);
    }

    void buildMerged(int startPos, int endPos,
            GsnapAlignmentSegment *merged, GsnapAlignmentSegment *seg,
            string *mergedRef, string *mergedQuery, string *mergedQuality) {
        for (int pos = startPos; pos < endPos; pos++) {
            char refChar = seg->segmentSequence->at(pos);
            if (refChar == '.') {
                refChar = 'C';
            }
            mergedRef->append(1, refChar);
            mergedQuery->append(1, merged->querySequence->at(pos));
            if (mergedQuality != NULL) {
                mergedQuality->append(1, merged->queryQuality->at(pos));
            }
        }
    }

    /**
     * When >= 1 segment are contiguous (all startType and endTypes are
     * STAR,END,INS,DEL,TERM) this should be called to merge them. This will
     * build the reference sequence. All of the work is done in the first
     * segment of segs and that segment will be returned, now this this
     * single segment is representative of all of the segments in segs.
     * This should be called even if these is just one segment.
     * @param segs the segements of the match. If the segments are reverseStrand,
     * seg->querySequence, seg->segmentSequence, seg->queryQuality have already
     * been reversed (complemented).
     * @return  the merged segments (the first of segs modified to be the
     * merged version).
     */
    GsnapAlignmentSegment *mergeSegments(
            vector<GsnapAlignmentSegment*> *segs, bool splicedSegment) {
        GsnapAlignmentSegment *merged = segs->at(0);
        int numSegs = segs->size();
        int origLength = merged->segmentSequence->size();
        string *mergedRef = new string("");
        string *mergedQuery = new string("");
        string *mergedQuality = NULL;
        if (merged->queryQuality != NULL) {
            mergedQuality = new string("");
        }
        bool reverseStrand = (merged->reverseStrand != 0);
        
        int deletes = 0;
        for (int segNum = 0; segNum < numSegs; segNum++) {
            GsnapAlignmentSegment *seg = segs->at(segNum);
            int startPos;
            int endPos;

            debug(
                cout << "c=";
                for (int i = 0; i < seg->segmentSequence->size(); i++) {
                    cout << (i % 10);
                }
                cout << endl;
                cout << "q=" << *(merged->querySequence) << endl;
                cout << "s=" << *(seg->segmentSequence) << endl;
            );

            if (segNum == 0) {
                if (splicedSegment) {
                    startPos = 0;
                } else {
                    startPos = seg->queryStart - seg->startClip;
                }
                endPos = seg->queryStart;
                debug(cout << "before" << endl;);
                buildMerged(startPos, endPos, merged, seg,
                    mergedRef, mergedQuery, mergedQuality);
            }

            startPos = seg->queryStart;
            endPos = seg->queryEnd + 1;
            debug(
                    cout << "----" << endl;
                    cout << "in" << endl;
            );
            buildMerged(startPos, endPos, merged, seg,
                mergedRef, mergedQuery, mergedQuality);

            debug(
                cout << "----" << endl;
                cout << "after" << endl;
            );
            if (seg->deletesSequence != NULL) {
                int numDeletes = seg->deletesSequence->size();
                mergedRef->append(*(seg->deletesSequence));
                mergedQuery->append(numDeletes, '-');
                if (mergedQuality != NULL) {
                    mergedQuality->append(numDeletes, GOBY_NO_QUAL);
                }
                deletes += numDeletes;
            } else if (seg->insertsSequence != NULL) {
                int numInserts = seg->insertsSequence->size();
                mergedRef->append(numInserts, '-');
                mergedQuery->append(*(seg->insertsSequence));
                if (mergedQuality != NULL) {
                    mergedQuality->append(merged->queryQuality->substr(
                        seg->queryEnd + 1, numInserts));
                }
            }

            if (segNum == numSegs - 1) {
                startPos = seg->queryEnd + 1;
                if (splicedSegment) {
                    // Gsnap splicing doesn't, at this time, support indels,
                    // so we don't have to worry about exending the end position.
                    endPos = origLength;
                } else {
                    endPos = startPos + seg->endClip;
                }
                buildMerged(startPos, endPos, merged, seg,
                    mergedRef, mergedQuery, mergedQuality);
                merged->queryEnd = seg->queryEnd + deletes;
                merged->endClip = seg->endClip;
                merged->endType = seg->endType;
                merged->endProb = seg->endProb;
            }
        }

        merged->referenceSequence = mergedRef;

        string *temp = merged->segmentSequence;
        merged->segmentSequence = mergedQuery;
        delete temp;

        temp = merged->queryQuality;
        merged->queryQuality = mergedQuality;
        if (temp != NULL) {
            delete temp;
        }

        debug(
            cout << "merged:" << endl;
            cout << "c=";
            for (int i = 0; i < merged->segmentSequence->size(); i++) {
                cout << (i % 10);
            }
            cout << endl;
            cout << "q=" << *(merged->querySequence) << endl;
            cout << "r=" << *(merged->referenceSequence) << endl;
            cout << "s=" << *(merged->segmentSequence) << endl;
        );

        
        return merged;
    }

    /**
     * Similar to merge segments above, this will build the reference, etc.
     * for a single segment, no merging required. 
     * @param seg the segment to process
     */
    void processSpliceSegment(GsnapAlignmentSegment *seg) {
        vector<GsnapAlignmentSegment*> tempVector;
        tempVector.push_back(seg);
        mergeSegments(&tempVector, true);
    }

    void outputSequenceVariations(
            CAlignmentsWriterHelper *writerHelper,
            GsnapAlignmentSegment *merged) {

        gobyAlignments_outputSequenceVariations(writerHelper,
                merged->referenceSequence->c_str(),
                merged->segmentSequence->c_str(),
                merged->queryQuality == NULL ? NULL : merged->queryQuality->c_str(),
                merged->queryStart, merged->queryEnd,
                merged->reverseStrand, -33 /*qualityShift*/,
                &(merged->matches), &(merged->subs),
                &(merged->inserts), &(merged->deletes));

        debug(
            cout << "matches=" << merged->matches << " " << 
                    "subs=" << merged->subs << " " << 
                    "inserts=" << merged->inserts << " " <<
                    "deletes=" << merged->deletes << endl;
        );
    }

    void writeAlignmentSegment(CAlignmentsWriterHelper *writerHelper,
            GsnapAlignmentSegment *mergedSeg, int fragmentIndex) {

        GsnapAlignment *alignment = writerHelper->gsnapAlignment;

        gobyAlignments_appendEntry(writerHelper);    // checked

        // This will output sequence variations && calculate actual subs, matches, index.
        outputSequenceVariations(writerHelper, mergedSeg);

        gobyAlEntry_setMultiplicity(writerHelper, 1);    // checked
        gobyAlEntry_setQueryIndex(writerHelper, alignment->queryIndex); // checked
        gobyAlEntry_setTargetIndex(writerHelper,
                mergedSeg->targetIndex);  //checked
        gobyAlEntry_setPosition(writerHelper,
                mergedSeg->targetStart - mergedSeg->startClip);  // Checked
        gobyAlEntry_setMatchingReverseStrand(writerHelper,
                mergedSeg->reverseStrand); // Checked
        gobyAlEntry_setQueryPosition(writerHelper, mergedSeg->startClip);  //
        gobyAlEntry_setScoreInt(writerHelper, mergedSeg->matches);
        gobyAlEntry_setNumberOfMismatches(writerHelper, mergedSeg->subs);
        gobyAlEntry_setNumberOfIndels(writerHelper,
                mergedSeg->inserts + mergedSeg->deletes);
        gobyAlEntry_setQueryAlignedLength(writerHelper,
                mergedSeg->referenceSequence->size() -
                mergedSeg->startClip - mergedSeg->endClip - 
                mergedSeg->deletes);
        gobyAlEntry_setTargetAlignedLength(writerHelper,
                mergedSeg->referenceSequence->size() -
                mergedSeg->startClip - mergedSeg->endClip - 
                mergedSeg->inserts);
        gobyAlEntry_setQueryLength(writerHelper,
                alignment->querySequence->size());
        gobyAlEntry_setFragmentIndex(writerHelper, fragmentIndex);
        gobyAlEntry_setMappingQuality(writerHelper, alignment->mapq);
    }
    
    void writePairAlignment(CAlignmentsWriterHelper *writerHelper) {
        //
        // TODO: FINISH THIS. Currently implementation isn't right.
        //
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        int size = alignment->alignmentEntries->size();
        for (int i = 0; i < size; i++) {
            GsnapAlignmentEntry *alignmentEntry =
                    alignment->alignmentEntries->at(i);
            if (alignmentEntry->alignmentSegments == NULL ||
                    alignmentEntry->alignmentSegments->size() == 0) {
                continue;
            }
            
            if (splicedSegment(alignmentEntry->alignmentSegments->at(0))) {
                // Splices not yet supported
                continue;
            } else {
                GsnapAlignmentSegment *mergedSeg = mergeSegments(
                        alignmentEntry->alignmentSegments, false);

                writeAlignmentSegment(writerHelper, mergedSeg, 0);
            }
        }

        size = alignment->alignmentEntriesPair->size();
        for (int i = 0; i < size; i++) {
            GsnapAlignmentEntry *alignmentEntry =
                    alignment->alignmentEntries->at(i);
            if (alignmentEntry->alignmentSegments == NULL ||
                    alignmentEntry->alignmentSegments->size() == 0) {
                continue;
            }
            
            if (splicedSegment(alignmentEntry->alignmentSegments->at(0))) {
                // Splices not yet supported
                continue;
            }
            
            GsnapAlignmentSegment *mergedSeg = mergeSegments(
                    alignmentEntry->alignmentSegments, false);

            writeAlignmentSegment(writerHelper, mergedSeg, 1);
        }
}

    void writeSingleAlignment(CAlignmentsWriterHelper *writerHelper) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        int size = alignment->alignmentEntries->size();
        for (int i = 0; i < size; i++) {
            GsnapAlignmentEntry *alignmentEntry =
                    alignment->alignmentEntries->at(i);
            if (alignmentEntry->alignmentSegments == NULL ||
                    alignmentEntry->alignmentSegments->size() == 0) {
                continue;
            }
            
            bool splicedAlignment = false;
            if (splicedSegment(alignmentEntry->alignmentSegments->at(0))) {
                
                splicedAlignment = true;
                    
                if (alignmentEntry->alignmentSegments->size() > 1)) {
                    // Splice but not a half splice... We'll output half
                    // splices as non-splices, see below.

                    GsnapAlignmentSegment *firstSegment;
                    GsnapAlignmentSegment *secondSegment;
                    if (alignmentEntry->alignmentSegments->at(0)->targetStart < 
                            alignmentEntry->alignmentSegments->at(1)->targetStart) {
                        firstSegment = alignmentEntry->alignmentSegments->at(0);
                        secondSegment = alignmentEntry->alignmentSegments->at(1);
                    } else {
                        firstSegment = alignmentEntry->alignmentSegments->at(1);
                        secondSegment = alignmentEntry->alignmentSegments->at(0);
                    }

                    processSpliceSegment(firstSegment);
                    processSpliceSegment(secondSegment);

                    writeAlignmentSegment(writerHelper, firstSegment, 0);
                    // TODO: How to determine novel from non-novel splices
                    gobyAlEntry_setSplicedFlags(writerHelper, 1 /* normal splice, non-novel*/);
                    gobyAlEntry_setSplicedForwardFragmentIndex(writerHelper, 1);
                    gobyAlEntry_setSplicedForwardPosition(writerHelper, secondSegment->targetStart);
                    gobyAlEntry_setSplicedForwardTargetIndex(writerHelper, secondSegment->targetIndex);
                    // Different values for splicing. No need to adjust for indels
                    // as Gsnap doesn't support indels in spliced alignments
                    gobyAlEntry_setQueryAlignedLength(writerHelper,
                            firstSegment->queryEnd - firstSegment->queryStart + 1);
                    gobyAlEntry_setTargetAlignedLength(writerHelper,
                            firstSegment->queryEnd - firstSegment->queryStart + 1);

                    writeAlignmentSegment(writerHelper, secondSegment, 1);
                    gobyAlEntry_setSplicedFlags(writerHelper, 1 /* normal splice, non-novel*/);
                    gobyAlEntry_setSplicedBackwardFragmentIndex(writerHelper, 0);
                    gobyAlEntry_setSplicedBackwardPosition(writerHelper, firstSegment->targetStart);
                    gobyAlEntry_setSplicedBackwardTargetIndex(writerHelper, firstSegment->targetIndex);
                    // Different values for splicing. No need to adjust for indels
                    // as Gsnap doesn't support indels in spliced alignments
                    gobyAlEntry_setQueryAlignedLength(writerHelper,
                            secondSegment->queryEnd - secondSegment->queryStart + 1);
                    gobyAlEntry_setTargetAlignedLength(writerHelper,
                            secondSegment->queryEnd - secondSegment->queryStart + 1);
                    // Next loop interation
                    continue;
                }
            }

            // Non spliced and half splices
            GsnapAlignmentSegment *mergedSeg = mergeSegments(
                    alignmentEntry->alignmentSegments, splicedAlignment);

            writeAlignmentSegment(writerHelper, mergedSeg, 0);
            
            if (splicedAlignment) {
                // half splice...
                // Different values for splicing. No need to adjust for indels
                // as Gsnap doesn't support indels in spliced alignments
                gobyAlEntry_setQueryAlignedLength(writerHelper,
                        mergedSeg->queryEnd - mergedSeg->queryStart + 1);
                gobyAlEntry_setTargetAlignedLength(writerHelper,
                        mergedSeg->queryEnd - mergedSeg->queryStart + 1);
            }
        }
    }

    /**
     * Write the Gsnap alignment from the Gsnap data structures into Goby
     * compact-alignment format.
     * @param writerHelper the Goby writer helper
     * @param alignmentStr the lines that make up a single Gsnap alignment
     */
    void writeGobyAlignment(CAlignmentsWriterHelper *writerHelper) {
        GsnapAlignment *alignment = writerHelper->gsnapAlignment;
        if (alignment->pairedEnd) {
            // Pair alignments not supported just yet
            // writePairAlignment(writerHelper);
        } else {
            writeSingleAlignment(writerHelper);
        }
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
    void gobyGsnap_parse(CAlignmentsWriterHelper *writerHelper, char *alignment) {
        parseToGsnapDataStructures(writerHelper, alignment);
        writeGobyAlignment(writerHelper);
    }
}
