/*
 * C structs to support GSNAP alignment / parsing.
 */

#ifndef C_GSNAP_STRUCTS_H_
#define C_GSNAP_STRUCTS_H_

/**
 * Structures for gsnap parsing in C++.
 */
#ifdef __cplusplus
    #include <vector>
    struct GsnapAlignmentSegment {
        char *readSequence;
        unsigned int readStart;
        unsigned int readEnd;
        int reverseStrand;
        char *targetIdentifier;
        unsigned int targetIndex;
        char *startType;
        int startClip;
        char *endType;
        int endClip;
        int matches;
        int subs;
    };

    struct GsnapAlignmentEntry {
        std::vector<GsnapAlignmentSegment*> *alignmentSegments;
    };

    struct GsnapAlignment {
        int lineNum;
        int pairedEnd;
        char *pairType;
        char *referenceSequence;
        char *referenceQuality;
        unsigned queryIndex;
        int numAlignmentEntries;
        std::vector<GsnapAlignmentEntry*> *alignmentEntries;
        int numAlignmentEntriesPair;
        std::vector<GsnapAlignmentEntry*> *alignmentEntriesPair;
        // If there are no entries, the below won't be populated.
        int alignScore;
        int mapq;
        int pairScore;
        int insertLength;
        char *pairSubType;
    };
#endif

#endif /* C_GSNAP_STRUCTS_H_ */
