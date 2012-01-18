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
        int readStart;
        int readEnd;
        int reverseStrand;
        char *targetIdentifier;
        char *targetIndex;
        char *startType;
        int startClip;
        char *endType;
        int endClip;
        int matches;
        int subs;
        int alignScore;
        int mapq;
        int pairScore;
        int insertLength;
    };

    struct GsnapAlignmentEntry {
        std::vector<GsnapAlignmentSegment*> *alignmentSegments;
    };

    struct GsnapAlignment {
        int pairedEnd;
        char *pairType;
        char *referenceSequence;
        char *referenceQuality;
        char *queryIdentifer;
        unsigned queryIndex;
        std::vector<GsnapAlignmentEntry*> *alignmentEntries;
        std::vector<GsnapAlignmentEntry*> *alignmentEntriesPair;
    };
#endif

#endif /* C_GSNAP_STRUCTS_H_ */
