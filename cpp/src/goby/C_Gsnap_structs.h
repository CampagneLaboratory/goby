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

    enum SegmentStartEndType {
        STARTENDTYPE_START,
        STARTENDTYPE_END,
        STARTENDTYPE_INS,
        STARTENDTYPE_DEL,
        STARTENDTYPE_TERM,
        STARTENDTYPE_DONOR,
        STARTENDTYPE_ACCEPTOR,
        STARTENDTYPE_UNKNOWN
    };

    enum SegmentSpliceDir {
            //"sense", "antisense"
        SPLICEDIR_NONE, 
        SPLICEDIR_SENSE, 
        SPLICEDIR_ANTISENSE, 
        SPLICEDIR_UNKNOWN 
    };

    enum SegmentSpliceType {
           //"consistent", "inversion", "scramble", "translocation"
        SPLICETYPE_NONE,
        SPLICETYPE_CONSISTENT,
        SPLICETYPE_INVERSION,
        SPLICETYPE_SCRAMBLE,
        SPLICETYPE_TRANSLOCATION,
        SPLICETYPE_UNKNOWN
    };

    enum PairType {
           //"concordant", "paired", "unpaired"
        PAIRTYPE_NONE,
        PAIRTYPE_CONCORDANT,
        PAIRTYPE_PAIRED,
        PAIRTYPE_UNPAIRED,
        PAIRTYPE_UNKNOWN
    };

    enum PairFlag {
        // the read is paired in sequencing, no matter whether it is mapped in a pair 
        PAIRFLAG_PAIRED                 =   1, // 0b0000000001
        // the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)  
        PAIRFLAG_PROPERLY_PAIRED        =   2, // 0b0000000010
        // the query sequence itself is unmapped 
        PAIRFLAG_READ_UNMAPPED          =   4, // 0b0000000100
        // the mate is unmapped  
        PAIRFLAG_MATE_UNMAPPED          =   8, // 0b0000001000
        // strand of the query (0 for forward; 1 for reverse strand) 
        PAIRFLAG_READ_REVERSE_STRAND    =  16, // 0b0000010000
        // strand of the mate  
        PAIRFLAG_MATE_REVERSE_STRAND    =  32, // 0b0000100000
        // the read is the ﬁrst read in a pair  
        PAIRFLAG_FIRST_IN_PAIR          =  64, // 0b0001000000
        // the read is the second read in a pair 
        PAIRFLAG_SECOND_IN_PAIR         = 128, // 0b0010000000
        // the alignment is not primary (a read having split hits may have multiple primary alignment records) 
        PAIRFLAG_NOT_PRIMARY_ALIGNMENT  = 256  // 0b0100000000
    };

    struct GsnapAlignmentSegment {
        string *referenceSequence;   // Will be generated during merge
        string *querySequence;       // The query sequence
        string *queryQuality;        // If !reverseStrand this will be null, use value in alignment
        string *segmentSequence;     // If reverseStrand, this will be reverseComplemented to be in the direction of the reference
        string *deletesSequence;     // If this sequence has deletes, they will be stored here otherwise NULL.
        string *insertsSequence;     // If this sequence has deletes, they will be stored here otherwise NULL.
        int queryStart;              // 0-based. Will be adjusted if reverseStrand
        int queryEnd;                // 0-based. Will be adjusted if reverseStrand
        bool reverseStrand;
        char *targetIdentifier;
        unsigned int targetIndex;    // 0-based
        unsigned int targetStart;    // 0-based
        unsigned int targetEnd;      // 0-based
        SegmentStartEndType startType;
        int startClip;
        double startProb;             // Only if spliced
        SegmentStartEndType endType;
        int endClip;
        double endProb;               // Only if spliced
        SegmentSpliceDir spliceDir;   // Only if spliced
        SegmentSpliceType spliceType; // Only if spliced
        unsigned int spliceDistance;  // Only if spliced
        int matches;  // Comes intially from gsnap, but we re-calculate
        int subs;     // Comes intially from gsnap, but we re-calculate
        int inserts;  // Number of inserts (calculated, doesn't come from Gsnap)
        int deletes;   // Number of deletes (calculated, doesn't come from Gsnap)
        bool merged;
        int fragmentIndex;
    };

    struct GsnapAlignmentEntry {
        std::vector<GsnapAlignmentSegment*> *alignmentSegments;
    };

    struct GsnapAlignment {
        int lineNum;
        bool pairedEnd;
        PairType pairType;
        string *querySequence;    // ONLY for parsing, use the version in Segment during output
        string *queryQuality;     // ONLY for parsing, use the version in Segment during output
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
        string *pairSubType;
        int nextFragmentIndex;
    };
#endif

#endif /* C_GSNAP_STRUCTS_H_ */
