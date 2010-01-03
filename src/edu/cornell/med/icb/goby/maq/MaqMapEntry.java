/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.maq;

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.lang.MutableString;

/**
 * Class to store a MAQ MAP entry.
 * @author Kevin Dorff
 */
public class MaqMapEntry {
    /** The max read length. */
    private final int maxReadLen;

    /** The sequence. */
    private short[] seq;          //uint8[MAX_READLEN]
    /** The size. */
    private short size;           //uint8
    /** The mapQual. */
    private short mapQual;       //uint8
    /** The info1. */
    private short info1;          //uint8
    /** The info2. */
    private short info2;          //uint8
    /** The two c values. */
    private short[] c;            //uint8[2]
    /** The flag. */
    private short flag;           //uint8
    /** The alt_qual. */
    private short altQual;       //uint8
    /** The seqid. */
    private long seqId;           //uint32
    /** The pos. */
    private long pos;             //uint32
    /** The dist. */
    private int dist;             //int
    /** The name. */
    private MutableString name;          //char[MAX_NAMELEN]
    /** The name index. */
    private int readIndex;        //(not stored)

    /** This is a value calculated from info1. */
    private short numMisMatches;

    /** The indexed identifier. */
    private IndexedIdentifier readNameIdentifiers;

    /**
     * Construct an empty MaqMapEntry. Will be populated by MaqDataInputHelper.
     * If you use this version, getNameIndex() will always return a -1.
     * @param maxReadLenVal the max reads length
     */
    public MaqMapEntry(final int maxReadLenVal) {
        this(maxReadLenVal, null);
    }

    /**
     * Construct an empty MaqMapEntry. Will be populated by MaqDataInputHelper.
     * If readNameIdentifiersVal == null, getNameIndex() will always return a -1
     * otherwise it will return the index value as registered in readNameIdentifiersVal
     * for the term.
     * @param maxReadLenVal the max read length
     * @param readNameIdentifiersVal the IndexedIdentifier to use to store the "name" value
     */
    public MaqMapEntry(
            final int maxReadLenVal, final IndexedIdentifier readNameIdentifiersVal) {
        this.seq = new short[0];
        this.maxReadLen = maxReadLenVal;
        this.c = new short[2];
        this.readNameIdentifiers = readNameIdentifiersVal;
    }

    /**
     * Get seq.
     * @return seq
     */
    public short[] getSeq() {
        return seq;
    }

    /**
     * Set seq.
     * @param seq seq
     */
    public void setSeq(final short[] seq) {
        this.seq = seq;
    }

    /**
     * Get size.
     * @return size
     */
    public short getSize() {
        return size;
    }

    /**
     * Set size.
     * @param size size
     */
    public void setSize(final short size) {
        this.size = size;
    }

    /**
     * Get mapQual.
     * @return mapQual
     */
    public short getMapQual() {
        return mapQual;
    }

    /**
     * Set mapQual.
     * @param mapQual mapQual
     */
    public void setMapQual(final short mapQual) {
        this.mapQual = mapQual;
    }

    /**
     * Get the number of mis-matches. This is a value calculated from info1 and set when you
     * set info1.
     * This value is not written to map files.
     * @return the number of mis matches
     */
    public short getNumMisMatches() {
        return numMisMatches;
    }

    /**
     * Get info1.
     * This SEEMS to be (based on the Eland converter) scored related
     * to the number of misses. It appears that if
     * numberOfMismatches = X
     *     info1 = ((numberOfMismatches << 4 | numberOfMismatches) & 0xf)
     * and you can back calculate with
     *     numberOfMismatches = info1 & 0xf
     * @return info1
     */
    public short getInfo1() {
        return info1;
    }

    /**
     * Set info1.
     * @param info1 info1
     */
    public void setInfo1(final short info1) {
        this.info1 = info1;
        this.numMisMatches = (short) (this.info1 & 0xf);
    }

    /**
     * Get info2.
     * @return info2
     */
    public short getInfo2() {
        return info2;
    }

    /**
     * Set info2.
     * @param info2 info2
     */
    public void setInfo2(final short info2) {
        this.info2 = info2;
    }

    /**
     * Get C.
     * @return C
     */
    public short[] getC() {
        return c;
    }

    /**
     * Set c.
     * @param c c
     */
    public void setC(final short[] c) {
        this.c = c;
    }

    /**
     * Get flag.
     * @return flag
     */
    public short getFlag() {
        return flag;
    }

    /**
     * Set flag.
     * @param flag flag
     */
    public void setFlag(final short flag) {
        this.flag = flag;
    }

    /**
     * Get altQual.
     * @return altQual
     */
    public short getAltQual() {
        return altQual;
    }

    /**
     * Set altQual.
     * @param altQual altQual
     */
    public void setAltQual(final short altQual) {
        this.altQual = altQual;
    }

    /**
     * Get seqid.
     * @return seqid
     */
    public long getReferenceSequenceId() {
        return seqId;
    }

    /**
     * Set seqid.
     * @param seqId seqid
     */
    public void setSeqId(final long seqId) {
        this.seqId = seqId;
    }

    /**
     * Get pos as stored in MAQ file >>encoded<<. If you want the ACTUAL
     * position you need to use getActualPosition(). Do not use directly
     * unless you know what you are doing.
     * @return pos shifted and contains forward / reverse information.
     * Do not use directly!
     */
    public long getPos() {
        return pos;
    }

    /**
     * Get the actual position. getPos() contains information about
     * forward / backward, etc. and is shifted. USE THIS VERSION to know the
     * actual position.
     * @return the actual position.
     */
    public long getActualPosition() {
        return (this.pos >> 1) + 1;
    }

    /**
     * Get if the sequence is matching the reverse strand, calculated using this.pos.
     * @return if reversed
     */
    public boolean isMatchingReverseStrand() {
        return ((this.pos & 1) > 0);
    }

    /**
     * Set pos.
     * @param pos pos
     */
    public void setPos(final long pos) {
        this.pos = pos;
    }

    /**
     * Get dist.
     * @return dist
     */
    public int getDist() {
        return dist;
    }

    /**
     * Set dist.
     * @param dist dist
     */
    public void setDist(final int dist) {
        this.dist = dist;
    }

    /**
     * Get name.
     * @return name
     */
    public MutableString getName() {
        return name;
    }

    /**
     * Get name index (the index in indexedIdentifier to the name). Not stored.
     * @return name index
     */
    public int getReadNameIndex() {
        return readIndex;
    }

    /**
     * Set name.
     * @param name name
     */
    public void setReadName(final MutableString name) {
        this.name = name;
        if (readNameIdentifiers != null) {
            this.readIndex = readNameIdentifiers.registerIdentifier(name);
        } else {
            this.readIndex = -1;
        }
    }

     /**
     * Set name.
     * @param name name
     */
    public void setReadName(final String name) {
       setReadName(new MutableString(name));
    }
    /**
     * User readable output HEADER for this entry.
     * @return User readable output HEADER for this entry
     */
    public static String formatHeader() {
        final StringBuilder sb = new StringBuilder();
        sb.append(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
            "name",
            "ref-name",
            "pos->>1+1",
            "plus-minus",
            "dist",
            "flag",
            "map-qual",
            "seq-last-char",
            "alt-qual",
            "info1-and-0xf-num-mismatches?",
            "info2",
            "c-0",
            "c-1",
            "size",
            "sequence",
            "quality"));
        return sb.toString();
    }

    /**
     * User readable format of this entry.
     * @param header the header associated with this entry
     * @return User readable format of this entry
     */
    public String format(final MaqMapHeader header) {
        final MaqMapEntry entry = this;
        final StringBuilder sb = new StringBuilder();
        final char plusMinus;
        if ((this.pos & 1) > 0) {
            plusMinus = '-';
        } else {
            plusMinus = '+';
        }
        sb.append(String.format("%s\t%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                this.getName(),
                header.getRefName((int) this.seqId),
                (this.pos >> 1) + 1,
                plusMinus,
                this.dist,
                this.flag,
                this.mapQual,
                this.seq[maxReadLen - 1],
                this.altQual,
                this.info1 & 0xf,
                this.info2,
                this.c[0],
                this.c[1],
                this.size));

        sb.append('\t');
        for (int j = 0; j < entry.size; j++) {
            if (getSeqEntry(j) == 0) {
                sb.append('n');
            } else {
                sb.append("ACGT".charAt(entry.seq[j] >> 6 & 3));
            }
        }

        sb.append('\t');
        for (int j = 0; j < entry.size; j++) {
            sb.append((char) ((getSeqEntry(j) & 0x3f) + 33));
        }

        return sb.toString();
    }

    /**
     * Get the sequence entry at position i, but return 0 if seq is too small.
     * @param i the index to retrieve
     * @return the value at index i or 0
     */
    private short getSeqEntry(final int i) {
        if (i < seq.length) {
            return seq[i];
        } else {
            return 0;
        }
    }
}
