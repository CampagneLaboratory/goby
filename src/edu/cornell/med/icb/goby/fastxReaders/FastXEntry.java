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

package edu.cornell.med.icb.goby.fastxReaders;

import it.unimi.dsi.lang.MutableString;

/**
 * A FASTA / FASTQ / similar entry.
 * TODO: We must accept ANY character in the quality after the first quality
 * TODO: character until the size matches the sequence size.
 * @author Kevin Dorff
 */
public class FastXEntry {
    /** The header symbol. */
    private char headerSymbol;

    /** The entry. */
    private MutableString entry;

    /** The entry without the header. */
    private MutableString entrySansHeader;

    /** The sequence header line WITH symbol. */
    private MutableString sequenceHeader;

    /** The sequence string. */
    private MutableString sequence;

    /** The quality header line (if FASTQ), WITH symbol. */
    private MutableString qualityHeader;

    /** The quality String. */
    private MutableString quality;

    /** True if the entry is complete. */
    private boolean entryComplete;

    /** Constructor. */
    public FastXEntry() {
        entry = new MutableString();
        entrySansHeader = new MutableString();
        sequenceHeader = new MutableString();
        sequence = new MutableString();
        qualityHeader = new MutableString();
        quality = new MutableString();
        reset();
    }

    /**
     * Reset this FASTA entry so it can be reused.
     */
    public void reset() {
        entry.length(0);
        entrySansHeader.length(0);
        sequenceHeader.length(0);
        sequence.length(0);
        qualityHeader.length(0);
        quality.length(0);
        entryComplete = false;
        headerSymbol = '\0';
    }

    /**
     * Add a line to the FASTX entry. Ignores empty lines.
     * @param line the line to add
     * @return if true is returned, the line was added file. If false is returned
     * the line is into the next record and was not used
     */
    public boolean addLine(final String line) {
        if (line.length() == 0) {
            return true;
        }
        final char symbol = line.charAt(0);
        boolean setSequenceHeader = false;
        if (sequenceHeader.length() == 0) {
            // First line, this is the sequence header
            headerSymbol = symbol;
            // header does not include the header symbol:
            sequenceHeader.append(line.subSequence(1, line.length()));
            setSequenceHeader = true;
        } else {
            if ((symbol == '@' || symbol == '>')
                    && (sequenceHeader.length() > 0)
                    && (qualityHeader.length() == 0)) {
                entryComplete = true;
                // Next FASTA entry (not FASTQ)
                return false;
            }
            if (symbol == '+' && qualityHeader.length() == 0) {
                // Quality header line
                qualityHeader.append(line);
            } else {
                if (qualityHeader.length() > 0) {
                    quality.append(line);
                } else {
                    sequence.append(line);
                }
            }
        }
        if ((quality.length() > 0) && (quality.length() >= sequence.length())) {
            // Next FASTQ
            entryComplete = true;
        }
        if (!setSequenceHeader) {
            if (entrySansHeader.length() > 0) {
                entrySansHeader.append('\n');
            }
            entrySansHeader.append(line);
        }
        if (entry.length() > 0) {
            entry.append('\n');
        }
        entry.append(line);
        return true;
    }

    /**
     * Get the FASTX entry.
     * @return the FASTX entry
     */
    public char getHeaderSymbol() {
        return headerSymbol;
    }

    /**
     * Get the FASTX entry.
     * @return the FASTX entry
     */
    public MutableString getEntry() {
        return entry;
    }

    /**
     * Get the FASTX sequence.
     * @return the FASTX entry sequence
     */
    public MutableString getSequence() {
        return sequence;
    }

    /**
     * Get the FASTQ quality.
     * @return the FASTX entry quality
     */
    public MutableString getQuality() {
        return quality;
    }

    /**
     * Get the FASTX entry without the header.
     * @return the FASTX entry without the header
     */
    public MutableString getEntrySansHeader() {
        return entrySansHeader;
    }

    /**
     * Get the FASTX entry header (first line of the entry not including
     * the header symbol (such as '>' or '@')).
     * @return the FASTX entry
     */
    public MutableString getEntryHeader() {
        return sequenceHeader;
    }

    /**
     * Get the read length for the FASTX entry.
     * @return the read length
     */
    public int getReadLength() {
        return sequence.length();
    }

    public boolean isEntryComplete() {
        return entryComplete;
    }

    /**
     * Copy the current object to a new one - you can use this if you need to store
     * this object in a list, etc.
     * @return this object copied to a new object.
     */
    public FastXEntry copy() {
        final FastXEntry copy = new FastXEntry();
        /** The entry. */
        copy.headerSymbol = headerSymbol;
        copy.entry = new MutableString(entry);
        copy.sequenceHeader = new MutableString(sequenceHeader);
        copy.sequence = new MutableString(sequence);
        copy.qualityHeader = new MutableString(qualityHeader);
        copy.quality = new MutableString(quality);
        copy.entryComplete = entryComplete;
        return copy;
    }
}
