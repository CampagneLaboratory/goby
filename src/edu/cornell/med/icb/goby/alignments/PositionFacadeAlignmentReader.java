/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.alignments;

import java.io.*;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * A facade over alignment reader that restricts entries read to a contiguous window of positions.
 *
 * @author Fabien Campagne
 *         Date: Dec 13, 2010
 *         Time: 10:52:30 AM
 */
public class PositionFacadeAlignmentReader implements Closeable,
        Iterator<Alignments.AlignmentEntry>,
        Iterable<Alignments.AlignmentEntry> {
    
    private int endPosition;
    private int endReferenceIndex;
    private int startPosition;
    private int startReferenceIndex;
    private AlignmentReader delegate;


    /**
     * A constructor that allows reading a slice of an alignment file contained exactly between a start
     * and an end location. Start and end locations are genomic/reference positions. Entries will be returned
     * that occur after the start position and before the end position.
     *
     * @param basename            Basename for the alignemnt.
     * @param startReferenceIndex Index of the reference for the start position.
     * @param startPosition       Position on the reference for the start position.
     * @param endReferenceIndex   Index of the reference for the end position.
     * @param endPosition         Position on the reference for the end position.
     * @throws IOException Thrown if an error occurs opening or reading the alignment file.
     */
    public PositionFacadeAlignmentReader(final String basename,
                                         final int startReferenceIndex,
                                         final int startPosition,
                                         final int endReferenceIndex,
                                         final int endPosition)
            throws IOException {

        AlignmentReaderImpl reader = null;

        try {
            reader = new AlignmentReaderImpl(basename);


            reader.readHeader();
            if (!reader.isIndexed())
                throw new UnsupportedOperationException("The alignment must be sorted and indexed to read slices of data by reference position.");
            reader.readIndex();

            final long startOffset = reader.getByteOffset(startReferenceIndex, startPosition, 0);
            final long endOffset = reader.getByteOffset(endReferenceIndex, endPosition + 1, 1);
            this.endPosition = endPosition;
            this.endReferenceIndex = endReferenceIndex;
            this.startPosition = startPosition;
            this.startReferenceIndex = startReferenceIndex;
            this.delegate = new AlignmentReaderImpl(startOffset, endOffset, basename);

        } finally {
            if (reader != null) {
                reader.close();
            }
        }

    }

    public void close() throws IOException {
        delegate.close();
    }

    Alignments.AlignmentEntry nextEntry = null;

    /**
     * Determines if the underlying reader has an entry within the position window boundaries.
     *
     * @return True or False.
     */
    public boolean hasNext() {
        if (nextEntry != null) return true;

        int entryTargetIndex;
        do {

            if (!delegate.hasNext()) return false;

            nextEntry = delegate.next();
            entryTargetIndex = nextEntry.getTargetIndex();
        } while (entryTargetIndex < startReferenceIndex ||
                (entryTargetIndex == startReferenceIndex && nextEntry.getPosition() < startPosition));

        // Early stop if we are past the end location:
        if (entryTargetIndex > endReferenceIndex ||
                (entryTargetIndex == endReferenceIndex && nextEntry.getPosition() > endPosition)) {
            nextEntry = null;
            return false;
        }
        return true;
    }

    public Alignments.AlignmentEntry next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        try {
            return nextEntry;

        } finally {

            nextEntry = null;
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Removing elements is not supported by AlignmentReader implementations.");
    }

    public Iterator<Alignments.AlignmentEntry> iterator() {
        return this;
    }
}
