/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.alignments;

import java.io.IOException;

/**
 * A factory that returns alignment reader. This interface can be subclassed to provide specific implementations of
 * the AlignmentReader interface. This is useful to provide implementations that perform some filtering on the fly.
 *
 * @author Fabien Campagne
 *         Date: Apr 9, 2011
 *         Time: 4:29:29 PM
 */
public interface AlignmentReaderFactory {

    /**
     * Returns an implementation of AlignmentReader.
     *
     * @param basename Basename of the alignment.
     * @return implementation of alignment reader, typically customized for some task.
     * @throws java.io.IOException If an error occurs opening a reader.
     */
    public AlignmentReader createReader(String basename) throws IOException;

    /**
     * Returns an array to hold AlignmentReader implementations.
     *
     * @param numElements Number of elements in the array.
     * @return an array of the alignment reader type.
     * @throws java.io.IOException If an error occurs opening a reader.
     */
    AlignmentReader[] createReaderArray(int numElements) throws IOException;

    /**
     * Create a reader for a specific slice of an alignment file contained exactly between a start
     * and an end location. Start and end locations are genomic/reference positions. Entries will be returned
     * that occur from the start position and up to the end position (start and end positions are inclusive).
     *
     * @param basename            Basename for the alignemnt.
     * @param startReferenceIndex Index of the reference for the start position.
     * @param startPosition       Position on the reference for the start position.
     * @param endReferenceIndex   Index of the reference for the end position.
     * @param endPosition         Position on the reference for the end position.
     * @return an alignment reader configured over the genomic slice/range.
     * @throws IOException Thrown if an error occurs opening or reading the alignment file.
     */
    AlignmentReader createReader(String basename,
                                 int startReferenceIndex, int startPosition,
                                 int endReferenceIndex, int endPosition) throws IOException;

    /**
     * Create a reader for a specific slice of an alignment file contained exactly between a start
     * and an end location. Start and end locations are genomic/reference positions. Entries will be returned
     * that occur from the start position and up to the end position in range (start and end positions are inclusive).
     * If range is null, the method defaults to createReader(basename) and opens the entire alignment.
     * @param basename Basename for the alignemnt.
     * @param range    Range/slice of the genome that the reader will be restricted to.
     * @return an alignment reader configured over the genomic range.
     * @throws IOException Thrown if an error occurs opening or reading the alignment file.
     */
    AlignmentReader createReader(String basename,
                                 GenomicRange range) throws IOException;

    /**
     * Create a reder for reading between the byte positions startOffset and endOffset.
     *
     * @param basename    Basename of the alignment to read.
     * @param startOffset Position in the file where reading will start (in bytes).
     * @param endOffset   Position in the file where reading will end (in bytes).
     * @return an alignment reader constrained to the offsets.
     * @throws IOException If an error occurs opening or reading the file.
     */
    AlignmentReader createReader(String basename, long startOffset, long endOffset) throws IOException;

    /**
     * Obtain a file slice for the specified genomic range, in the specified alignment. When range is null, returns
     * a slice corresponding to the complete file.
     * @param basename Basename of the alignment for the slice refers to.
     * @param range    Genomic range for which a slice should be generated.
     * @return a suitable FileSlice.
     */
    FileSlice getSlice(String basename, GenomicRange range) throws IOException;
}
