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

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.Properties;

/**
 * @author Fabien Campagne
 *         Date: Apr 10, 2011
 *         Time: 12:43:51 PM
 */
public interface AlignmentReader extends Closeable, Iterator<Alignments.AlignmentEntry>, Iterable<Alignments.AlignmentEntry> {
    /**
     * Returns whether this alignment is sorted. Entries in a sorted alignment appear in order of
     * increasing target index and position.
     *
     * @return True if this position is sorted by position. False otherwise.
     */
    boolean isSorted();

    /**
     * Returns true if this alignment is indexed. When this method returns true,
     * a file called 'basename'.index is expected.
     *
     * @return
     */
    boolean isIndexed();

    /**
     * Returns the basename for the alignment being read, or null if the basename is unknown.
     * If a path was part of basename provided to the constructor, it is returned.
     *
     * @return basename for the alignment being read
     */
    String basename();

    /**
     * Returns true if the input has more entries.
     *
     * @return true if the input has more entries, false otherwise.
     */
    boolean hasNext();

    /**
     * Returns the next alignment entry from the input stream.
     *
     * @return the alignment read entry from the input stream.
     */
    Alignments.AlignmentEntry next();

    /**
     * Skip all entries that have position before (targetIndex,position). This method will use the alignment index
     * when it is available to skip directly to the closest chunk start before the entry identified by targetIndex
     * and position. Please note that this method skips entries, it does not go back to previously visited entries.
     * This is useful to simply iterate through an alignment visiting only certain regions.
     * If you need full random access in an alignment, you need to call reposition(ref,pos) followed by skipTo(ref,pos).
     *
     * @param targetIndex The index of the target sequence to skip to.
     * @param position    The position on the target sequence.
     * @return The next entry, at position or past position (if not entry at position is found).
     * @throws java.io.IOException If an error occurs reading the alignment header. The header is accessed to check that the alignment is sorted.
     */
    Alignments.AlignmentEntry skipTo(int targetIndex, int position) throws IOException;

    /**
     * Reposition the reader to a new target sequence and start position. Since 1.9.6 this method supports
     * full random access and will reposition to locations earlier or after the current location. Goby 1.9.5-
     * would only reposition to locations after the current location. Please be advised that calling reposition
     * relocates to the start of the chunk that contains the entry with targetIndex and position, it does not
     * garantee that a subsequent call to next() will return an entry at location targetIndex, position. You must
     * call skipTo after reposition if you require this behavior.
     *
     * @param targetIndex Index of the target sequence to reposition to.
     * @param position    Position in the target sequence to reposition to.
     * @throws java.io.IOException If an error occurs repositioning.
     */
    void reposition(int targetIndex, int position) throws IOException;

    /**
     * This operation is not supported.
     */
    void remove();

    /**
     * Read the header of this alignment.
     *
     * @throws java.io.IOException If an error occurs.
     */
    void readHeader() throws IOException;

    /**
     * Read the index. The header is also loaded.
     *
     * @throws java.io.IOException If an error occurs loading the index or header.
     */
    void readIndex() throws IOException;

    /**
     * {@inheritDoc}
     */
    void close();

    Iterator<Alignments.AlignmentEntry> iterator();

    Properties getStatistics();

    int getNumberOfAlignedReads();


    /**
     * Returns a sample of locations covered by this alignment.
     *
     * @param modulo Modulo to avoid sampling every position in the genome.
     * @return A set of positions that do occur in the genome, rounded to the specified modulo value (absoluteLocation-(absoluteLocation % modulo)).
     *         * @throws IOException
     */
    ObjectList<ReferenceLocation> getLocations(int modulo) throws IOException;

    boolean isQueryLengthStoredInEntries();

    /**
     * Return the name of the aligner that produced this alignment.
     *
     * @return the name of the aligner that produced this alignment.
     */
    String getAlignerName();

    /**
     * Return the version of the aligner that produced this alignment.
     *
     * @return the version of the aligner that produced this alignment.
     */
    String getAlignerVersion();

    int getSmallestSplitQueryIndex();

    int getLargestSplitQueryIndex();

    IndexedIdentifier getTargetIdentifiers();

    int[] getTargetLength();

    int getNumberOfTargets();

    int getNumberOfQueries();


    boolean isConstantQueryLengths();

    int getConstantQueryLength();


    IndexedIdentifier getQueryIdentifiers();

    /**
     * Obtain the byte offset at the start of the chunk with entries at or after the given genomic position.
     *
     * @param startReferenceIndex Index of the reference sequence for the genomic position.
     * @param startPosition       genomic coordinate in the reference sequence.
     * @return byte offset in the entries file.
     */
    long getStartByteOffset(int startReferenceIndex, int startPosition);

    /**
     * Indicate whether query indices are small indices. When true, create a PermutationReader to reconstitute
     * original query indices from the stored small indices.
     *
     * @return True when a permutation was created, false otherwise.
     */
    boolean getQueryIndicesWerePermuted();

    /**
     * Obtain the byte offset at the end of the chunk with entries at or before the given genomic end position.
     * The start position is taken as an argument to ensure that the end offset contains at least one chunk.
     * You must provide the same start arguments that you use to call getStartByteOffset. When called in this
     * way,getStartByteOffset and getEndyteOffset provide a slice of the alignment that contain entries whose
     * positions p are startPosition <= p < endPosition
     *
     * @param startReferenceIndex Index of the reference sequence for the start genomic position.
     * @param startPosition       start genomic coordinate in the reference sequence.
     * @param endReferenceIndex   Index of the reference sequence for the end genomic position.
     * @param endPosition         end genomic coordinate in the reference sequence.
     * @return byte offset in the entries file.
     */
    long getEndByteOffset(int startReferenceIndex, int startPosition, int endReferenceIndex, int endPosition);
}
