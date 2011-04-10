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

import it.unimi.dsi.fastutil.objects.ObjectList;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.Properties;

import edu.cornell.med.icb.identifier.IndexedIdentifier;

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
     * and position.
     *
     * @param targetIndex The index of the target sequence to skip to.
     * @param position    The position on the target sequence.
     * @return The next entry, at position or past position (if not entry at position is found).
     * @throws java.io.IOException If an error occurs reading the alignment header. The header is accessed to check that the alignment is sorted.
     */
    Alignments.AlignmentEntry skipTo(int targetIndex, int position) throws IOException;

    /**
     * Reposition the reader to a new target sequence and start position.
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
     * Returns whether this read has query length information.
     *
     * @return True or false.
     */
    boolean hasQueryLengths();

    /**
     * Returns a sample of locations covered by this alignment.
     * @param modulo Modulo to avoid sampling every position in the genome.
     * @return  A set of positions that do occur in the genome, rounded to the specified modulo value (absoluteLocation-(absoluteLocation % modulo)).
     *  * @throws IOException
     */
    ObjectList<ReferenceLocation> getLocations(int modulo) throws IOException;

    boolean isQueryLengthStoredInEntries();

    /**
     * Return the name of the aligner that produced this alignment.
    * @return the name of the aligner that produced this alignment.
    */
    String getAlignerName();

    /**
     * Return the version of the aligner that produced this alignment.
     * @return the version of the aligner that produced this alignment.
     */
    String getAlignerVersion();

    int getSmallestSplitQueryIndex();

    int getLargestSplitQueryIndex();

    IndexedIdentifier getTargetIdentifiers();

    int[] getTargetLength();

    int getNumberOfTargets();

    int getNumberOfQueries();

    int[] getQueryLengths();
}
