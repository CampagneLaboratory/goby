/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.Properties;

/**
 * Read over a set of alignments. This aligner concatenates entries from the input alignment.
 * Reference sequences must match exactly across the input alignments.
 * Query are assumed to be entirely distinct and will be treated as independent observations (e.g.,
 * reads from multiple independent samples). To this effect, alignment entries read from
 * different input basenames, which would otherwise share an identical query index,
 * are renumbered with distinct query indices.
 *
 * @author Fabien Campagne
 *         Date: May 20, 2009
 *         Time: 5:06:01 PM
 */
public class ConcatAlignmentReader extends AbstractAlignmentReader
        implements Iterator<Alignments.AlignmentEntry>, Iterable<Alignments.AlignmentEntry>, Closeable {

    AlignmentReader[] readers;
    IntSet readersWithMoreEntries;
    int[] numQueriesPerReader;
    int[] queryIndexOffset;

    private int activeIndex;
    private boolean adjustQueryIndices = true;
    private int numberOfAlignedReads;

    /**
     * Construct an alignment reader over a set of alignments.
     * Please note that the constructor access the header of each individual alignment to
     * check reference sequence identity and obtain the number of queries in each input alignment.
     *
     * @param basenames Basenames of the individual alignemnts to combine.
     * @throws IOException If an error occurs reading the header of the alignments.
     */
    public ConcatAlignmentReader(final String... basenames) throws IOException {
        readers = new AlignmentReader[basenames.length];
        readersWithMoreEntries = new IntArraySet();
        int readerIndex = 0;
        for (final String basename : basenames) {
            readers[readerIndex] = new AlignmentReader(basename);
            readersWithMoreEntries.add(readerIndex);
            readerIndex++;
        }
        numQueriesPerReader = new int[basenames.length];
        queryIndexOffset = new int[basenames.length];
        readHeader();
    }

    /**
     * Read the header of this alignment.
     *
     * @throws java.io.IOException If an error occurs.
     */
    @Override
    public final void readHeader() throws IOException {
        if (!isHeaderLoaded()) {
            final IntSet targetNumbers = new IntArraySet();
            int readerIndex = 0;
            numberOfQueries = 0;

            for (final AlignmentReader reader : readers) {
                reader.readHeader();

                targetNumbers.add(reader.getNumberOfTargets());
                final int numQueriesForReader = reader.getNumberOfQueries();
                numQueriesPerReader[readerIndex] = numQueriesForReader;
                if (adjustQueryIndices) {
                    numberOfQueries += numQueriesForReader;
                } else {
                    numberOfQueries = Math.max(numberOfQueries, numQueriesForReader);
                }
                readerIndex++;
                numberOfAlignedReads += reader.getNumberOfAlignedReads();
            }
            if (targetNumbers.size() != 1) {
                throw new IllegalArgumentException("The number of targets must match exactly across the input basenames. Found " + targetNumbers.toString());
            } else {
                this.numberOfTargets = targetNumbers.iterator().nextInt();
            }
            // they are all the same, use the first one:
            targetIdentifiers = readers[0].getTargetIdentifiers();

            queryLengths = new int[numberOfQueries];

            for (int i = 1; i < queryIndexOffset.length; i++) {
                final int offset = numQueriesPerReader[i - 1] + queryIndexOffset[i - 1];
                queryIndexOffset[i] = offset;
                final int[] localQueryLenth = readers[i].getQueryLength();
                final int j = 0;
                if (localQueryLenth != null) {
                    for (final int length : localQueryLenth) {

                        queryLengths[j + offset] = length;
                    }
                }
            }

            setHeaderLoaded(true);
        }
    }

    @Override
    public final Alignments.AlignmentEntry nextAlignmentEntry() {
        if (!hasNext()) {
            throw new IllegalStateException("next() cannot be called when hasNext() would have returned false.");
        } else {
            final Alignments.AlignmentEntry alignmentEntry = readers[activeIndex].next();
            final int newQueryIndex = mergedQueryIndex(alignmentEntry.getQueryIndex());
            if (adjustQueryIndices) {
                return alignmentEntry.newBuilderForType().mergeFrom(alignmentEntry).setQueryIndex(newQueryIndex).build();
            } else {
                return alignmentEntry;
            }
        }
    }

    private int mergedQueryIndex(final int queryIndex) {
        return queryIndexOffset[activeIndex] + queryIndex;
    }

    @Override
    public final boolean hasNextAligmentEntry() {

        while (!readersWithMoreEntries.isEmpty()) {
            activeIndex = readersWithMoreEntries.iterator().nextInt();
            final AlignmentReader reader = readers[activeIndex];
            final boolean hasNext = reader.hasNext();
            if (!hasNext) {
                readersWithMoreEntries.remove(activeIndex);
            } else {
                return true;
            }

        }
        return false;
    }

    /**
     * Iterator over alignment entries.
     */
    public final Iterator<Alignments.AlignmentEntry> iterator() {
        return this;
    }

    public final boolean hasNext() {
        return hasNextAligmentEntry();
    }


    public final Alignments.AlignmentEntry next() {
        return nextAlignmentEntry();
    }


    /**
     * This operation is not supported.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a reader.");
    }


    public void setAdjustQueryIndices(final boolean adjustQueryIndices) {
        this.adjustQueryIndices = adjustQueryIndices;
    }

    /**
     * Obtain statistics about this alignment as a Java property instance.
     *
     * @return
     */
    public Properties getStatistics() {
        int index = 1;
        Properties result = new Properties();
        for (AlignmentReader reader : this.readers) {
            Properties localProps = reader.getStatistics();

            for (Object key : localProps.keySet()) {
                result.put("part" + index + "." + key.toString(), localProps.get(key));
            }
            index++;
        }
        return result;
    }

    public int getNumberOfAlignedReads() {
        return numberOfAlignedReads;
    }

    /**
     * Close the underlying readers.
     * @throws IOException
     */
    public void close() throws IOException {
        for (AlignmentReader reader : readers) {
            reader.close();
        }
    }
}
