/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;

/**
 * An an implementation of AlignmentWriter that buffers entries up to a certain capacity and ensures
 * entries in the buffer are sorted by targetIndex/position before they are written to disk.
 *
 * @author Fabien Campagne
 *         Date: 4/17/12
 *         Time: 5:31 PM
 */
public class BufferedSortingAlignmentWriter implements AlignmentWriter {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(BufferedSortingAlignmentWriter.class);

    private final AlignmentWriter delegate;
    private final ObjectHeapPriorityQueue<Alignments.AlignmentEntry> heap;
    private static final int DEFAULT_CAPACITY = 1000;
    private static final Comparator<? super Alignments.AlignmentEntry> GENOMIC_POSITION_COMPARATOR = new AlignmentPositionComparator();
    private int capacity;
    /**
     * The maximum target index seen so far.
     */
    private int frontTargetIndex;
    /**
     * The maximum position seen so far on the targetIndex.
     */
    private int frontPosition;

    public BufferedSortingAlignmentWriter(final AlignmentWriter destination, final int capacity) {
        this.capacity = capacity;
        this.delegate = destination;
        this.heap = new ObjectHeapPriorityQueue<Alignments.AlignmentEntry>(capacity, GENOMIC_POSITION_COMPARATOR);
    }

    public BufferedSortingAlignmentWriter(final AlignmentWriter destination) {
        this(destination, DEFAULT_CAPACITY);
    }

    @Override
    public void appendEntry(final Alignments.AlignmentEntry entry) throws IOException {

        while (heap.size() > capacity) {
            final Alignments.AlignmentEntry queueEntry = heap.dequeue();
            checkFront(queueEntry);
            delegate.appendEntry(queueEntry);
        }
        heap.enqueue(entry);
    }

    private void checkFront(Alignments.AlignmentEntry entry) {
        final int targetIndex = entry.getTargetIndex();
        final int position = entry.getPosition();
        if (targetIndex < frontTargetIndex || targetIndex == frontTargetIndex && position < frontPosition) {
            // we detected an entry that occurs before the front of dequeued entries. We failed to restore sort order
            // we mark the destination writer as non-sorted and inform the end user with a warning.
            delegate.setSorted(false);
            LOG.warn("Local sorting strategy failed to restore sort order. The destination has been marked as unsorted. You must sort the output manually to improve compression.");
        }
        if (frontTargetIndex != targetIndex) {
            frontPosition = 0;
        }
        frontTargetIndex = Math.max(frontTargetIndex, targetIndex);
        frontPosition = Math.max(frontPosition, position);
    }

    @Override
    public void close() throws IOException {
        while (!heap.isEmpty()) {
            final Alignments.AlignmentEntry queueEntry = heap.dequeue();
            checkFront(queueEntry);
            delegate.appendEntry(queueEntry);
        }
        delegate.close();
    }

    @Override
    public void printStats(PrintStream out) {
        delegate.printStats(out);
    }

    @Override
    public void putStatistic(String description, double value) {
        delegate.putStatistic(description, value);
    }

    @Override
    public void putStatistic(String description, int value) {
        delegate.putStatistic(description, value);
    }

    @Override
    public void putStatistic(String description, String value) {
        delegate.putStatistic(description, value);
    }

    @Override
    public void setAlignerName(String alignerName) {
        delegate.setAlignerName(alignerName);
    }

    @Override
    public void setAlignerVersion(String alignerVersion) {
        delegate.setAlignerVersion(alignerVersion);
    }

    @Override
    public void setNumQueries(int numQueries) {
        delegate.setNumQueries(numQueries);
    }

    @Override
    public void setNumTargets(int numTargets) {
        delegate.setNumTargets(numTargets);
    }

    @Override
    public void setPermutation(boolean state) {
        delegate.setPermutation(state);
    }

    @Override
    public void setQueryIdentifiers(IndexedIdentifier queryIdentifiers) {
        delegate.setQueryIdentifiers(queryIdentifiers);
    }

    @Override
    public void setQueryIdentifiersArray(String[] queryIdentifiersArray) {
        delegate.setQueryIdentifiersArray(queryIdentifiersArray);
    }

    @Override
    public void setReadOriginInfo(ObjectArrayList<Alignments.ReadOriginInfo.Builder> readOriginInfoBuilderList) {
        delegate.setReadOriginInfo(readOriginInfoBuilderList);
    }

    @Override
    public void setSorted(boolean sortedState) {
        delegate.setSorted(sortedState);
    }

    @Override
    public void setTargetIdentifiers(IndexedIdentifier targetIdentifiers) {
        delegate.setTargetIdentifiers(targetIdentifiers);
    }

    @Override
    public void setTargetIdentifiersArray(String[] targetIdentifiersArray) {
        delegate.setTargetIdentifiersArray(targetIdentifiersArray);
    }

    @Override
    public void setTargetLengths(int[] targetLengths) {
        delegate.setTargetLengths(targetLengths);
    }
}
