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
    private final AlignmentWriter delegate;
    private final ObjectHeapPriorityQueue<Alignments.AlignmentEntry> heap;
    private static final int DEFAULT_CAPACITY = 1000;
    private static final Comparator<? super Alignments.AlignmentEntry> GENOMIC_POSITION_COMPARATOR = new AlignmentPositionComparator();
    private int capacity;

    public BufferedSortingAlignmentWriter(final AlignmentWriter destination, final int capacity) {
        this.capacity = capacity;
        this.delegate=destination;
        this.heap = new ObjectHeapPriorityQueue<Alignments.AlignmentEntry>(capacity, GENOMIC_POSITION_COMPARATOR);
    }

    public BufferedSortingAlignmentWriter(final AlignmentWriter destination) {
        this(destination, DEFAULT_CAPACITY);
    }

    @Override
    public void appendEntry(final Alignments.AlignmentEntry entry) throws IOException {
        while (heap.size() > capacity) {
            delegate.appendEntry(heap.dequeue());
        }
        heap.enqueue(entry);
    }


    @Override
    public void close() throws IOException {
        while (!heap.isEmpty()) {
            delegate.appendEntry(heap.dequeue());
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
