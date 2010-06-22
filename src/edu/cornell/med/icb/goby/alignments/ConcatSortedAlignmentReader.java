/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectArrayPriorityQueue;
import it.unimi.dsi.fastutil.AbstractPriorityQueue;

import java.util.NoSuchElementException;
import java.util.ArrayList;
import java.util.Comparator;
import java.io.IOException;

/**
 * Concatenates sorted alignments while preserving entry sort order across inputs. The result is a sorted
 * alignment.
 *
 * @author Fabien Campagne
 *         Date: Jun 22, 2010
 *         Time: 10:42:34 AM
 */
public class ConcatSortedAlignmentReader extends ConcatAlignmentReader {

    AbstractPriorityQueue<Bucket> entryHeap;
    private static AlignmentPositionComparator comparator = new AlignmentPositionComparator();
    private boolean[] nextLoadedForReader;
    private Bucket[] buckets;

    public ConcatSortedAlignmentReader(String... basenames) throws IOException {
        super(basenames);
        init(basenames);
    }

    private void init(String... basenames) {
        nextLoadedForReader = new boolean[basenames.length];


        final Comparator<Bucket> bucketComparator = new Comparator<Bucket>() {
            public int compare(Bucket bucket, Bucket bucket1) {
                return comparator.compare(bucket.entry, bucket1.entry);
            }
        };
        // For small number of input alignments, use the array priority queue.
        // For more than 10 input alignments, use the heap implementation.
        // This heuristic tries to pick the data structure implementation that will yield the best
        // performance.
        entryHeap = basenames.length < 10 ? new ObjectArrayPriorityQueue<Bucket>(bucketComparator) :
                new ObjectHeapPriorityQueue<Bucket>(bucketComparator);
        buckets = new Bucket[basenames.length];
        for (int i = 0; i < buckets.length; i++) {
            buckets[i] = new Bucket(null, i);
        }
    }

    public ConcatSortedAlignmentReader(boolean adjustQueryIndices, String... basenames) throws IOException {
        super(adjustQueryIndices, basenames);
        init(basenames);
    }

    class Bucket {
        Alignments.AlignmentEntry entry;
        int readerIndex;

        public Bucket(Alignments.AlignmentEntry alignmentEntry, int index) {
            this.entry = alignmentEntry;
            this.readerIndex = index;
        }
    }

    /**
     * The minimum entry, according to position sort order.
     */
    Alignments.AlignmentEntry minEntry = null;
    /**
     * The index of the reader that provided minEntry.
     */
    int minReaderIndex;

    /**
     * Returns true if the input has more entries.
     *
     * @return true if the input has more entries, false otherwise.
     */
    public boolean hasNext() {
        if (hasNext) return true;
        for (int readerIndex : readersWithMoreEntries) {
            if (!nextLoadedForReader[readerIndex]) {
                // the reader at position readerIndex was used in the previous next
                activeIndex = readerIndex;

                final AlignmentReader reader = readers[activeIndex];
                final boolean hasNext = reader.hasNext();
                if (!hasNext) {
                    // reader has no more entries. Remove from further consideration
                    readersWithMoreEntries.remove(activeIndex);

                } else {

                    final Alignments.AlignmentEntry alignmentEntry = reader.next();
                    nextLoadedForReader[readerIndex] = true;
                    final Bucket bucket = buckets[readerIndex];
                    bucket.entry = alignmentEntry;
                    entryHeap.enqueue(bucket);
                    //    System.out.println("entryHeap.size()" + entryHeap.size());
                }
            }
        }
        return (hasNext = (!entryHeap.isEmpty()));

    }

    /**
     * Returns the next alignment entry from the input stream.
     *
     * @return the alignment read entry from the input stream.
     */
    public Alignments.AlignmentEntry next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        } else {

            Bucket bucket = entryHeap.dequeue();
            nextLoadedForReader[bucket.readerIndex] = false;

            //   minEntry = null;
            hasNext = false;
            final Alignments.AlignmentEntry alignmentEntry = bucket.entry;
            final int newQueryIndex = mergedQueryIndex(alignmentEntry.getQueryIndex());
            if (adjustQueryIndices) {
                return alignmentEntry.newBuilderForType().mergeFrom(alignmentEntry).setQueryIndex(newQueryIndex).build();
            } else {
                return alignmentEntry;
            }

        }
    }

    boolean hasNext;
}