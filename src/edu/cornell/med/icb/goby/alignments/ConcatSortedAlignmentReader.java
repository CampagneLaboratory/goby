/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

import it.unimi.dsi.fastutil.AbstractPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectArrayPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;

import java.io.IOException;
import java.util.Comparator;
import java.util.NoSuchElementException;

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

    public ConcatSortedAlignmentReader(final String... basenames) throws IOException {
        super(basenames);
        init(basenames);
    }

    private void init(final String... basenames) {
        nextLoadedForReader = new boolean[basenames.length];


        final Comparator<Bucket> bucketComparator = new Comparator<Bucket>() {
            public int compare(final Bucket bucket, final Bucket bucket1) {
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

    public ConcatSortedAlignmentReader(final boolean adjustQueryIndices, final String... basenames) throws IOException {
        super(adjustQueryIndices, basenames);
        init(basenames);
    }

    /**
     * Skip all entries that have position before (targetIndex,position). This method will use the alignment index
     * to skip directly to the closest chunk start before the entry identified by targetIndex and position.
     *
     * @param targetIndex The index of the target sequence to skip to.
     * @param position    The position on the target sequence.
     * @return The next entry, at position or past position (if not entry at position is found).
     * @throws IOException If an error occurs reading the alignment header. The header is accessed to check that the alignment is sorted.
     */
    public final Alignments.AlignmentEntry skipTo(final int targetIndex, final int position) throws IOException {
        // remove entries from heap if they are located before the skipTo location:
        {
            Bucket bucket;
            while (!entryHeap.isEmpty()) {
                bucket = entryHeap.first();
                if (bucket.entry.getTargetIndex() < targetIndex || bucket.entry.getPosition() < position) {
                    // the first entry in the heap has locaton before the skip to location. We remove it.
                    Bucket removed=entryHeap.dequeue();
                    nextLoadedForReader[removed.readerIndex]=false;
                } else {
                    // the first entry in the heap is at or after the skipTo location. We are done cleaning up the heap.
                    break;
                }
            }
        }
        // populate the heap with the next entry at or past the skipTo position:
        for (final int readerIndex : readersWithMoreEntries) {
            if (!nextLoadedForReader[readerIndex]) {
                // the reader at position readerIndex was used in the previous next
                activeIndex = readerIndex;

                final AlignmentReader reader = readers[activeIndex];
                final Alignments.AlignmentEntry alignmentEntry = reader.skipTo(targetIndex, position);
                if (alignmentEntry == null) {
                    // reader has no more entries. Remove from further consideration
                    readersWithMoreEntries.remove(activeIndex);

                } else {


                    nextLoadedForReader[readerIndex] = true;
                    final Bucket bucket = buckets[readerIndex];
                    bucket.entry = alignmentEntry;
                    bucket.readerIndex=readerIndex;
                    entryHeap.enqueue(bucket);
                    //    System.out.println("entryHeap.size()" + entryHeap.size());
                }
            }
        }
        // return the next entry from the heap:
        if (entryHeap.isEmpty()) return null;
        
        final Bucket bucket = entryHeap.dequeue();
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

    class Bucket {
        Alignments.AlignmentEntry entry;
        int readerIndex;

        public Bucket(final Alignments.AlignmentEntry alignmentEntry, final int index) {
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
    @Override
    public boolean hasNext() {
        if (hasNext) {
            return true;
        }
        for (final int readerIndex : readersWithMoreEntries) {
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
     * Get the index of the alignmnet reader that provided the entry. The entry returned by the previous call to next()
     * originated from the AlignmentReader at index 'readerIndex' in the list provided to the constructor.
     *
     * @return readerIndex.
     */
    public int getReaderIndex() {
        return activeIndex;
    }

    /**
     * Returns the next alignment entry from the input stream.
     *
     * @return the alignment read entry from the input stream.
     */
    @Override
    public Alignments.AlignmentEntry next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        } else {

            final Bucket bucket = entryHeap.dequeue();
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
