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

package edu.cornell.med.icb.goby.alignments.processors;

import com.sun.tools.javac.util.Position;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.fastutil.AbstractPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectArrayPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;

import java.io.IOException;
import java.util.Comparator;

/**
 * An alignment processor that re-sorts entries in a window by increasing genomic position. This processor can be
 * used to locally sort entry start positions after realignment of some entries around indels. This is required
 * because local realignment around indels can change the alignment start position for those entries realigned
 * to the left).
 *
 * @author Fabien Campagne
 *         Date: 6/2/11
 *         Time: 10:49 AM
 */
public class LocalSortProcessor implements AlignmentProcessorInterface {
    AbstractPriorityQueue<Alignments.AlignmentEntry> entryHeap;
    AlignmentProcessorInterface delegate;
    /**
     * We store at most alignments for 30 consecutive genomic positions. 30 is chosen because it is much larger than
     * the maximum indel size that aligners can detect. Realignment cannot change the start position by more than the
     * indel length, so there is no need for a larger window.
     */
    private int windowLength = 30;
    private int processedCount = 0;
    private int modifiedCount = 0;
    final Comparator<? super Alignments.AlignmentEntry> comparator = new Comparator<Alignments.AlignmentEntry>() {
        @Override
        public int compare(final Alignments.AlignmentEntry alignmentEntry, final Alignments.AlignmentEntry alignmentEntry1) {
            final int targetDiff = alignmentEntry.getTargetIndex() - alignmentEntry1.getTargetIndex();
            if (targetDiff != 0) {
                return targetDiff;
            } else {
                final int value = alignmentEntry.getPosition() - alignmentEntry1.getPosition();
                return value;
            }
        }
    };
    private int lastCurrentPosition;
    private int lastTargetIndex;

    public LocalSortProcessor(final AlignmentProcessorInterface delegate) {
        this.delegate = delegate;
        // heap will hold up to 10,000 elements initially, but will grow as needed:

        entryHeap =  new ObjectHeapPriorityQueue<Alignments.AlignmentEntry>(10000, comparator);
    }

    /**
     * Target index for the front of the window.
     */
    int currentTargetIndex = Integer.MAX_VALUE;
    /**
     * Position (zero-based) for the front of the window.
     */
    int currentPosition = 0;
    boolean finished = false;

    @Override
    public Alignments.AlignmentEntry nextRealignedEntry(final int targetIndex, final int position) throws IOException {
        if (finished && entryHeap.isEmpty()) {
            return null;
        }
        boolean mustLoadPool = false;
        if (!finished) {
            if (entryHeap.isEmpty()) {
                // nothing seen yet, load the pool
                mustLoadPool = true;
            } else {
                //determine if the pool has enough entry within windowSize:
                Alignments.AlignmentEntry firstEntry = entryHeap.first();

                // windowLength is zero at the beginning

                mustLoadPool = entryHeap.isEmpty() || currentTargetIndex == firstEntry.getTargetIndex() &&
                        currentPosition - windowLength < firstEntry.getPosition();
            }
        }
        if (mustLoadPool) {
            Alignments.AlignmentEntry entry;
          final int initialCurrentPosition=currentPosition;
            do {
                entry = delegate.nextRealignedEntry(targetIndex, position);
                if (entry != null) {
                    currentTargetIndex = Math.min(currentTargetIndex, entry.getTargetIndex());

                    if (currentTargetIndex != entry.getTargetIndex()) {
                        if (!entryHeap.isEmpty()) {
                            pushEntry(entry);
                            break;
                        }
                        //  assert entryHeap.isEmpty() : "entryHeap must be empty when enqueing the first entry of a new target";
                    }
                    pushEntry(entry);
                } else {
                    finished = true;
                }
            } while (entry != null && entry.getPosition() < initialCurrentPosition + windowLength);

        }
        if (entryHeap.isEmpty()) {
            finished = true;
            return null;
        } else {

            Alignments.AlignmentEntry entry = entryHeap.dequeue();
            //    System.out.println("dequeuing: " + entry);
            ++processedCount;
            // update the position of the front:

            currentPosition = entry.getPosition();
            currentTargetIndex = entry.getTargetIndex();

            return entry;
        }

    }

    private void pushEntry(final Alignments.AlignmentEntry entry) {
    //    System.out.println("Position=" + entry.getPosition());
        if (entry.getTargetIndex() == lastTargetIndex && entry.getPosition() <lastCurrentPosition) {
        //    System.out.println("entry to be enqueued needs to be resorted.");
            // entry to be enqueued will need to be resorted.
            ++modifiedCount;
        }
        //  currentPosition = Math.min(currentPosition, entry.getPosition());

        //   System.out.println("enqueuing " + entry);
        entryHeap.enqueue(entry);
        lastCurrentPosition=Math.max(lastCurrentPosition,entry.getPosition());
        lastTargetIndex=Math.max(lastTargetIndex,entry.getTargetIndex());
    }

    @Override
    public void setGenome(RandomAccessSequenceInterface genome) {
        delegate.setGenome(genome);
    }

    @Override
    public int getModifiedCount() {
   //     System.out.printf("Have locally sorted %d entries%n",modifiedCount);
        return delegate.getModifiedCount();
    }

    @Override
    public int getProcessedCount() {
        return delegate.getProcessedCount();

    }
}
