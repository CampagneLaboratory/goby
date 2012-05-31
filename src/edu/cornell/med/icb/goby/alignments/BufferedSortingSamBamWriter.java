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

import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.Comparator;

/**
 * A helper class that buffers samRecords up to a certain capacity and ensures
 * records in the buffer are sorted by genomic position before they are written to disk.
 *
 * @author Fabien Campagne
 *         Date: 5/26/12
 *         Time: 4:51 PM
 */
public class BufferedSortingSamBamWriter implements SAMFileWriter {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(BufferedSortingSamBamWriter.class);

    private final SAMFileWriter delegate;
    private final ObjectHeapPriorityQueue<SAMRecord> heap;
    private static final int DEFAULT_CAPACITY = 1000;
    private static final Comparator<net.sf.samtools.SAMRecord> GENOMIC_POSITION_COMPARATOR = new SamRecordGenomicPositionComparator();
    private int capacity;
    /**
     * The maximum target index seen so far.
     */
    private int frontTargetIndex;
    /**
     * The maximum position seen so far on the targetIndex.
     */
    private int frontPosition;
    private boolean check = true;

    public BufferedSortingSamBamWriter(final SAMFileWriter destination, final int capacity) {
        this.capacity = capacity;
        this.delegate = destination;
        this.heap = new ObjectHeapPriorityQueue<net.sf.samtools.SAMRecord>(capacity, GENOMIC_POSITION_COMPARATOR);
    }

    public BufferedSortingSamBamWriter(final SAMFileWriter destination) {
        this(destination, DEFAULT_CAPACITY);
    }
       @Override
    public void addAlignment(net.sf.samtools.SAMRecord entry) {

        while (heap.size() > capacity) {
            final net.sf.samtools.SAMRecord queueEntry = heap.dequeue();
            checkFront(queueEntry);
            delegate.addAlignment(queueEntry);
        }
        heap.enqueue(entry);
    }

    private void checkFront(SAMRecord entry) {
        if (check) {
            final int targetIndex = entry.getReferenceIndex();
            final int position = entry.getAlignmentStart();
            if (targetIndex < frontTargetIndex || targetIndex == frontTargetIndex && position < frontPosition) {
                // we detected an entry that occurs before the front of dequeued entries. We failed to restore sort order
                // we mark the destination writer as non-sorted and inform the end user with a warning.

                LOG.warn("Local sorting strategy failed to restore sort order. The destination will be unsorted. You must sort the output manually to restore sort order.");
                // no need to check anymore:
                check = false;
            }
            if (frontTargetIndex != targetIndex) {
                frontPosition = 0;
            }
            frontTargetIndex = Math.max(frontTargetIndex, targetIndex);
            frontPosition = Math.max(frontPosition, position);
        }
    }





    @Override
    public SAMFileHeader getFileHeader() {
        return delegate.getFileHeader();
    }
    @Override
    public void close()  {
        while (!heap.isEmpty()) {
            final net.sf.samtools.SAMRecord queueEntry = heap.dequeue();
            checkFront(queueEntry);
            delegate.addAlignment(queueEntry);
        }
        delegate.close();
    }


}
