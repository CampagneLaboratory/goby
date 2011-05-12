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

import edu.cornell.med.icb.goby.algorithmic.data.UnboundedFifoPool;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.alignments.processors.AlignmentProcessorInterface;
import edu.cornell.med.icb.goby.alignments.*;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.objects.ObjectListIterator;

import java.io.IOException;

/**
 * Support to realign reads on the fly in the proximity of indels.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2011
 *         Time: 11:58:07 AM
 */
public class RealignmentProcessor implements AlignmentProcessorInterface {
    int windowLength = 0;
    int windowStartPosition = Integer.MAX_VALUE;
    UnboundedFifoPool<Alignments.AlignmentEntry> entriesInWindow = new UnboundedFifoPool<Alignments.AlignmentEntry>();
    int currentTargetIndex = -1;

    protected IntArraySet positionsWithSpanningIndel = new IntArraySet();
    final SkipToIterator iterator;
    private RandomAccessSequenceInterface genome;

    public RealignmentProcessor(final ConcatSortedAlignmentReader sortedReaders) {
        iterator = new SkipToSortedReader(sortedReaders);
    }

    public RealignmentProcessor(final ObjectListIterator<Alignments.AlignmentEntry> entryIterator) {
        iterator = new SkipToListIterator(entryIterator);
    }

    public Alignments.AlignmentEntry nextRealignedEntry(final int targetIndex, final int position) throws IOException {

        Alignments.AlignmentEntry entry;
        do {
            entry = iterator.skipTo(targetIndex, position);
            if (entry != null) {
                pushEntryToPool(position, entry);

            }
        }
        while (entry != null && entry.getPosition() < windowStartPosition + windowLength);

        if (entriesInWindow.isEmpty()) return null;
        Alignments.AlignmentEntry returnedEntry = entriesInWindow.remove();
        if (returnedEntry.getTargetIndex() != currentTargetIndex) {
            // this is the first entry in a new reference sequence.
            positionsWithSpanningIndel.clear();
        }

        if (positionsWithSpanningIndel.size() > 0) {
            returnedEntry = realign(returnedEntry);
        }
        // advance the windowStartPosition
        int previousWindowStart = windowStartPosition;
        windowStartPosition = Math.max(windowStartPosition, returnedEntry.getPosition());
        // remove indel locations that are now outside the new window position
        if (previousWindowStart != windowStartPosition) {
            for (int pos = previousWindowStart; pos < windowStartPosition - 1; ++pos) {
                positionsWithSpanningIndel.rem(pos);
            }
        }
        return returnedEntry;

    }

    private Alignments.AlignmentEntry realign(Alignments.AlignmentEntry entry) {
        return entry;
    }

    public void pushEntryToPool(int position, Alignments.AlignmentEntry entry) {
        windowStartPosition = Math.min(windowStartPosition, position);
        windowStartPosition = Math.min(windowStartPosition, entry.getPosition());
        // set window length to twice the longest read length.
        windowLength = Math.max(windowLength, entry.getQueryLength() * 2);

        // detect if the entry contains an indel. Update the window indel state accordingly
        for (int i = 0; i < entry.getSequenceVariationsCount(); ++i) {
            final Alignments.SequenceVariation var = entry.getSequenceVariations(i);
            if (isIndel(var)) {
                final int startPosition = var.getPosition() + entry.getPosition();
                final int lastPosition = var.getPosition() + entry.getPosition() +
                        Math.max(var.getFrom().length(), var.getTo().length());

                for (int p = startPosition; p < lastPosition; p++) {
                    positionsWithSpanningIndel.add(p);
                }
            }
        }
        entriesInWindow.add(entry);
    }

    private boolean isIndel(Alignments.SequenceVariation var) {
        return (var.getFrom().indexOf('-') != 0 ||
                var.getTo().indexOf('-') != 0);
    }


    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }
}
