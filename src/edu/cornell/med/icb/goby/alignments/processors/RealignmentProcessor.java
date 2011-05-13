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
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.ConcatSortedAlignmentReader;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
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
    int currentTargetIndex = -1;
    private int numTargets;
    /**
     * The targetIndex that was active and for which we may still have entries stored in the corresponding pool:
     */
    private int previousActiveTargetIndex;

    protected class InfoForTarget {
        int targetIndex;
        public IntArraySet positionsWithSpanningIndel = new IntArraySet();
        UnboundedFifoPool<Alignments.AlignmentEntry> entriesInWindow = new UnboundedFifoPool<Alignments.AlignmentEntry>();

        InfoForTarget(int targetIndex) {
            this.targetIndex = targetIndex;
        }
    }

    ObjectArrayList<InfoForTarget> targetInfo = new ObjectArrayList<InfoForTarget>();

    final SkipToIterator iterator;
    private RandomAccessSequenceInterface genome;

    public RealignmentProcessor(final ConcatSortedAlignmentReader sortedReaders) {
        iterator = new SkipToSortedReader(sortedReaders);
        numTargets = sortedReaders.getNumberOfTargets();
        targetInfo = new ObjectArrayList<InfoForTarget>(numTargets);
    }

    public RealignmentProcessor(final ObjectListIterator<Alignments.AlignmentEntry> entryIterator) {
        iterator = new SkipToListIterator(entryIterator);
    }

    public Alignments.AlignmentEntry nextRealignedEntry(final int targetIndex, final int position) throws IOException {

        int minTargetIndex = Integer.MAX_VALUE;
        // push entries to the pool until we have a pool windowLength wide:
        Alignments.AlignmentEntry entry;
        do {
            entry = iterator.skipTo(targetIndex, position);
            if (entry != null) {
                minTargetIndex = Math.min(minTargetIndex, entry.getTargetIndex());
                final int entryTargetIndex = entry.getTargetIndex();
                final InfoForTarget frontInfo = reallocateTargetInfo(entryTargetIndex);
                pushEntryToPool(frontInfo, position, entry);
                if (entryTargetIndex != minTargetIndex) {
                    windowStartPosition = entry.getPosition();
                }
            }
        }
        while (entry != null && entry.getPosition() < windowStartPosition + windowLength);
        // check if we still have entries in the previously active target:

        final InfoForTarget backInfo; // the pool info we will use to dequeue the entry at the back of the window
        final InfoForTarget altBackInfo = previousActiveTargetIndex >= 0 ? targetInfo.get(previousActiveTargetIndex) : null;
        if (altBackInfo == null || altBackInfo.entriesInWindow.isEmpty()) {
            // we are done with the previous active target, update to current target:
            final int currentTargetIndex = targetInfo.size() - 1;
            previousActiveTargetIndex = currentTargetIndex;
            // set backInfo to use the new target
            backInfo = targetInfo.get(currentTargetIndex);
        } else {
            // set backInfo on the previously active target since it still contains entries to dequeue:
            backInfo = altBackInfo;
        }

        // now find the entry at the left of the realignment window on the active target:
        if (backInfo.entriesInWindow.isEmpty()) return null;
        Alignments.AlignmentEntry returnedEntry = backInfo.entriesInWindow.remove();

        if (backInfo.positionsWithSpanningIndel.size() > 0) {
            returnedEntry = realign(returnedEntry);
        }
        // advance the windowStartPosition
        int previousWindowStart = windowStartPosition;
        windowStartPosition = Math.max(windowStartPosition, returnedEntry.getPosition());
        // remove indel locations that are now outside the new window position
        if (previousWindowStart != windowStartPosition) {
            for (int pos = previousWindowStart; pos < windowStartPosition - 1; ++pos) {
                backInfo.positionsWithSpanningIndel.rem(pos);
            }
        }
        return returnedEntry;

    }

    private InfoForTarget reallocateTargetInfo(int targetIndex) {
        int intermediateTargetIndex = targetInfo.size() - 1;
        while (targetIndex >= numTargets) {
            previousActiveTargetIndex = targetInfo.size() - 1;
            targetInfo.add(new InfoForTarget(++intermediateTargetIndex));
            numTargets = targetInfo.size();
        }
        return targetInfo.get(targetIndex);

    }

    private Alignments.AlignmentEntry realign(Alignments.AlignmentEntry entry) {
        return entry;
    }

    public void pushEntryToPool(RealignmentProcessor.InfoForTarget tinfo, int position, Alignments.AlignmentEntry entry) {
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
                    tinfo.positionsWithSpanningIndel.add(p);
                }
            }
        }
        tinfo.entriesInWindow.add(entry);
    }

    private boolean isIndel(Alignments.SequenceVariation var) {
        return (var.getFrom().indexOf('-') != 0 ||
                var.getTo().indexOf('-') != 0);
    }


    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }
}
