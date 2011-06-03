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

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.ConcatSortedAlignmentReader;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.WarningCounter;
import it.unimi.dsi.fastutil.ints.IntArrayFIFOQueue;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectListIterator;
import it.unimi.dsi.lang.MutableString;
import org.apache.log4j.Logger;

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

    int currentTargetIndex = -1;
    private int numTargets;

    private int processedCount;
    private int numEntriesRealigned;

    @Override
    public int getModifiedCount() {
        return numEntriesRealigned;
    }

    @Override
    public int getProcessedCount() {
        return processedCount;
    }

    /**
     * The targetIndex that was active and for which we may still have entries stored in the corresponding pool:
     */
    private int previousActiveTargetIndex = -1;
    /**
     * The FIFO queue that holds target indices that have entries pooled in tinfo:
     */
    private IntArrayFIFOQueue activeTargetIndices = new IntArrayFIFOQueue();
    private WarningCounter genomeNull = new WarningCounter(2);


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

    int enqueuedCount = 0;

    public Alignments.AlignmentEntry nextRealignedEntry(final int targetIndex, final int position) throws IOException {

        boolean mustLoadPool;
        if (activeTargetIndices.isEmpty()) {
            // nothing seen yet, load the pool
            mustLoadPool = true;
        } else {
        //determine if the pool has enough entry within windowSize:
            InfoForTarget backTargetInfo = targetInfo.get(activeTargetIndices.firstInt());
            // windowLength is zero at the beginning
            int wl=windowLength ==0?1000:windowLength;
            mustLoadPool=backTargetInfo.entriesInWindow.isEmpty() ||
                    backTargetInfo.maxEntryPosition<backTargetInfo.windowStartPosition+wl;
        }

        if (mustLoadPool) {
            // fill the window pool only if we don't have enough entries to return already.
            int windowStartPosition = Integer.MAX_VALUE;
            int minTargetIndex = Integer.MAX_VALUE;
            // push entries to the pool until we have a pool windowLength wide:
            Alignments.AlignmentEntry entry;
            do {
                entry = iterator.skipTo(targetIndex, position);

                if (entry != null) {

                    minTargetIndex = Math.min(minTargetIndex, entry.getTargetIndex());
                    final int entryTargetIndex = entry.getTargetIndex();
                    // push new targetIndices to activeTargetIndices:
                    if (activeTargetIndices.isEmpty() || activeTargetIndices.lastInt() != entryTargetIndex) {

                        activeTargetIndices.enqueue(entryTargetIndex);
                    }
                    final InfoForTarget frontInfo = reallocateTargetInfo(entryTargetIndex);
                    // System.out.printf("windowStartPosition=%,d %n",windowStartPosition);

                    pushEntryToPool(frontInfo, position, entry);

                    windowStartPosition = frontInfo.windowStartPosition;
                }
            }
            while (entry != null && entry.getPosition() < windowStartPosition + windowLength);
        }
        // check if we still have entries in the previously active target:
        int backTargetIndex = activeTargetIndices.firstInt();
        if (targetInfo.get(backTargetIndex).entriesInWindow.isEmpty()) {
            activeTargetIndices.dequeueInt();

            //     targetInfo.remove(0);
            targetInfo.get(0).clear();
            if (activeTargetIndices.isEmpty()) {
                // no more targets, we are done.
                return null;
            }
            backTargetIndex = activeTargetIndices.firstInt();
        }

        final InfoForTarget backInfo = targetInfo.get(backTargetIndex); // the pool info we will use to dequeue the entry at the back of the window
        //  System.out.printf("back is holding %d entries %n", backInfo.entriesInWindow.size());
        // now find the entry at the left of the realignment window on the active target:
        if (backInfo.entriesInWindow.isEmpty()) {
            return null;
        }
        Alignments.AlignmentEntry returnedEntry = backInfo.entriesInWindow.remove();

        if (backInfo.positionsWithSpanningIndel.size() > 0) {
            returnedEntry = realign(returnedEntry, backInfo);
        }
        // advance the windowStartPosition
        int previousWindowStart = backInfo.windowStartPosition;
        int windowStartPosition = Math.max(backInfo.windowStartPosition, returnedEntry.getPosition());
        // remove indel locations that are now outside the new window position
        if (previousWindowStart != windowStartPosition) {
            int lastPosition = windowStartPosition - 1;

            backInfo.removeIndels(previousWindowStart, lastPosition);
        }

        ++processedCount;
        return returnedEntry;

    }

    private InfoForTarget reallocateTargetInfo(int targetIndex) {
        int intermediateTargetIndex = targetInfo.size() - 1;
        while (intermediateTargetIndex <= targetIndex) {

            targetInfo.add(new InfoForTarget(++intermediateTargetIndex));
            numTargets = targetInfo.size();
        }
        return targetInfo.get(targetIndex);

    }

    private final boolean[] directions = new boolean[]{true, false};

    private Alignments.AlignmentEntry realign(Alignments.AlignmentEntry entry, InfoForTarget tinfo) {
        int currentBestScore = 0;
        ObservedIndel bestScoreIndel = null;
        boolean bestScoreDirection = false;
        for (ObservedIndel indel : tinfo.potentialIndels) {
            for (boolean direction : directions) {

                final int realignedScore = score(entry, indel, direction, currentBestScore, genome);
                if (realignedScore > currentBestScore) {
                    currentBestScore = realignedScore;
                    bestScoreIndel = indel;
                    bestScoreDirection = direction;
                }
            }
        }
        if (currentBestScore == 0) {
            return entry;
        } else {
            // actually modify entry to realign through the indel:
            ++numEntriesRealigned;
            return realign(entry, bestScoreIndel, bestScoreDirection, currentBestScore);
        }

    }

    private Alignments.AlignmentEntry realign(Alignments.AlignmentEntry entry,
                                              ObservedIndel indel,
                                              boolean shiftForward, int scoreDelta) {
        // use entry as prototype:
        Alignments.AlignmentEntry.Builder builder = Alignments.AlignmentEntry.newBuilder(entry);
        // update the score to reflect the realignment:
        builder.setScore(entry.getScore() + scoreDelta);
        int indelLength = indel.length();
        int entryPosition = entry.getPosition();
        final int originalEntryPosition = entryPosition;
        if (!shiftForward && indel.isReadInsertion()) {
            // shifting to the left a read insertion, must shift alignment start position to the right
            entryPosition = entry.getPosition() - indelLength;
            builder.setPosition(entryPosition);

            //  when entry position changes, we need to update the pool to reflect the new entry sort order
            // this is accomplished by wrapping this processor with a LocalSortProcessor instance.
        }
        int indelOffsetInAlignment = indel.getStart() - entryPosition;

        int varCount = entry.getSequenceVariationsCount();
        int targetIndex = entry.getTargetIndex();
        int score = 0;
        int direction = shiftForward ? 1 : -1;
        if (genome == null) {
            genomeNull.warn(LOG, "Genome must not be null outside of Junit tests.");
            return entry;
        }
        /*
         *Reference positions for which the alignment does not agree with the reference, 0-based:
         */
        IntArraySet variantPositions = new IntArraySet();
        // determine if rewrittenVariations become compatible with reference when indel is introduced in this alignment:
        // increase the score by 1 for every base that beomes compatible.
        ObjectArrayList<Alignments.SequenceVariation> rewrittenVariations = new ObjectArrayList<Alignments.SequenceVariation>();
        for (int i = 0; i < varCount; i++) {
            Alignments.SequenceVariation var = entry.getSequenceVariations(i);
            // check if var becomes compatible with reference when var's varPosition is shifted by the length of the indel in the specified shiftForward
            // newGenomicPosition is zero-based
            final int newGenomicPosition = var.getPosition() + (direction * indelLength) + originalEntryPosition - 1;
            for (int j = 0; j < var.getTo().length(); ++j) {

                final char toBase = var.getTo().charAt(j);
                final int index = newGenomicPosition + j;
                if (index < 0 || index > genome.getLength(targetIndex)) {
                    score += -10;
                } else {
                    final boolean compatible = genome.get(targetIndex, newGenomicPosition + j) == toBase;
                    if (!compatible) {
                        // we keep only sequence rewrittenVariations that continue to be  incompatible with the reference after inserting the indel:

                        rewrittenVariations.add(var);
                    }

                    variantPositions.add(var.getPosition() + entryPosition + j - 1);
                }
            }

        }

        // Determine which previously unchanged bases become incompatible with the reference when the the indel is introduced
        // startAlignment and endAlignment are zero-based
        int startAlignment = shiftForward ? entryPosition + indelOffsetInAlignment : entryPosition;
        int endAlignment = shiftForward ? entry.getTargetAlignedLength() + entryPosition : indelOffsetInAlignment + entryPosition + (direction * indelLength);
        //  String pre = getGenomeSegment(genome, targetIndex, startAlignment, endAlignment);
        //  String post = getGenomeSegment(genome, targetIndex, startAlignment + (direction * indelLength), endAlignment + (direction * indelLength));
        //   System.out.printf(" pre and post alignments: %n%s\n%s%n", pre, post);
        // pos is zero-based:
        for (int pos = startAlignment; pos < endAlignment; pos++) {
            // both variantPositions and pos are zero-based:
            if (!variantPositions.contains(pos)) {
                // this base matched the reference sequence:
                final int realignedPos = pos + (direction * indelLength);
                // if the realigned varPosition lies outside of the reference penalize heavily with -10, otherwise
                // count -1 for every new mismatch introduced by the indel:
                assert realignedPos >= 0 : "realignedPos cannot be negative for best indel.";


                final char fromBase = genome.get(targetIndex, realignedPos);
                final char toBase = genome.get(targetIndex, pos);
                final boolean compatible = fromBase == toBase;

                if (!compatible) {
                    Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();
                    // varPosition is one-based while realignedPos and entryPos are zero-based:
                    final int varPosition = direction * (realignedPos - entryPosition) + 1;
                    varBuilder.setPosition(varPosition);
                    varBuilder.setFrom(Character.toString(fromBase));
                    varBuilder.setTo(Character.toString(toBase));

                    int readIndex = entry.getMatchingReverseStrand() ?
                            entry.getQueryLength() - indelOffsetInAlignment + (shiftForward ? 1 : indelLength) :
                            varPosition;
                    varBuilder.setReadIndex(readIndex);
                    rewrittenVariations.add(varBuilder.build());

                }
            }
        }
        // finally, add the indel into the revised alignment:
        Alignments.SequenceVariation.Builder varBuilder = Alignments.SequenceVariation.newBuilder();

        //  fix varPosition for negative strand, var positions are one-based:
        final int varPosition = shiftForward ? indelOffsetInAlignment + 1 : indel.getStart() - entryPosition + 1;
        varBuilder.setPosition(varPosition);
        varBuilder.setFrom(indel.from);
        varBuilder.setTo(indel.to);
        // we set readIndex to the left most varPosition before the read gap, by convention, see  http://tinyurl.com/goby-sequence-variations (read deletion)
        int readIndex = entry.getMatchingReverseStrand() ?
                entry.getQueryLength() - indelOffsetInAlignment + (shiftForward ? 1 : indelLength) :
                varPosition;
        varBuilder.setReadIndex(readIndex);
        rewrittenVariations.add(varBuilder.build());
        builder = builder.clearSequenceVariations();
        for (Alignments.SequenceVariation var : rewrittenVariations) {
            builder = builder.addSequenceVariations(var);
        }
        final Alignments.AlignmentEntry alignmentEntry = builder.build();
        //    System.out.printf("realigned queryIndex=%d%n", alignmentEntry.getQueryIndex());
        return alignmentEntry;
    }

    /**
     * Score the realignment of an entry with respect to a potential indel.
     *
     * @param entry            The entry to score as if it was realigned with respect to the indel
     * @param indel            The indel under consideration.
     * @param shiftForward     Whether the indel should be introduced by shifting bases forward
     * @param currentBestScore The current maximum score over a set of indels under consideration.  @return The score obtained when realigning the entry with respect to the provided indel.
     * @param genome           The genome to use to lookup reference bases
     * @return The score that would be observed if the indel was inserted into the alignment represented by entry.
     */

    public final int score(final Alignments.AlignmentEntry entry, final ObservedIndel indel, final boolean shiftForward, final int currentBestScore,
                           final RandomAccessSequenceInterface genome) {
        int entryPosition = entry.getPosition();
        int indelOffsetInAlignment = indel.getStart() - entryPosition;
        int indelLength = indel.length();
        int varCount = entry.getSequenceVariationsCount();
        int targetIndex = entry.getTargetIndex();
        int score = 0;
        int direction = shiftForward ? 1 : -1;

        if (genome == null) {
            genomeNull.warn(LOG, "Genome must not be null outside of JUnit tests.");
            return Integer.MIN_VALUE;
        }
        final int targetLength = genome.getLength(targetIndex);
        /*
         *Reference positions for which the alignment does not agree with the reference, 0-based:
         */
        IntArraySet variantPositions = new IntArraySet();
        // determine if variations become compatible with reference when indel is introduced in this alignment:
        // increase the score by 1 for every base that beomes compatible.
        for (int i = 0; i < varCount; i++) {
            Alignments.SequenceVariation var = entry.getSequenceVariations(i);
            // check if var becomes compatible with reference when var's position is shifted by the length of the indel in the specified shiftForward
            // newGenomicPosition is zero-based 
            // final int originalGenomicPosition = var.getPosition() + entryPosition - 1;
            final int newGenomicPosition = var.getPosition() + (direction * indelLength) + entryPosition - 1;
            for (int j = 0; j < var.getTo().length(); ++j) {

                final char toBase = var.getTo().charAt(j);
                final int index = newGenomicPosition + j;
                if (index < 0 || index > genome.getLength(targetIndex)) {
                    score += -10;
                } else {
                    final boolean compatible = genome.get(targetIndex, newGenomicPosition + j) == toBase;

                    score += compatible ? 1 : 0;
                    // store which reference positions are different from the reference:

                    variantPositions.add(var.getPosition() + entryPosition + j - 1);
                }
            }

        }
        // Determine which previously unchanged bases become incompatible with the reference when the the indel is introduced
        // Consider the span of reference between the indel insertion point and the end of the reference alignment going in the direction of extension.
        // startAlignment and endAlignment are zero-based
        int startAlignment = shiftForward ? entryPosition + indelOffsetInAlignment : entryPosition;
        int endAlignment = shiftForward ? entry.getTargetAlignedLength() + entryPosition : indelOffsetInAlignment + entryPosition + (direction * indelLength);


        //       String pre = getGenomeSegment(genome, targetIndex, startAlignment, endAlignment);
//

//         String post = getGenomeSegment(genome, targetIndex, startAlignment + indelLength, endAlignment + indelLength);
        //   System.out.printf(" pre and post alignments: %n%s\n%s%n", pre, post);
        // pos is zero-based:
        for (int pos = startAlignment; pos < endAlignment; pos++) {
            // both variantPositions and pos are zero-based:
            if (!variantPositions.contains(pos)) {
                // this base matched the reference sequence:
                final int realignedPos = pos + (direction * indelLength);
                // if the realigned position lies outside of the reference penalize heavily with -10, otherwise
                // count -1 for every new mismatch introduced by the indel:

                if (realignedPos < 0 || realignedPos >= targetLength) {
                    score += -10;
                } else {
                    final char refBase = genome.get(targetIndex, pos);
                    final char newRefBase = genome.get(targetIndex, realignedPos);

                    score += (refBase == newRefBase) ? 0 : -1;
                }
            }
            // System.out.printf("indelOffsetInAlignment: %d shiftForward: %b score: %d%n", indelOffsetInAlignment, shiftForward, score);
        }

        //   System.out.printf("indelOffsetInAlignment: %d shiftForward: %b score: %d%n", indelOffsetInAlignment, shiftForward, score);
        return score;
    }

    private String getGenomeSegment(RandomAccessSequenceInterface genome, int targetIndex, int startAlignment, int endAlignment) {
        MutableString sequence = new MutableString();
        for (int pos = startAlignment; pos < endAlignment; pos++) {
            sequence.append(genome.get(targetIndex, pos));
        }
        return sequence.toString();
    }

    public void pushEntryToPool(InfoForTarget tinfo, int position, Alignments.AlignmentEntry entry) {
        // the window start is only decreased in the pushing step, never increased.
        final int entryPosition = entry.getPosition();
        tinfo.windowStartPosition = Math.min(tinfo.windowStartPosition, entryPosition);
        tinfo.maxEntryPosition=Math.max(tinfo.maxEntryPosition,entryPosition);
        // set window length to twice the longest read length.
        windowLength = Math.max(windowLength, entry.getQueryLength() * 2);

        // detect if the entry contains an indel. Update the window indel state accordingly
        for (int i = 0; i < entry.getSequenceVariationsCount(); ++i) {
            final Alignments.SequenceVariation var = entry.getSequenceVariations(i);
            if (isIndel(var)) {
                // start and last position are zero-based:           A--CAC start=1 end=3
                final int startPosition = var.getPosition() + entryPosition - 1;
                final int lastPosition = var.getPosition() + entryPosition +
                        Math.max(var.getFrom().length(), var.getTo().length()) - 1;

                tinfo.addIndel(startPosition, lastPosition, var.getFrom(), var.getTo());

            }
        }

        tinfo.entriesInWindow.add(entry);
        ++enqueuedCount;
    }

    private boolean isIndel(Alignments.SequenceVariation var) {
        return (var.getFrom().indexOf('-') >= 0 ||
                var.getTo().indexOf('-') >= 0);
    }


    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    private static final Logger LOG = Logger.getLogger(RealignmentProcessor.class);
}
