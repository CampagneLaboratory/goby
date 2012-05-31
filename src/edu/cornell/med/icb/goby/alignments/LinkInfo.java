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

import it.unimi.dsi.bits.Fast;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;
import java.util.Collections;

/**
 * @author Fabien Campagne
 *         Date: 3/7/12
 *         Time: 5:35 PM
 */
class LinkInfo {
    private final Alignments.RelatedAlignmentEntry EMPTY_RELATED_ENTRY_LINK = Alignments.RelatedAlignmentEntry.newBuilder().build();

    private final IntList hasLinks = new IntArrayList();
    private final IntList deltaPositions = new IntArrayList();
    private final IntList deltaTargetIndices = new IntArrayList();
    private final IntList fragmentIndices = new IntArrayList();
    private int previousPositionPositiveStrand;
    private int previousPositionNegativeStrand;
    private int previousTargetIndex;
    int index = 0;
    private final String label;
    private final String HAS_LINK_LABEL;
    private final String POS_LABEL;
    private final String TARGETS_LABEL;
    private final String FRAGMENTS_LABEL;
    private final AlignmentCollectionHandler handler;
    private int fragmentIndex;
    private Int2ObjectMap<IntArrayList> queryIndex2PositionList;
    private Int2ObjectMap<IntArrayList> queryIndex2FragmentIndexList;


    LinkInfo(final AlignmentCollectionHandler handler, final String label) {
        this.handler = handler;
        this.label = label;
        this.queryIndex2PositionList = handler.queryIndexToPositionList;
        this.queryIndex2FragmentIndexList = handler.queryIndex2FragmentIndices;
        HAS_LINK_LABEL = "has" + label;
        POS_LABEL = label + "-pos";
        TARGETS_LABEL = label + "-targets";
        FRAGMENTS_LABEL = label + "-fragments";
    }

    public void reset() {
        hasLinks.clear();
        deltaPositions.clear();
        deltaTargetIndices.clear();
        fragmentIndices.clear();
        index = 0;
        deltaPositionIndex = 0;
        fragmentIndex = 0;
        previousPositionPositiveStrand = 0;
        previousPositionNegativeStrand = 0;
    }

    public Alignments.RelatedAlignmentEntry code(final boolean hasLink, final Alignments.AlignmentEntry entry,
                                                 final Alignments.RelatedAlignmentEntry relatedLink) {

        hasLinks.add(hasLink ? 1 : 0);

        //return null;
        if (!hasLink) {
            return null;
        }
        if (handler.enableDomainOptimizations) {
            if (canOptimize(entry, relatedLink)) {
                return null;
            }
        }
        final int position = relatedLink.getPosition();
        final int targetIndex = relatedLink.getTargetIndex();
        final Alignments.RelatedAlignmentEntry.Builder result = Alignments.RelatedAlignmentEntry.newBuilder(relatedLink);
        boolean justResetPos = false;
        if (index > 0 && targetIndex == previousTargetIndex) {
            result.clearPosition();
            result.clearTargetIndex();

        } else {
            justResetPos = true;
        }
        boolean entryMatchingReverseStrand = entry.getMatchingReverseStrand();
        if (!justResetPos) {

            final int deltaPairPos = position - (entryMatchingReverseStrand ? previousPositionNegativeStrand : previousPositionPositiveStrand);
            deltaPositions.add(deltaPairPos);
            deltaTargetIndices.add(targetIndex - previousTargetIndex);
        }
        if (entryMatchingReverseStrand) {
            previousPositionNegativeStrand = position;
        } else {
            previousPositionPositiveStrand = position;
        }
        previousTargetIndex = targetIndex;
        fragmentIndices.add(relatedLink.getFragmentIndex());
        result.clearFragmentIndex();
        index++;
        final Alignments.RelatedAlignmentEntry built = result.buildPartial();
        if (EMPTY_RELATED_ENTRY_LINK.equals(built)) {
            // no other filed in this link, remove from transformed.
            return null;
        } else {
            return built;
        }

    }

    private boolean canOptimize(Alignments.AlignmentEntry entry, Alignments.RelatedAlignmentEntry relatedLink) {
        final IntArrayList positionList = queryIndex2PositionList.get(entry.getQueryIndex());
        final IntArrayList fragmentList = queryIndex2FragmentIndexList.get(entry.getQueryIndex());
        assert positionList != null : "positionList must be found for queryIndex=" + entry.getQueryIndex();
        final int position = entry.getPosition();
        final int linkPosition = relatedLink.getPosition();
        int countForPos = 0;
        int countForLinkPos = 0;
        int index = 0;
        int thisEntryIndex = -1;
        int linkIndex = -1;
        final int entryFragmentIndex = entry.getFragmentIndex();
        final int linkFragmentIndex = relatedLink.getFragmentIndex();
        for (final int pos : positionList) {

            if (pos == position && fragmentList.get(index)== entryFragmentIndex) {
                countForPos++;
                thisEntryIndex = index;
            }
            if (pos == linkPosition && fragmentList.get(index)== linkFragmentIndex) {
                countForLinkPos++;
                linkIndex = index;
            }
            index += 1;
        }

        //final int thisEntryIndex = Collections.binarySearch(positionList, position);

        if (countForPos > 1 || countForLinkPos > 1 || entry.getTargetIndex() != relatedLink.getTargetIndex()
                || linkIndex < 0 ||thisEntryIndex<0) {
            handler.linkOffsetOptimization.add(AlignmentCollectionHandler.MISSING_VALUE);
            return false;
        } else {
            final int linkOffset = Fast.int2nat(linkIndex - thisEntryIndex);
            // linkIndex=Fast.nat2int(linkOffset)+thisEntryIndex
            handler.linkOffsetOptimization.add(linkOffset);
            return true;
        }


    }

    private Alignments.RelatedAlignmentEntry wasOptimized(Alignments.AlignmentEntry.Builder entry, Alignments.RelatedAlignmentEntry source) {

        final int linkOffset = handler.getNextLinkOptimizationOffset();
        if (linkOffset == AlignmentCollectionHandler.MISSING_VALUE) {
            return null;
        } else {

            Alignments.RelatedAlignmentEntry.Builder result = source != null ? Alignments.RelatedAlignmentEntry.newBuilder(source) :
                    Alignments.RelatedAlignmentEntry.newBuilder();
            //    final int thisEntryIndex = Collections.binarySearch(list, entry.getPosition());
            //   final int linkIndex = Fast.nat2int(linkOffset) + thisEntryIndex;
            // store linkOffset temporarily in the optimized index field. We will resolve this to actual position and fragment index
            // after all entries have been decoded.
            result.setOptimizedIndex(Fast.nat2int(linkOffset));
            result.setTargetIndex(entry.getTargetIndex());
            return result.buildPartial();
        }

    }

    int deltaPositionIndex;

    public Alignments.RelatedAlignmentEntry decode(final int originalIndex, final Alignments.AlignmentEntry.Builder entry,
                                                   final Alignments.RelatedAlignmentEntry source) {
        if (hasLinks.getInt(originalIndex) == 1) {
            final Alignments.RelatedAlignmentEntry optimizedAway = wasOptimized(entry, source);
            if (optimizedAway != null) {
                return optimizedAway;
            }
            final boolean entryMatchingReverseStrand = entry.hasMatchingReverseStrand() ? entry.getMatchingReverseStrand() : false;

            Alignments.RelatedAlignmentEntry.Builder result = source != null ? Alignments.RelatedAlignmentEntry.newBuilder(source) :
                    Alignments.RelatedAlignmentEntry.newBuilder();

            result.setFragmentIndex(fragmentIndices.getInt(fragmentIndex++));
            boolean justResetPos = false;
            int position = 0;
            int targetIndex = 0;
            if (source != null && (source.hasPosition() || source.hasTargetIndex())) {
                position = source.getPosition();
                targetIndex = source.getTargetIndex();
                justResetPos = true;
                //  System.out.println("just reset pos");

            }
            if (index > 0 && !justResetPos) {

                // keep track of previous positions for each strand the entry is on. This helps because these positions
                // are typically shifted by about insert size
                final int previousPosition = entryMatchingReverseStrand ? previousPositionNegativeStrand : previousPositionPositiveStrand;
                position = deltaPositions.getInt(deltaPositionIndex) + previousPosition;
                targetIndex = deltaTargetIndices.getInt(deltaPositionIndex) + previousTargetIndex;
                deltaPositionIndex++;
                //System.out.println("getting pos from deltas, position="+position);

            }
            result.setPosition(position);
            result.setTargetIndex(targetIndex);
            if (entryMatchingReverseStrand) {
                previousPositionNegativeStrand = position;
            } else {
                previousPositionPositiveStrand = position;
            }
            previousTargetIndex = targetIndex;
            ++index;
            return result.build();
        } else {
            ++index;
            return null;
        }
    }


    public void write(final OutputBitStream out) throws IOException {

        handler.writeArithmetic(whenDebug(HAS_LINK_LABEL), hasLinks, out);
        //    System.out.printf("delta-links-%s n=%d %s %n", label, deltaPositions.size(), deltaPositions.toString());
        handler.writeArithmetic(whenDebug(POS_LABEL), deltaPositions, out);
        //    handler.writeRiceCoding(whenDebug(POS_LABEL), deltaPositions, out);
        handler.writeArithmetic(whenDebug(TARGETS_LABEL), deltaTargetIndices, out);
        handler.writeArithmetic(whenDebug(FRAGMENTS_LABEL), fragmentIndices, out);
    }

    public void read(final int numEntriesInChunk, final InputBitStream bitInput) throws IOException {
        handler.decodeArithmetic(whenDebug(HAS_LINK_LABEL), numEntriesInChunk, bitInput, hasLinks);
        handler.decodeArithmetic(whenDebug(POS_LABEL), numEntriesInChunk, bitInput, deltaPositions);
        //    System.out.printf("delta-links-%s n=%d %s %n", label, deltaPositions.size(), deltaPositions.toString());
        handler.decodeArithmetic(whenDebug(TARGETS_LABEL), numEntriesInChunk, bitInput, deltaTargetIndices);
        handler.decodeArithmetic(whenDebug(FRAGMENTS_LABEL), numEntriesInChunk, bitInput, fragmentIndices);
    }

    private String whenDebug(String s) {
        return handler.debug(1) ? s : "";
    }

}
