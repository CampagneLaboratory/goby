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
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.IOException;

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
    private int previousPosition;
    private int previousTargetIndex;
    int index = 0;
    private final String label;
    private final String HAS_LINK_LABEL;
    private final String POS_LABEL;
    private final String TARGETS_LABEL;
    private final String FRAGMENTS_LABEL;
    private final AlignmentCollectionHandler handler;
    private int fragmentIndex;


    LinkInfo(final AlignmentCollectionHandler handler, final String label) {
        this.handler = handler;
        this.label = label;
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
        deltaPositionIndex=0;
        fragmentIndex=0;
    }

    public Alignments.RelatedAlignmentEntry code(final boolean hasLink, final Alignments.RelatedAlignmentEntry relatedLink) {
        hasLinks.add(hasLink ? 1 : 0);
        if (!hasLink) {
            return null;
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
        if (!justResetPos) {

            deltaPositions.add(Fast.int2nat(position - previousPosition));
            deltaTargetIndices.add(Fast.int2nat(targetIndex - previousTargetIndex));
        }

        previousPosition = position;
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
    int deltaPositionIndex;
    public Alignments.RelatedAlignmentEntry decode(int originalIndex, final Alignments.RelatedAlignmentEntry source) {
        if (hasLinks.getInt(originalIndex) == 1) {
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

                position = Fast.nat2int(deltaPositions.getInt(deltaPositionIndex)) + previousPosition;
                targetIndex = Fast.nat2int(deltaTargetIndices.getInt(deltaPositionIndex)) + previousTargetIndex;
                deltaPositionIndex++;
                //System.out.println("getting pos from deltas, position="+position);

            }
            result.setPosition(position);
            result.setTargetIndex(targetIndex);
            previousPosition = position;
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
