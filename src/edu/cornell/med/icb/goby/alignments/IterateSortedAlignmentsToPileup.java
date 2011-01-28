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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.lang.MutableString;

import java.io.PrintWriter;

/**
 * @author Fabien Campagne
 *         Date: Jan 28, 2011
 *         Time: 10:44:15 AM
 */
public class IterateSortedAlignmentsToPileup extends IterateSortedAlignmentsListImpl2 {
    private PrintWriter outWriter;
    private String[] basenameIds;

    private IntArrayList indices = new IntArrayList();
    private ObjectList<Sequence> sequenceBuffers = new ObjectArrayList<Sequence>();
    private IntArrayList alignmentQueryIndices = new IntArrayList();
    private int priorPosition;
    private int startFlapStart;
    // set isPrintIds to facilitate debugging and visual comparisons of sequences in an IDE. 
    private boolean isPrintIds = true;

    private class Sequence {
        String basename;
        int alignmentQueryIndex;
        MutableString bases = new MutableString();
    }

    int maxVariationLength = -1;

    public void processPositions(int referenceIndex, int position, ObjectArrayList<PositionBaseInfo2> positionBaseInfos) {

        for (PositionBaseInfo2 info : positionBaseInfos) {

            int alignmentQueryIndex = info.alignmentEntryQueryIndex;
            int bufferIndex = find(info.readerIndex, alignmentQueryIndex);
            if (info.matchesReference) {
                sequenceBuffers.get(bufferIndex).bases.append(//info.position + ":" +
                        ".");
            } else {
            maxVariationLength=Math.max(info.variationLength, maxVariationLength);
                sequenceBuffers.get(bufferIndex).bases.append(//info.position + ":" +
                        info.to
                        //                + " "
                );
            }
        }
        priorPosition++;
    }

    private int find(int basenameIndex, int alignmentQueryIndex) {
        int bufferIndex = indices.indexOf(alignmentQueryIndex);
        if (bufferIndex == -1) {
            indices.add(alignmentQueryIndex);
            bufferIndex = indices.indexOf(alignmentQueryIndex);
            final Sequence sequence = new Sequence();
            sequence.basename = basenameIds[basenameIndex];
            sequence.alignmentQueryIndex = alignmentQueryIndex;
            for (int i = 0; i < priorPosition; i++) {
                sequence.bases.append('-');
            }
            sequenceBuffers.add(sequence);
            alignmentQueryIndices.add(alignmentQueryIndex);
        }
        return bufferIndex;
    }

    public void initialize(PrintWriter outWriter, String basenameIds[], int flapStartSize) {
        this.outWriter = outWriter;
        this.basenameIds = basenameIds;
        this.startFlapStart = flapStartSize;
        maxVariationLength = -1;
    }

    public void finish() {
        // write alignment to output, reads are grouped by basename.

        for (String basename : basenameIds) {
            int i = 0;
            for (Sequence seq : sequenceBuffers) {
                if (seq.basename.equals(basename)) {
                    final MutableString bases = sequenceBuffers.get(i).bases;
                    if (bases.length() > startFlapStart) {
                        String id = String.format(">%s read %d\n",
                                AlignmentReader.getBasename(seq.basename),
                                seq.alignmentQueryIndex);
                        // we remove the bases from 0 to the end of the flap, just before the longest
                        // variation we have seen in this window. This will always include the variation of
                        // interest in the printed window, and will make sure we show all the bases of the
                        // variations that overlap this window.
                        outWriter.print(String.format("%s%s%n", isPrintIds ? id : "",
                                bases.substring(startFlapStart - maxVariationLength)));
                    }
                }
                i++;
            }
        }
    }
}
