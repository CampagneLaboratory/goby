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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import edu.cornell.med.icb.goby.alignments.processors.ObservedIndel;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.lang.MutableString;

/**
 * Implements the method of Krawitz et al to determine a span of equivalent indel regions, given
 * an observed indel reported by an aligner.
 * Citation: Peter Krawitz et al Bioinformatics 2011. Microindel detection in short-read data.
 *
 * @author Fabien Campagne
 *         Date: 6/7/11
 *         Time: 6:00 PM
 */
public class EquivalentIndelRegionCalculator {
    public static final int FLANK_SIZE = 4;
    RandomAccessSequenceInterface genome;

    public EquivalentIndelRegionCalculator(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    MutableString p = new MutableString();

    /**
     * Determine the span of equivalent indel regions for the observed indel. Indel are often non uniquely described
     * by the information provided by observed indel.
     *
     * @param indel the observed indel
     * @return the span of equivalent indel regions
     */
    public EquivalentIndelRegion determine(final int referenceIndex, final ObservedIndel indel) {
        final EquivalentIndelRegion result = new EquivalentIndelRegion();
        result.startPosition = indel.getStart();
        result.endPosition = indel.getEnd();

        int ip = indel.getStart();
        p.setLength(0);
        boolean insertion = insertionInRead(indel);
        p.append(insertion ? indel.to() : indel.from());

        // extend left:
        int leftExtensions = 0;
        final int indelSize = p.length();
        final int lastBaseIndex = indelSize - 1;


        while (p.charAt(lastBaseIndex + leftExtensions % indelSize) == genome.get(referenceIndex, result.startPosition)) {
            leftExtensions++;
            result.startPosition = indel.getStart() - leftExtensions;
        }
        int rightExtensions = 0;

        while (p.charAt(rightExtensions % indelSize) == genome.get(referenceIndex, result.endPosition)) {
            rightExtensions++;
            result.endPosition = indel.getEnd() + rightExtensions;
        }
        from.setLength(0);
        to.setLength(0);
        final MutableString toFill;
        final MutableString gaps;

        if (insertion) {
            toFill = to;
            gaps = from;
        } else {
            toFill = from;
            gaps = to;
        }
        for (int i = result.startPosition + 1; i < result.startPosition + leftExtensions; ++i) {
            char c = genome.get(referenceIndex, i);
            toFill.append(c);
            gaps.append(c);
        }
        toFill.append(insertion ? indel.to() : indel.from());
        for (int i = 0; i < indelSize; i++) {
            gaps.append('-');

        }
        for (int i = result.endPosition - rightExtensions; i < result.endPosition; ++i) {
            final char c = genome.get(referenceIndex, i);
            toFill.append(c);
            gaps.append(c);
        }
        if (insertion) {
            result.from = gaps.toString();
            result.to = toFill.toString();
        } else {
            result.from = toFill.toString();
            result.to = gaps.toString();
        }
        flankingLeft.setLength(0);
        for (int i = result.startPosition - FLANK_SIZE + 1; i >= 0 && i < result.startPosition + 1; ++i) {
            flankingLeft.append(genome.get(referenceIndex, i));
        }
        result.flankLeft = flankingLeft.toString();
        flankingRight.setLength(0);

        final int maxRefLength = genome.getLength(referenceIndex);
        for (int i = result.endPosition; i < maxRefLength && i < result.endPosition + FLANK_SIZE; ++i) {
            flankingRight.append(genome.get(referenceIndex, i));
        }
        result.flankRight = flankingRight.toString();
        return result;

    }

    MutableString from = new MutableString();
    MutableString to = new MutableString();
    MutableString flankingLeft = new MutableString();
    MutableString flankingRight = new MutableString();

    private boolean insertionInRead(final ObservedIndel indel) {
        final String from = indel.from();
        final String to = indel.to();
        if (from.indexOf('-') >= 0) {
            return true;
        } else if (to.indexOf('-') >= 0) {
            return false;
        }
        throw new InternalError("indel must either be an insertionInRead or a deletion");
    }

}
