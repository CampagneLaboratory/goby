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
import org.apache.log4j.Logger;


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

    RandomAccessSequenceInterface genome;
    private static final String GAPS = "----------------------------------------------------------------";
    private int flankRightSize = 4;
    public int flankLeftSize = 4;
    private static final Logger LOG = Logger.getLogger(EquivalentIndelRegionCalculator.class);

    /**
     * Set an array with the permutation to convert alignment reference indices to genome reference indices.
     * Each element of the array at index alignmentTargetIndex holds the genomeTargetIndex. One element
     * per reference sequence in the genome and alignment.
     *
     * @param alignmentToGenomeTargetIndices permutation.
     */
    public void setReferenceIndexPermutation(final int[] alignmentToGenomeTargetIndices) {
        this.alignmentToGenomeTargetIndices = alignmentToGenomeTargetIndices;
    }

    private int[] alignmentToGenomeTargetIndices;

    /**
     * Set the number of bases to record in flank-right.
     *
     * @param flankRightSize number of flanking bases on the right.
     */
    public void setFlankRightSize(final int flankRightSize) {
        this.flankRightSize = flankRightSize;
    }

    /**
     * Set the number of bases to record in flank-left.
     *
     * @param flankLeftSize number of flanking bases on the left.
     */
    public void setFlankLeftSize(final int flankLeftSize) {

        this.flankLeftSize = flankLeftSize;
    }

    /**
     * Return the number of bases flanking this indel to hte right.
     *
     * @return flank right size.
     */
    public int getFlankRightSize() {
        return flankRightSize;
    }

    /**
     * Return the number of bases flanking this indel to the left.
     *
     * @return flank left size.
     */
    public int getFlankLeftSize() {
        return flankLeftSize;
    }

    public EquivalentIndelRegionCalculator(
            final RandomAccessSequenceInterface genome,
            final int[] alignmentToGenomeTargetIndices) {
        init(genome, alignmentToGenomeTargetIndices);

    }

    public EquivalentIndelRegionCalculator(final RandomAccessSequenceInterface genome) {
        final int[] permutation = new int[genome.size()];
        for (int i = 0; i < genome.size(); i++) {
            permutation[i] = i;
        }
        init(genome, permutation);

    }

    private void init(RandomAccessSequenceInterface genome, int[] alignmentToGenomeTargetIndices) {
        this.genome = genome;
        if (genome == null) {
            throw new IllegalArgumentException("genome cannot be null");
        }
        if (alignmentToGenomeTargetIndices == null) {
            throw new IllegalArgumentException("permutation cannot be null");
        }

        this.alignmentToGenomeTargetIndices = alignmentToGenomeTargetIndices;
    }

    MutableString p = new MutableString();

    /**
     * Determine the span of equivalent indel regions for the observed indel. Indels are often non uniquely described
     * by the information provided by observed indel.
     *
     * @param referenceIndex Index of the reference sequence where the indel is located.
     * @param indel          the observed indel
     * @return the span of equivalent indel regions , or null, if the indel positions are outside the boundaries of the genome sequence.
     */
    public EquivalentIndelRegion determine(final int referenceIndex, final ObservedIndel indel) {
        final EquivalentIndelRegion result = new EquivalentIndelRegion();
        result.startPosition = indel.getStart();
        result.endPosition = indel.getEnd();
        result.referenceIndex = referenceIndex;
        result.readIndices.add( indel.readIndex);
        p.setLength(0);
        final boolean insertion = insertionInRead(indel);
        p.append(insertion ? indel.to() : indel.from());
        final int genomeReferenceIndex = alignmentToGenomeTargetIndices[referenceIndex];
        // extend left:
        int leftExtensions = 0;
        final int indelSize = p.length();
        final int lastBaseIndex = indelSize - 1;
        int rewindLeft = 0;
        if (result.startPosition > genome.getLength(genomeReferenceIndex)) {
            LOG.warn(String.format("Cannot determine sequence at position %d of reference-index %d ",
                    result.startPosition,
                    referenceIndex));
            return null;
        }

        while (result.startPosition >= 1 &&
                p.charAt(lastBaseIndex - leftExtensions + rewindLeft) == genome.get(genomeReferenceIndex, result.startPosition)) {
            leftExtensions++;
            result.startPosition = indel.getStart() - leftExtensions;
            if (lastBaseIndex - leftExtensions + rewindLeft < 0) {
                rewindLeft += indelSize;
            }
            //       debug("extending left", result);
        }
        int rightExtensions = 0;
        int rewindRight = 0;
        int chromosomeSize = genome.getLength(genomeReferenceIndex);
        while (result.endPosition < chromosomeSize &&
                p.charAt(rightExtensions + rewindRight) == genome.get(genomeReferenceIndex, result.endPosition)) {
            rightExtensions++;
            result.endPosition = indel.getEnd() + rightExtensions;
            if (rightExtensions + rewindRight >= indelSize) {
                rewindRight -= indelSize;
            }
            //      debug("extending right:", result);
        }


        from.setLength(0);
        to.setLength(0);
        toFill.setLength(0);

        genome.getRange(genomeReferenceIndex, result.startPosition + 1, result.endPosition - result.startPosition - 1, from);

        if (insertion)

        {

            // construct the read sequence in the insertion region of the eir:
            to.append(roll(leftExtensions, indel.to()));
            to.append(from);

            from.insert(0, GAPS.subSequence(0, indelSize));

        } else

        {
            // construct the read sequence in the deletion region of the eir:
            final int length = from.length();
            to.append(GAPS.subSequence(0, indelSize));
            to.append(from.subSequence(Math.min(indelSize, length), length));

        }

        result.from = from.toString();
        result.to = to.toString();
        /*  gaps.setLength(0);
        gaps.append(toFill);
        gaps.delete(0, indelSize);
        gaps.insert(0, GAPS.subSequence(0, indelSize));
        */

        //     debug("from: ", result);
        flankingLeft.setLength(0);
        genome.getRange(genomeReferenceIndex, result.startPosition - flankLeftSize + 1, flankLeftSize, flankingLeft);

        result.flankLeft = flankingLeft.toString();

        flankingRight.setLength(0);
        genome.getRange(genomeReferenceIndex, result.endPosition, flankRightSize, flankingRight);

        final int maxRefLength = genome.getLength(genomeReferenceIndex);

        result.flankRight = flankingRight.toString();
        //              debug("flanks: ", result);
        return result;

    }

    final MutableString rollBuffer = new MutableString();

    private final MutableString roll(int leftExtensions, String from) {
        rollBuffer.setLength(0);
        rollBuffer.append(from);
        final int length = from.length();
        final int lastCharIndex = length - 1;
        for (int i = 0; i < leftExtensions; i++) {

            rollBuffer.insert(0, rollBuffer.charAt(lastCharIndex));
            rollBuffer.setLength(length);
        }
        return rollBuffer;
    }

    final MutableString toFill = new MutableString();
    final MutableString gaps = new MutableString();

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
