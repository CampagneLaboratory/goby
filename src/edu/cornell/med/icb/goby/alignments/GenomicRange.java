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

import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;


/**
 * Defines a range/slice of a genome.
 *
 * @author Fabien Campagne
 *         Date: 10/19/11
 *         Time: 12:44 PM
 */
public class GenomicRange {
    public int startReferenceIndex;
    public String startChromosome;
    public int startPosition;
    public int endReferenceIndex;
    public String endChromosome;
    public int endPosition;
    private DoubleIndexedIdentifier ids;

    public void resolveChromosomeIndices(DoubleIndexedIdentifier referenceIds) {
        startReferenceIndex = referenceIds.getIndex(startChromosome);
        endReferenceIndex = referenceIds.getIndex(endChromosome);
    }

    /**
     *
     * @param chromosome   chromosome.
     * @param segmentStart zero-based position of the start of the segment on chromosome.
     * @param segmentEnd   zero-based position of the end of the segment on chromosome.
     * @return
     */
    public boolean fullyContains(String chromosome, int segmentStart, int segmentEnd) {
        assert ids != null : " identifiers must have been provided. Call setTargetIds to provide. ";
        int chromosomeIndex = ids.getIndex(chromosome);
        if (chromosomeIndex == -1) return false;
        if (chromosomeIndex < startReferenceIndex) return false;
        if (chromosomeIndex > endReferenceIndex) return false;
        if (chromosomeIndex == startReferenceIndex && segmentStart < startPosition) {
            return false;
        }
        if (chromosomeIndex == endReferenceIndex && segmentEnd > endPosition) {
            return false;
        }
        return true ;


    }


    public void setTargetIds(IndexedIdentifier targetIdentifiers) {
        ids = new DoubleIndexedIdentifier(targetIdentifiers);
        startReferenceIndex = ids.getIndex(startChromosome);
        endReferenceIndex = ids.getIndex(endChromosome);
    }
}
