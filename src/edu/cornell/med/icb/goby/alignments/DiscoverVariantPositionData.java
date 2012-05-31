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

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;

/**
 * Stores information collected about each genomic position inspected by IterateSortedAlignmentsImpl (used by
 * DiscoverSequenceVariantsMode).
 *
 * @author Fabien Campagne
 *         Date: 6/6/11
 *         Time: 3:27 PM
 */
public class DiscoverVariantPositionData extends ObjectArrayList<PositionBaseInfo> {
    private static final long serialVersionUID = 9212001398502402859L;
    private ObjectArraySet<EquivalentIndelRegion> candidateIndels;
    private int position;
    private ObjectArraySet<EquivalentIndelRegion> failedIndels;
    private static final ObjectArraySet<EquivalentIndelRegion> EMPTY_SET = new ObjectArraySet<EquivalentIndelRegion>();

    public int getZeroBasedPosition() {
        return position;
    }

    public DiscoverVariantPositionData() {
        super();
        position = -1;
    }

    public DiscoverVariantPositionData(final int position) {
        super();
        this.position = position;
    }

    /**
     * This method is called if a candidate indel is observed whose start position overlaps with position.
     *
     * @param candidateIndel the candidate indels observed with the same startPosition == position
     */
    public void observeCandidateIndel(final EquivalentIndelRegion candidateIndel) {
        if (candidateIndels == null) {
            candidateIndels = new ObjectArraySet<EquivalentIndelRegion>();
        }
        if (!candidateIndels.contains(candidateIndel)) {
            candidateIndels.add(candidateIndel);
        } else {
            for (final EquivalentIndelRegion eir : candidateIndels) {
                if (eir.equals(candidateIndel)) {
                    eir.incrementFrequency();
                }
            }
        }
    }

    public ObjectArraySet<EquivalentIndelRegion> getIndels() {
      // switch to disable indels calls
        // return EMPTY_SET;
       return candidateIndels;
    }

    /**
     * Mark an indel observation as failing genotype filters.
     *
     * @param indel the candidate indel that failed tests.
     */
    public void failIndel(final EquivalentIndelRegion indel) {
        if (candidateIndels!=null) {
            candidateIndels.remove(indel);
        }
        if (failedIndels == null) {
            failedIndels = new ObjectArraySet<EquivalentIndelRegion>();
        }
        failedIndels.add(indel);
    }

    /**
     * Test if this set of genotype observation includes indels.
     *
     * @return True when the genotype observed include indels.
     */
    public boolean hasCandidateIndels() {
        return candidateIndels != null && !candidateIndels.isEmpty();
    }

    public char getReferenceBase() {
        return get(0).from;
    }

    public ObjectArraySet<EquivalentIndelRegion> getFailedIndels() {
        if (failedIndels == null) {
            return EMPTY_SET;
        } else {
            return failedIndels;
        }
    }
}
