/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import org.junit.Test;
import static org.junit.Assert.assertTrue;

/**
 * @author Fabien Campagne
 *         Date: Jun 21, 2010
 *         Time: 6:06:45 PM
 */
public class TestAlignmentPositionComparator {
    @Test
    public void compare() {

        AlignmentPositionComparator comparator = new AlignmentPositionComparator();
        Alignments.AlignmentEntry.Builder builderA = newInstance();
        Alignments.AlignmentEntry.Builder builderB = newInstance();
        Alignments.AlignmentEntry.Builder builderC = newInstance();

        //initialized required fields:

        // A < B by target
        // A < C by position
        // B > C by target
        builderA.setTargetIndex(1);
        builderB.setTargetIndex(2);
        builderC.setTargetIndex(1);

        builderA.setPosition(7);
        builderB.setPosition(6);
        builderC.setPosition(8);

        Alignments.AlignmentEntry entryA = builderA.build();
        Alignments.AlignmentEntry entryB = builderB.build();
        Alignments.AlignmentEntry entryC = builderC.build();
        
        assertTrue(comparator.compare(entryA, entryB) < 0);
        assertTrue(comparator.compare(entryA, entryC) < 0);
        assertTrue(comparator.compare(entryB, entryC) > 0);
    }

    private Alignments.AlignmentEntry.Builder newInstance() {
        return Alignments.AlignmentEntry.newBuilder().setQueryIndex(0).setMatchingReverseStrand(false);
    }
}
