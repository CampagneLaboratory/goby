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

//import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import edu.cornell.med.icb.goby.alignments.processors.ObservedIndel;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 6/7/11
 *         Time: 6:11 PM
 */
public class TestEquivalentIndelRegionCalculator {
    RandomAccessSequenceTestSupport genome;
    private String[] sequences = {
            "ACTCAAAGACT",  // will delete one A in the three consecutive As
            "AAACAGATCCCACA"  // will insert AG between C and AG
    };
    private EquivalentIndelRegionCalculator equivalentIndelRegionCalculator;

    @Test
    public void emptyTest() {
    }

    @Before
    public void setUp() throws Exception {

        genome = new RandomAccessSequenceTestSupport(sequences);
        equivalentIndelRegionCalculator = new EquivalentIndelRegionCalculator(genome);
    }

    @Test
    public void testDetermine() throws Exception {
        ObservedIndel indel = new ObservedIndel(4, 5, "-", "A");
        EquivalentIndelRegion result = equivalentIndelRegionCalculator.determine(0, indel);
        assertEquals(0, result.referenceIndex);
        assertEquals(3, result.startPosition);
        assertEquals(7, result.endPosition);
        assertEquals("-AA", result.from);
        assertEquals("AAA", result.to);
        assertEquals("ACTC", result.flankLeft);
        assertEquals("GACT", result.flankRight);
    }

    @Test
    public void testInsertAG() throws Exception {
        ObservedIndel indel = new ObservedIndel(3, 4, "--", "AG");
        EquivalentIndelRegion result = equivalentIndelRegionCalculator.determine(1, indel);
        assertEquals(0, result.referenceIndex);
        assertEquals(3, result.startPosition);
        assertEquals(7, result.endPosition);
        assertEquals("--AGA", result.from);
        assertEquals("AGAGA", result.to);
        assertEquals("AAAC", result.flankLeft);
        assertEquals("TCCC", result.flankRight);
        assertEquals("AAAC--AGATCCC", result.fromInContext());
        assertEquals("AAACAGAGATCCC", result.toInContext());
        // "AAAC  AGATCCC"
        // "AAACAGAGATCCC"
    }


}
