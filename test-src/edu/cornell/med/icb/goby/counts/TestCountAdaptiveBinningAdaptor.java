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

package edu.cornell.med.icb.goby.counts;

import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * @author Fabien Campagne
 *         Date: 6/11/11
 *         Time: 12:48 PM
 */
public class TestCountAdaptiveBinningAdaptor {
    CountsReaderTestSupport source1;
    private CountBinningAdapterI binAdapter;
    private CountsReaderTestSupport source2;


    @Test
    public void testSource1() throws IOException {
        source1 = new CountsReaderTestSupport("(10,0)(2,2)(10,0)(4,6)");  //  2*2+4*6/(4+2)=28/6=4.67 is average count over entire range.
        binAdapter = new CountAdaptiveBinningAdaptor(source1);   // bin size overlaps entire source range.

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(0, binAdapter.getPosition());
        assertEquals(10, binAdapter.getLength());
        assertEquals(0, binAdapter.getAverage(), .1);
        assertEquals(0, binAdapter.getMax(), .1);

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(10, binAdapter.getPosition());
        assertEquals(2, binAdapter.getLength());
        assertEquals(2, binAdapter.getAverage(), .1);
        assertEquals(2, binAdapter.getMax(), .1);

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(12, binAdapter.getPosition());
        assertEquals(10, binAdapter.getLength());
        assertEquals(0, binAdapter.getAverage(), .1);
        assertEquals(0, binAdapter.getMax(), .1);

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(22, binAdapter.getPosition());
        assertEquals(4, binAdapter.getLength());
        assertEquals(6, binAdapter.getAverage(), .1);
        assertEquals(6, binAdapter.getMax(), .1);

        assertFalse(binAdapter.hasNextTransition());

    }

    @Test
    public void testSource2() throws IOException {
        source2 = new CountsReaderTestSupport("(10,0)(2,2)(1,3)(4,1)(10,0)(4,6)");
        binAdapter = new CountAdaptiveBinningAdaptor(source2);   // bin size overlaps entire source range.

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(0, binAdapter.getPosition());
        assertEquals(10, binAdapter.getLength());
        assertEquals(0, binAdapter.getAverage(), .1);
        assertEquals(0, binAdapter.getMax(), .1);

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(10, binAdapter.getPosition());
        assertEquals(7, binAdapter.getLength());
        assertEquals((2d*2+1*3+4*1  )/(2d+1+4), binAdapter.getAverage(), .1);
        assertEquals(3.0, binAdapter.getMax(), .1);

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(17, binAdapter.getPosition());
        assertEquals(10, binAdapter.getLength());
        assertEquals(0, binAdapter.getAverage(), .1);
        assertEquals(0, binAdapter.getMax(), .1);

        assertTrue(binAdapter.hasNextTransition());
        binAdapter.nextTransition();
        assertEquals(27    , binAdapter.getPosition());
        assertEquals(4, binAdapter.getLength());
        assertEquals(6, binAdapter.getAverage(), .1);
        assertEquals(6, binAdapter.getMax(), .1);

        assertFalse(binAdapter.hasNextTransition());

    }

}
