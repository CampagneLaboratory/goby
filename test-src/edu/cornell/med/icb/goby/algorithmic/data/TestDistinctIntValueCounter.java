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

package edu.cornell.med.icb.goby.algorithmic.data;

import it.unimi.dsi.bits.LongArrayBitVector;
import it.unimi.dsi.fastutil.longs.LongSortedSet;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.util.BitSet;
import java.util.Random;

/**
 * @author Fabien Campagne
 *         Date: Mar 15, 2011
 *         Time: 12:43:45 PM
 */
public class TestDistinctIntValueCounter {
    @Test
    public void testCount1() {
        DistinctIntValueCounterBitSet counter = new DistinctIntValueCounterBitSet();
        counter.observe(0);
        counter.observe(1);
        counter.observe(2);
        counter.observe(3);
        assertEquals(4, counter.count());
        counter.observe(4);
        assertEquals(5, counter.count());
        // value 5 is missing

        counter.observe(6);
        assertEquals(6, counter.count());
        // missing values: 7, 8, 9
        counter.observe(10);
        assertEquals(7, counter.count());
    }

    @Test
    public void testCountNoZero() {
        DistinctIntValueCounterBitSet counter = new DistinctIntValueCounterBitSet();

        counter.observe(1);
        counter.observe(2);
        counter.observe(3);
        assertEquals(3, counter.count());
        counter.observe(4);
        assertEquals(4, counter.count());
        // value 5 is missing

        counter.observe(6);
        assertEquals(5, counter.count());
        // missing values: 7, 8, 9
        counter.observe(10);
        assertEquals(6, counter.count());
    }

    @Test
    public void testSetLargeIndex() {
        DistinctIntValueCounterBitSet counter = new DistinctIntValueCounterBitSet();

        counter.observe(10);
        counter.observe(2129225);
        counter.observe(3129225);
        // we need to be able to store up to Integer.MAX_VALUE query indices:
        counter.observe(Integer.MAX_VALUE);

        assertEquals(4, counter.count());

    }
}



