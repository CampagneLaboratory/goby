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

package edu.cornell.med.icb.goby.algorithmic.data.ranges;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.Test;

import java.util.Collections;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 10/19/11
 *         Time: 11:57 AM
 */
public class RangesTest {
    @Test
    public void testFindNonOverlapping() {

        Range a = new Range();
        a.min = 1;
        a.max = 5;
        Range b = new Range();
        b.min = 8;
        b.max = 10;
        Range c = new Range();
        c.min = 15;
        c.max = 25;

        ObjectArrayList<Range> list = new ObjectArrayList<Range>();
        list.add(a);
        list.add(b);
        list.add(c);
        Collections.sort(list);
        Range resultBefore = Ranges.findNonOverlapping(list, 9, 1, -1);

        assertEquals(new Range(1, 5), resultBefore);
        Range resultAfter = Ranges.findNonOverlapping(list, 9, 1, 1);

        assertEquals(new Range(15, 25), resultAfter);

        Range resultBeforeFirst = Ranges.findNonOverlapping(list, 0, 0, -1);
        assertEquals(new Range(0, 0), resultBeforeFirst);
    }

    @Test
    public void testOverlap() {

        Range a = new Range();
        a.min = 1;
        a.max = 5;
        Range b = new Range();
        b.min = 8;
        b.max = 10;
        Range c = new Range();
        c.min = 15;
        c.max = 25;

        Ranges ranges = new Ranges();
        ranges.add(a, 0);
        ranges.add(b, 0);
        ranges.add(c, 0);
        ranges.order();

        assertEquals(new Range(6, 6), ranges.findNonOverlappingRange(0, 6));

        assertEquals(new Range(5, 8), ranges.findNonOverlappingRange(0, 9));
        assertEquals(new Range(5, 8), ranges.findNonOverlappingRange(0, 10));
        assertEquals(new Range(5, 8), ranges.findNonOverlappingRange(0, 8));
    }
}
