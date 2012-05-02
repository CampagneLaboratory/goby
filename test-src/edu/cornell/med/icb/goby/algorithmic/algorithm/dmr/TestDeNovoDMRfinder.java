/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.shorts.Short2LongLinkedOpenHashMap;
import org.junit.Assert;
import org.junit.Test;

import java.lang.reflect.Array;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: nyasha
 * Date: 4/9/12
 * Time: 2:20 PM
 */
public class TestDeNovoDMRfinder {


    @Test
    public void testNoDMRfound() {


    }

  //  @Test
    public void testOneDMRinMiddleInterval() {

        // Meth counts group1 ->    0  0   0    10   10  10  10  0   0   0
        // Meth counts group2 ->    0  0   0    0   0   0   0   0   0   0
        // Unmeth counts group1->   10 10 10    0  0   0   0   10  10  10
        // Unmeth counts group2 ->  14  15  8   12  9   8   5   12  10  12
        // expected result --> DMR between indices: 3 and 6 inclusive

        final SlidingCountArray group1MethylatedCounts = new SlidingCountArray(10);
        final int[] a = {0, 0, 0, 10, 20, 30, 40, 40, 40, 40};
        group1MethylatedCounts.setCumC(a);

        final SlidingCountArray group2MethylatedCounts = new SlidingCountArray(10);
        int[] b = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        group2MethylatedCounts.setCumC(b);

        final SlidingCountArray group1UnMethylatedCounts = new SlidingCountArray(10);
        int[] c = {10, 20, 30, 30, 30, 30, 30, 40, 50, 60};
        group1UnMethylatedCounts.setCumC(c);

        final SlidingCountArray group2UnMethylatedCounts = new SlidingCountArray(10);
        int[] d = {14, 29, 37, 49, 58, 66, 71, 83, 93, 105};
        group2MethylatedCounts.setCumC(d);

        final SlidingCountArray[] methylatedCountsArray = {group1MethylatedCounts, group2MethylatedCounts};
        final ObjectArrayList<SlidingCountArray> methylatedCounts = new ObjectArrayList(methylatedCountsArray);
        final SlidingCountArray[] unMethylatedCountsArray = {group1UnMethylatedCounts, group2UnMethylatedCounts};
        final ObjectArrayList<SlidingCountArray> unMethylatedCounts = new ObjectArrayList(unMethylatedCountsArray);

        final GroupComparison groupComp = new GroupComparison("Group1", "Group2", 0, 1, 0);

        DeNovoDMRfinder finder = new DeNovoDMRfinder(0.05, methylatedCounts, unMethylatedCounts);

        assertEquals("The number of DMRs found: ", 1, finder.search(groupComp).size());
        assertEquals("DMR coordinates: start is ", 3, finder.dmrResultList.get(0).start);
        assertEquals("DMR coordinates: end is ", 6, finder.dmrResultList.get(0).end);

    }

   // @Test
    public void testOneDMRinMiddleofMiddleInterval() {

        // expected result --> DMR between indices: 7 and 10 inclusive

        final SlidingCountArray group1MethylatedCounts = new SlidingCountArray(18);
        final int[] a = {0, 0, 0, 0, 0, 0, 0, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40};
        group1MethylatedCounts.setCumC(a);

        final SlidingCountArray group2MethylatedCounts = new SlidingCountArray(18);
        int[] b = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        group2MethylatedCounts.setCumC(b);

        final SlidingCountArray group1UnMethylatedCounts = new SlidingCountArray(18);
        int[] c = {10, 20, 30, 40, 50, 60, 70, 70, 70, 70, 70, 80, 90, 100, 110, 120, 130, 140};
        group1UnMethylatedCounts.setCumC(c);

        final SlidingCountArray group2UnMethylatedCounts = new SlidingCountArray(18);
        int[] d = {14, 29, 37, 49, 58, 73, 81, 93, 102, 110, 115, 127, 137, 149, 165, 221, 246, 267};
        group2MethylatedCounts.setCumC(d);

        final SlidingCountArray[] methylatedCountsArray = {group1MethylatedCounts, group2MethylatedCounts};
        final ObjectArrayList<SlidingCountArray> MethylatedCounts = new ObjectArrayList(methylatedCountsArray);
        final SlidingCountArray[] UnMethylatedCountsArray = {group1UnMethylatedCounts, group2UnMethylatedCounts};
        final ObjectArrayList<SlidingCountArray> UnMethylatedCounts = new ObjectArrayList(UnMethylatedCountsArray);

        final GroupComparison groupComp = new GroupComparison("Group1", "Group2", 0, 1, 0);

        DeNovoDMRfinder finder = new DeNovoDMRfinder(0.05, MethylatedCounts, UnMethylatedCounts);

        assertEquals("The number of DMRs found: ", 1, finder.search(groupComp).size());
        assertEquals("DMR coordinates: start is ", 7, finder.dmrResultList.get(0).start);
        assertEquals("DMR coordinates: end is ", 10, finder.dmrResultList.get(0).end);

    }

   // @Test
    public void testTwoDMRsinMiddleInterval() {
        // expected result --> DMR between indices: 11 and 14 inclusive
        //                 --> DMR between indices: 16 and 19 inclusive

        final SlidingCountArray group1MethylatedCounts = new SlidingCountArray(30);
        final int[] a = {0,0,0,0,0,0,0,0,0,0,0,56,71,127,175,175,205,261,289,328,328,328,328,328,328,328,328,328,328};
        group1MethylatedCounts.setCumC(a);

        final SlidingCountArray group2MethylatedCounts = new SlidingCountArray(30);
        int[] b = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        group2MethylatedCounts.setCumC(b);

        final SlidingCountArray group1UnMethylatedCounts = new SlidingCountArray(30);
        int[] c = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,230,253,263,278,398,410,422,440};
        group1UnMethylatedCounts.setCumC(c);

        final SlidingCountArray group2UnMethylatedCounts = new SlidingCountArray(30);
        int[] d = {14,29,37,49,58,73,81,93,102,110,115,127,137,149,165,221,246,267,323,349,375,387,402,414,459,474,489,490,511};
        group2MethylatedCounts.setCumC(d);

        final SlidingCountArray[] methylatedCountsArray = {group1MethylatedCounts, group2MethylatedCounts};
        final ObjectArrayList<SlidingCountArray> MethylatedCounts = new ObjectArrayList(methylatedCountsArray);
        final SlidingCountArray[] UnMethylatedCountsArray = {group1UnMethylatedCounts, group2UnMethylatedCounts};
        final ObjectArrayList<SlidingCountArray> UnMethylatedCounts = new ObjectArrayList(UnMethylatedCountsArray);

        final GroupComparison groupComp = new GroupComparison("Group1", "Group2", 0, 1, 0);

            DeNovoDMRfinder finder = new DeNovoDMRfinder(0.05, MethylatedCounts, UnMethylatedCounts);

        assertEquals("The number of DMRs found: ", 2, finder.search(groupComp).size());
        assertEquals("DMR coordinates: start is ", 11, finder.dmrResultList.get(0).start);
        assertEquals("DMR coordinates: end is ", 14, finder.dmrResultList.get(0).end);
        assertEquals("DMR coordinates: start is ", 16, finder.dmrResultList.get(1).start);
        assertEquals("DMR coordinates: end is ", 19, finder.dmrResultList.get(1).end);


    }

}
