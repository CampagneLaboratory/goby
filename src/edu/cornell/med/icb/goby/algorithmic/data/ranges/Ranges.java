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

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

/**
 * Stores list of ranges for a set of chromosomes.
 *
 * @author Fabien Campagne
 *         Date: 10/19/11
 *         Time: 11:35 AM
 */
public class Ranges {
    private Int2ObjectMap<ObjectArrayList<Range>> map;

    public Ranges() {
        this.map = new Int2ObjectOpenHashMap<ObjectArrayList<Range>>();
    }

    public void add(Range range, int referenceIndex) {
        ObjectArrayList<Range> list = map.get(referenceIndex);
        if (list == null) {
            list = new ObjectArrayList<Range>();
            map.put(referenceIndex, list);
        }

        list.add(range);

    }

    public void order() {
        for (ObjectArrayList<Range> value : map.values()) {
            Collections.sort(value);
        }
    }

    public Range findNonOverlappingRange(int referenceIndex, int position) {
        ObjectArrayList<Range> result = new ObjectArrayList<Range>();
        ObjectArrayList<Range> list = map.get(referenceIndex);
      final Range query = new Range();
        query.min = position;
        query.max = position;

        if (list==null || list.isEmpty()) {
            // no annotations, so query does not overlap..
            return query;
        }

        final int index = Collections.binarySearch(list, query);
        Range nonOverlaping = new Range();
        if (position == 0) {
            return new Range(0, 0);
        }
        int ip;
        Range closest;
        if (index < 0) {

            // position does not match min exactly, we get the insertion point instead:
            ip = -index - 1;
            closest = findClosest(list, position, ip);
        } else {
            ip = index;
            closest=list.get(ip);
        }

        if (!overlaps(position, closest)) {
            return query;
        } else {
            // position overlaps with the closest range. Find the next closest, non-overlapping range to define the
            // non overlapping interval closest to position:
            Range a = findNonOverlapping(list, position, ip, -1);
            Range b = findNonOverlapping(list, position, ip, +1);
            double distancePositionToA = Math.sqrt(Math.pow(a.max - closest.min, 2));
            double distancePositionToB = Math.sqrt(Math.pow(b.min - closest.max, 2));
            if (distancePositionToA > distancePositionToB) {
                // B is closest to position:
                nonOverlaping.min = closest.max;
                nonOverlaping.max = b.min;
            } else { // A is closest to position:
                nonOverlaping.min = a.max;
                nonOverlaping.max = closest.min;
            }
            return nonOverlaping;
        }


    }

    /**
     * Returns true if position overlaps range, completely or partially.
     *
     * @param position
     * @param range
     * @return
     */
    private static boolean overlaps(int position, Range range) {
        return (position >= range.min && position <= range.max);
    }

    protected static Range findNonOverlapping(ObjectArrayList<Range> list, int position, int ip, int offset) {
        boolean overlapping = true;
        Range range = null;
        int j = ip + offset;
        int maxJ = list.size() - 1;

        while (overlapping && j >= 0 && j <= maxJ) {

            range = list.get(j);
            overlapping = overlaps(position, range);
            j += offset;
        }
        if (!overlapping) {
            return range;
        } else {
            if (j < 0) {
                Range first = list.get(0);
                int min = Math.max(0, first.min - 1000);
                return new Range(min, min);
            } else if (j > maxJ) {
                Range last = list.get(maxJ);
                if (!overlaps(position, last)) return last;
                else return new Range(last.max + 1000, last.max + 1000);
            }
            assert false : " should never happen.";
            return null;
        }
    }

    protected static Range findClosest(ObjectArrayList<Range> list, int position, int ip) {
        Range a = findClosest(list, position, ip, -1);
        Range b = findClosest(list, position, ip, +1);
        double distancePositionToA = Math.sqrt(Math.pow(a.max - position, 2));
        double distancePositionToB = Math.sqrt(Math.pow(b.min - position, 2));
        if (distancePositionToA > distancePositionToB) {
            return b;
        } else {
            return a;
        }
    }

    protected static Range findClosest(ObjectArrayList<Range> list, int position, int ip, int offset) {
        boolean found = false;
        Range range = null;
        int j = ip+offset;
        int maxJ = list.size() - 1;

        if (j >= 0 && j <= maxJ) {

            range = list.get(j);
            found = true;
        }
        if (found) {
            return range;
        } else {
            if (j < 0) {
                Range first = list.get(0);
                int min = Math.max(0, first.min - 1000);
                return new Range(min, min);
            } else if (j > maxJ) {
                Range last = list.get(maxJ);
                if (!overlaps(position, last)) return last;
                else return new Range(last.max + 1000, last.max + 1000);
            }
            assert false : " should never happen.";
            return null;
        }
    }

    public ObjectArrayList<Range> forReference(int targetIndex) {
        return map.get(targetIndex);
    }
}
