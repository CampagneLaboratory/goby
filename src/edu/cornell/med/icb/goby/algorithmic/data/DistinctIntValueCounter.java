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

import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntCollection;

import java.util.BitSet;

/**
 * Counts the number of distinct integers observed. A naive implementation would store n distinct integers
 * in a set and return the size of this set. This does not scale when the number of integers observed grows very large.
 * This implementation is memory efficient because it stores only missing values. The number of distinct integers can
 * be calculated knowing the minimum value observed, the maximum and the number of missing values in the sequence.
 *
 * @author Fabien Campagne
 *         Date: Mar 15, 2011
 *         Time: 12:35:03 PM
 */
public class DistinctIntValueCounter {
    private class MissingValues {
        BitSet missingValues = new BitSet();
      
        public void add(int value) {
            missingValues.set(value);
        }

        public void remove(int value) {
            missingValues.clear(value);

        }

        public int size() {
            return missingValues.cardinality();
        }
    }

    MissingValues missingValues = new MissingValues();

    //IntSet missingValues = new IntOpenHashSet();
    int minValue = Integer.MAX_VALUE;
    int maxValue = Integer.MIN_VALUE;

    public void observe(IntCollection values) {
        for (int i : values) observe(i);
    }

    public void observe(int[] values) {
        for (int i : values) observe(i);
    }

    public void observe(int value) {
        if (value < maxValue) {
            // consider value could have been missing in previous observations:

        }
        int previousMax = maxValue;
        minValue = Math.min(value, minValue);
        maxValue = Math.max(value, maxValue);
        if (value == maxValue && previousMax != Integer.MIN_VALUE) {
            // value is the new max. Determine which missing values
            for (int i = previousMax + 1; i < value; i++) {
                missingValues.add(i);

                //  missingValues.add(i);
            }
        }
    }

    /**
     * Return the number of distinct integer values in the observed sequence.
     *
     * @return
     */
    public int count() {
        return maxValue - minValue - missingValues.size() + (minValue == 0 ? 1 : 0);
    }
}
