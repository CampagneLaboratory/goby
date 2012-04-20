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
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.longs.LongSortedSet;

import java.util.BitSet;

/**
 * Counts the number of distinct integers observed. A naive implementation would store n distinct integers
 * in a set and return the size of this set. This does not scale when the number of integers observed grows very large.
 * This implementation uses a BitSet and requires one bit per query index.
 *
 * @author Fabien Campagne
 *         Date: Mar 15, 2011
 *         Time: 12:35:03 PM
 */
public class DistinctIntValueCounterBitSet implements DistinctIntValueCounterInterface {

  //  BitSet values = new BitSet();
    private LongArrayBitVector bitVector=LongArrayBitVector.ofLength(1L<<(32-6));
    private LongSortedSet set=bitVector.asLongSet();
    public void observe(IntCollection values) {
        for (int i : values) observe(i);
    }

    public final void observe(final int[] values) {

        for (int i : values) observe(i);
    }

    public final void observe(final int value) {
        set.add(value);
     }

        /**
     * Return the number of distinct integer values in the observed sequence.
     *
     * @return
     */
    public int count() {

        return set.size();
    }
}
