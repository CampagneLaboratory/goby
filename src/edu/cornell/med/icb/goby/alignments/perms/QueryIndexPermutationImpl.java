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

package edu.cornell.med.icb.goby.alignments.perms;

import edu.cornell.med.icb.goby.alignments.Alignments;
import it.unimi.dsi.fastutil.ints.Int2IntAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;

/**
 * Builds a permutation of query indices to small values.
 *
 * @author Fabien Campagne
 *         Date: 3/5/12
 *         Time: 5:10 PM
 */
public class QueryIndexPermutationImpl implements QueryIndexPermutation {

    private int smallestIndex = Integer.MAX_VALUE;
    private int biggestSmallIndex = Integer.MIN_VALUE;

    @Override
    public void reset() {
        smallestIndex = Integer.MAX_VALUE;
        biggestSmallIndex = Integer.MIN_VALUE;
        smallIndexCounter = 0;
        queryIndexPermutation.clear();
        queryIndexPermutation.defaultReturnValue(-1);
    }

    public QueryIndexPermutationImpl() {
        reset();
    }

    @Override
    public Alignments.AlignmentEntry makeSmallIndices(final Alignments.AlignmentEntry entry) {
        final Alignments.AlignmentEntry.Builder merged = Alignments.AlignmentEntry.newBuilder(entry);
        makeSmallIndices(merged);
        return merged.build();
    }

    @Override
    public void makeSmallIndices(final Alignments.AlignmentEntry.Builder entry) {
        final int smallIndex = getSmallIndex(entry.getQueryIndex());
        entry.setQueryIndex(smallIndex);
        smallestIndex = Math.min(smallestIndex, smallIndex);
        biggestSmallIndex = Math.max(biggestSmallIndex, smallIndex);

    }

    @Override
    public int getSmallestIndex() {
        if (smallestIndex == Integer.MAX_VALUE) {
            return 0;
        }
        return smallestIndex;
    }

    @Override
    public int getBiggestSmallIndex() {
        return biggestSmallIndex;
    }

    @Override
    public void setSmallestIndex(int value) {
        smallestIndex = value;
    }

    @Override
    public void setBiggestSmallIndex(int value) {
        biggestSmallIndex = value;
    }

    @Override
    public int permutate(int queryIndex) {
        int smallIndex = getSmallIndex(queryIndex);
        smallestIndex = Math.min(smallestIndex, smallIndex);
        biggestSmallIndex = Math.max(biggestSmallIndex, smallIndex);
        return smallIndex;
    }

    private int smallIndexCounter = 0;
    Int2IntMap queryIndexPermutation = new Int2IntAVLTreeMap();

    private int getSmallIndex(final int queryIndex) {
        final int result = queryIndexPermutation.get(queryIndex);
        // TODO! remove this statement and flush to disk instead:
        if (queryIndexPermutation.size() > 100000) {
            queryIndexPermutation.clear();
        }
        if (result == -1) {
            final int smallIndex = smallIndexCounter++;
            queryIndexPermutation.put(queryIndex, smallIndex);
            return smallIndex;
        } else {
            return result;
        }
    }
}
