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

/**
 * A class that does not permutation query indices to small values.
 *
 * @author Fabien Campagne
 *         Date: 3/5/12
 *         Time: 5:32 PM
 */
public class NoOpPermutation implements QueryIndexPermutationInterface {
    private int smallestIndex = Integer.MAX_VALUE;
    private int biggestSmallIndex = Integer.MIN_VALUE;


    public void reset() {
        smallestIndex = Integer.MAX_VALUE;
        biggestSmallIndex = Integer.MIN_VALUE;


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
        smallestIndex = Math.min(smallestIndex, queryIndex);
        biggestSmallIndex = Math.max(biggestSmallIndex, queryIndex);

        return queryIndex;
    }

    @Override
    public void setPruneLimit(byte x) {

    }

    @Override
    public void close() {

    }

    @Override
    public Alignments.AlignmentEntry makeSmallIndices(final Alignments.AlignmentEntry entry) {
        permutate(entry.getQueryIndex());
        return entry;
    }

    @Override
    public void makeSmallIndices(final Alignments.AlignmentEntry.Builder entry) {
        final int queryIndex = entry.getQueryIndex();
        smallestIndex = Math.min(smallestIndex, queryIndex);
        biggestSmallIndex = Math.max(biggestSmallIndex, queryIndex);

    }

    @Override
    public int getSmallestIndex() {
        return smallestIndex;
    }

    @Override
    public int getBiggestSmallIndex() {
        return biggestSmallIndex;
    }
}
