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

import java.io.IOException;

/**
 * An interface for implementations that replace query indices with small values.
 *
 * @author Fabien Campagne
 *         Date: 3/5/12
 *         Time: 5:31 PM
 */
public interface QueryIndexPermutationInterface {

    Alignments.AlignmentEntry makeSmallIndices(Alignments.AlignmentEntry entry);

    void makeSmallIndices(Alignments.AlignmentEntry.Builder entry);

    int getSmallestIndex();

    int getBiggestSmallIndex();

    void setSmallestIndex(int value);

    void setBiggestSmallIndex(int value);

    /**
     * Permutate a query index and return the smaller value.
     * @param index
     * @return
     */
    int permutate(int index);

    /**
     * Query indices will pruned from the permutation map after they have been requested x times. Pruned indices are
     * written to disk.
     * @param x
     */
    public void setPruneLimit(byte x);

    /**
     * Flush all information to disk as needed and release resources.
     */
    public void close();
}
