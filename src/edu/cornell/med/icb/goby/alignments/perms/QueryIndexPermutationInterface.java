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
 * An interface for implementations that replace query indices with small values. Some of these implementations will
 * keep a record of all associations made, others will just discard them.
 *
 * @author Fabien Campagne
 *         Date: 3/5/12
 *         Time: 5:31 PM
 */
public interface QueryIndexPermutationInterface {
    /**
     * Give the entry a small index according to the entry query occurrence limit, or to the global prune limit.
     *
     * @param entry alignment entry.
     * @return entry with permuted query_index value.
     */
    Alignments.AlignmentEntry makeSmallIndices(Alignments.AlignmentEntry entry);

    /**
     * Give the entry a small index according to the entry query occurrence limit, or to the global prune limit.
     * The query index will be changed in place.
     *
     * @param entry alignment entry builder.
     */
    void makeSmallIndices(Alignments.AlignmentEntry.Builder entry);

    /**
     * Get the smallest value of small index returned so far.
     *
     * @return smallest value of small index returned so far.
     */
    int getSmallestIndex();

    /**
     * Get the largest value of small index returned so far.
     *
     * @return largest value of small index returned so far.
     */
    int getBiggestSmallIndex();

    /**
     * Set the smallest value of small index returned so far.
     *
     * @param value value for set.
     */
    void setSmallestIndex(int value);

    /**
     * Set the largest value of small index returned so far.
     * @param value value for set.
     */
    void setBiggestSmallIndex(int value);

    /**
     * Permutate a query index and return the smaller value. The global prune limit is used as
     * maxQueryIndexOccurrence for each index.
     *
     * @param queryIndex the large randomly distributed query index to permutate to a monotonically increasing value.
     * @return the small index associated with queryIndex, or -1 is this query index has been requested more than the global prune limit.
     */
    int permutate(int queryIndex);

    /**
     * Permutate a query index and return the smaller value.
     *
     * @param queryIndex
     * @param maxQueryIndexOccurrence the maximum number of times the queryIndex can be requested before it is pruned and written to disk.
     * @return the small index associated with queryIndex.
     */
    int permutate(int queryIndex, int maxQueryIndexOccurrence);

    /**
     * Query indices will pruned from the permutation map after they have been requested x times. Pruned indices are
     * written to disk.
     *
     * @param x query index occurrence limit.
     */
    public void setPruneLimit(byte x);

    /**
     * Flush all information to disk as needed and release resources.
     */
    public void close();
}
