/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.counts.CountsWriter;
import it.unimi.dsi.fastutil.ints.IntList;

import java.io.IOException;
import java.util.Map;

/**
 * @author Fabien Campagne
 *         Date: May 16, 2010
 *         Time: 12:55:53 PM
 */
public interface ComputeCountInterface {
    /**
     * This method must be called before calling the populate method. It initializes data structures.
     */
    void startPopulating();

    /**
     * Accumulate start and end counts to produce cumulative start and end count.
     * Pre-condition: the data structures starts and ends must have been populated (see method populate).
     * Post-condition: the data structures starts and end now contain the cumulative start and end counts.
     * It is bad design to reuse the same data structure to store different information, but is there a
     * significant performance advantage in this case?
     */
    void accumulate();

    /**
     * Calculate base counts and write the result to the specified CountsWriter.
     *
     * @param writer A CountsWriter to write the counts to.
     */
    void baseCount(final CountsWriter writer) throws IOException;

    /**
     * Calculates base counts. Stores the result in the countPerBase and countKey maps.
     */
    void baseCount();

    /**
     * Returns a map from position on the reference sequence to count per base.
     *
     * @return
     */
    Map getCountPerBase();

    /**
     * Return a list of positions on the reference sequence that have non-zero counts.
     */
    IntList getCountKeys();

    /**
     * Populate this data structure with a read.
     *
     * @param startPosition position where the read starts on the reference.
     * @param endPosition   position where the read ends on the reference.
     * @param forwardStrand True when the read matches the forward strand.
     * @param queryIndex    Index of the read/query.
     */
    void populate(int startPosition, int endPosition, boolean forwardStrand, int queryIndex);

    /**
     * This implementation ignores strand, but some sub-classes need this information.
     *
     * @param startPosition Start position of a read.
     * @param endPosition   End position of a read.
     * @param forwardStrand True when the read matches the forward strand, false otherwise.
     */

    void populate(final int startPosition, final int endPosition, final boolean forwardStrand);


}
