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

package edu.cornell.med.icb.goby.counts;

import java.io.Closeable;
import java.io.IOException;

/**
 * Interface for implementations that support writing counts.
 * @author Fabien Campagne
 *         Date: 10/29/11
 *         Time: 1:21 PM
 */
public interface CountsWriterI extends Closeable {

    /**
     * Return the number of bits written to the output.
     * @return  the number of bits written to the output.
     */
    long getNumberOfBitsWritten();
    /**
     * Return the number of transitions written to the output.
     * @return the number of transitions written to the output.
     */
    int getNumberOfTransitions();

    /**
     * Append a transition.
     * @param count The transition brings count to this value.
     * @param lengthConstant The counts are constant, with count, for this number of bases.
     * @throws IOException If an error occurs writing this transition.
     */
    void appendCount(int count, int lengthConstant) throws IOException;

    /**
     * {@inheritDoc}
     */
    void close() throws IOException;

    /**
     * The total number of bases seen in the counts data we wrote.
     * This is defined as the sum of count*length over all transitions written by this reader. A normalization
     * factor for count data can be defined as   getNumberOfBasesSeen()/  getNumberOfSitesSeen() : this represents
     * the average coverage per site observed.
     *
     * @return number of bases seen.
     */
    long getNumberOfBasesSeen();

    /**
     * The total number of sites observed at which count!=0.
     *
     * @return number of sites seen.
     */
    long getNumberOfSitesSeen();

    int getInitialCount();
}
