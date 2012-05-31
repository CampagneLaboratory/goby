/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.counts;

import java.io.Closeable;
import java.io.IOException;

/**
 * Minimum contract for all implementations that support iterating over counts.
 *
 * @author Fabien Campagne
 *         Date: Jun 15, 2009
 *         Time: 6:37:21 PM
 */
public interface CountsReaderI extends Closeable {
    /**
     * Return the zero-based position along the sequence where the count is observed. For instance
     * getPosition() will return 1 on the first transition of the toy histogram shown below (and 3
     * on the second transition):
     * <p/>
     * xx
     * 0123
     *
     * @return
     */
    int getPosition();

    /**
     * Determines if the reader has data about another transition.
     *
     * @return True when a call to nextTransition() will succeed, False otherwise.
     * @throws IOException
     */
    boolean hasNextTransition() throws IOException;

    /**
     * Advance to the next transition. After this method has been called successfully, position,
     * length, deltaCount and currentCount are available through getters of this reader.
     *
     * @throws IOException
     */
    void nextTransition() throws IOException;

    /**
     * Return the count at the given position.
     *
     * @return
     */
    int getCount();

    /**
     * Advance up to or past the specified position. The reader is advanced until the
     * position returned by {@link #getPosition()} is at least equal, or greater to
     * the specified position.
     *
     * @param position
     * @throws IOException
     */
    void skipTo(int position) throws IOException;

    /**
     * Reposition the reader on a genomic position. In contrast to skipTo, this method reposition to any position,
     * including position that occured before the position returned last by getPosition().
     *
     * @param position Position to seek to.
     * @throws IOException If an error occurs accessing the index or counts data.
     */
    void reposition(int position) throws java.io.IOException;

    /**
     * The length of the region/peak where the count is observed.
     *
     * @return The length of the region
     */
    int getLength();
}
