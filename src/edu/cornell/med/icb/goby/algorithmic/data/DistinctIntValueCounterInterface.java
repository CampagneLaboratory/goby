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

import it.unimi.dsi.fastutil.ints.IntCollection;

/**
 * Counts the number of distinct integers observed.
 *
 * @author Fabien Campagne
 *         Date: Mar 15, 2011
 *         Time: 5:18:26 PM
 */
public interface DistinctIntValueCounterInterface {
    /**
     * Observe a set of values.
     */
    public void observe(IntCollection values);

    /**
     * Observe a set of values.
     */
    public void observe(int[] values);

    /**
     * Observe a value.
     *
     * @param value
     */
    public void observe(int value);

    /**
     * Return the number of distinct integer values in the observed sequence.
     *
     * @return
     */
    public int count();
}
