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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 4:11 PM
 */
public interface BinningStrategy  extends Serializable {
    public int getBinIndex(double covariate);

    /**
     * Get the lower bound on the covariate in the given bin.
     *
     * @param binIndex
     * @return
     */
    public int getLowerBound(int binIndex);

    /**
     * Get the upper bound on the covariate in the given bin (exclusive).
     *
     * @param binIndex
     * @return
     */
    public int getUpperBound(int binIndex);

    /**
     * Get the midPoint on the covariate in the given bin.
     *
     * @param binIndex
     * @return
     */
    public int getMidpoint(int binIndex);

    /**
     * Return a short name for this strategy.
     *
     * @return
     */
    public String getName();
}
