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

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;

/**
 * Provider that estimates a p-value over a region.
 * @author Fabien Campagne
 *         Date: 5/6/12
 *         Time: 2:53 PM
 */
public interface PVAlueProvider {
    /**
     *
     * @param start Start coordinate of the region under test.
     * @param end   End coordinate of the region under test.
     * @param groupComp groups under test.
     * @return p-value estimated over the region, measuring the chance that the difference observed  in the region is due to chance.
     */
     public double getPValue(final int start, final int end, final GroupComparison groupComp);

}
