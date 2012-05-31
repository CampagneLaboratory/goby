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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

/**
 * Keep the max p-value observed.
 * @author Fabien Campagne
 *         Date: 2/23/12
 *         Time: 2:59 PM
 */
public class MaxCombinator implements EvidenceCombinator {
    double max = 0;
    int numP;

    @Override
    public void observe(double pValue) {
        max = Math.max(pValue, max);
        numP += 1;
    }

    @Override
    public void reset() {
        numP = 0;
        max = 0.0;
    }

    @Override
    public double adjust() {
        if (numP > 0.0) {
            return max;
        } else {
            // no p-value was observed.
            return 1.0;
        }
    }
}
