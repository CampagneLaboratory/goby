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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.goby.algorithmic.algorithm.EvidenceCombinator;

/**
 * Keep the min p-value observed.
 * @author Fabien Campagne
 *         Date: 5/1/12
 *         Time: 5:34 PM
 */
public class MinCombinator implements EvidenceCombinator {
    double min = 1;
    int numP;

    @Override
    public void observe(double pValue) {
        min = Math.min(pValue, min);
        numP += 1;
    }

    @Override
    public void reset() {
        numP = 0;
        min = 1.0;
    }

    @Override
    public double adjust() {
        if (numP > 0.0) {
            return min;
        } else {
            // no p-value was observed.
            return 1.0;
        }
    }
}