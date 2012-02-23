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
 * A combinator that sums individual p-values and clips the result at 1.0.
 *
 * @author Fabien Campagne
 *         Date: 2/23/12
 *         Time: 2:47 PM
 */
public class SummedCombinator implements EvidenceCombinator {
    double sum;
    boolean accumulate;
    int numP = 0;

    @Override
    public void observe(double pValue) {

        if (accumulate)
            if (sum + pValue <= 1) {
                sum += pValue;
            } else {
                //we stop accumulating if adding the p-value would result in q>1.
                accumulate = false;
            }
    }

    @Override
    public void reset() {
        sum = 0;
        accumulate = true;
        numP = 0;
    }

    @Override
    public double adjust() {
        if (numP > 0) {
            return Math.min(1.0, sum);
        } else {
            return 1.0;
        }
    }
}
