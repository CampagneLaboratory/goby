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
 * Implements the QFAST algorithm, described in  Combining multiple p-values: application to sequence
 * homology searches. Bailey TL, Gribskov M. Bioinformatics 1998.
 *
 * @author Fabien Campagne
 *         Date: 2/23/12
 *         Time: 12:32 PM
 */
public class QFast implements EvidenceCombinator {
    int numP = 0;
    double logProduct = 0;

    @Override
    public void observe(double pValue) {

        logProduct += StrictMath.log(pValue);
        numP += 1;

    }

    @Override
    public void reset() {
        numP = 0;
        logProduct = 0;
    }

    /**
     * Return the adjusted p-value after a series of p-values have been observed.
     * @return  the QFAST adjusted p-value.
     */
    @Override
    public double adjust() {
       if (numP>0) {
        final double product = StrictMath.exp(logProduct);
        return qfast(numP, product);
       } else {
           return 1.0;
       }
    }

    public static double qfast(final int numP, final double productP) {
        if (productP == 0.0) {
            return 0.0;
        }

        if (numP > 1) {
            double t = productP;
            double q = productP;
            final double x = -StrictMath.log(productP);
            for (int i = 1; i < numP - 1; i++) {
                t = t * x / i;
                q = q + t;
            }
            return q;
        } else {
            return productP;
        }
    }
}
