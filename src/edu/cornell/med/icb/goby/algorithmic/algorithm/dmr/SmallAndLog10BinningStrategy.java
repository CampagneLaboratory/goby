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

/**
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 4:14 PM
 */
public final class SmallAndLog10BinningStrategy implements BinningStrategy {
    private static final long serialVersionUID = 9153261064066006838L;

    @Override
    public int getBinIndex(final double covariate) {
        if (covariate < 100) {
            return 0;
        } else {
            return (int) (StrictMath.log10(covariate)-1);
        }
    }

    @Override
    public int getLowerBound(final int binIndex) {
        if (binIndex == 0) {
            return 0;
        }
        else {
            return (int) StrictMath.pow(10, binIndex + 1);
        }
    }

    @Override
    public int getUpperBound(final int binIndex) {
        if (binIndex == 0) {
            return 99;
        }
        else {
            return (int) StrictMath.pow(10, binIndex + 2);
        }
    }

    @Override
    public int getMidpoint(int binIndex) {
        final int lowerBound = getLowerBound(binIndex);
        final int upperBound = getUpperBound(binIndex);
        return (upperBound - lowerBound) / 2 + lowerBound;
    }
}
