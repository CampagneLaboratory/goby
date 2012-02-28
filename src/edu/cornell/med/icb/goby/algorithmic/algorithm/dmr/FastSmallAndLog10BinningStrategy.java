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
 * A fast implementation of SmallAndLog10BinningStrategy that assumes the covariate is not larger than Integer.MAX_INT.
 * This implementation uses a divide and conquer strategy
 *
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 4:14 PM
 */
public final class FastSmallAndLog10BinningStrategy implements BinningStrategy {
    private static final long serialVersionUID = 9153261064066006838L;


    @Override
    /**
     * Returns an index between 0 and 8, inclusive.
     */
    public int getBinIndex(final double covariate) {
        final int n = (int) covariate;
        int index;
        if (n < 100000) {
            // 5 or less
            if (n < 100) {
                return 0;
            } else {
                // 3 or 4 or 5
                if (n < 1000) {
                    return 3 - 2;
                } else {
                    // 4 or 5
                    if (n < 10000) {
                        return 4 - 2;
                    } else {
                        return 5 - 2;
                    }
                }
            }
        } else {
            // 6 or more
            if (n < 10000000) {
                // 6 or 7
                if (n < 1000000) {
                    return 6 - 2;
                } else {
                    return 7 - 2;
                }
            } else {
                // 8 to 10
                if (n < 100000000) {
                    return 8 - 2;
                } else {
                    // 9 or 10
                    if (n < 1000000000) {
                        return 9 - 2;
                    } else {
                        return 10 - 2;
                    }
                }
            }
        }

    }

    @Override
    public int getLowerBound(final int binIndex) {
        if (binIndex == 0) {
            return 0;
        } else {
            return (int) StrictMath.pow(10, binIndex + 1);
        }
    }

    @Override
    public int getUpperBound(final int binIndex) {
        if (binIndex == 0) {
            return 99;
        } else {
            return (int) StrictMath.pow(10, binIndex + 2);
        }
    }

    @Override
    public int getMidpoint(int binIndex) {
        final int lowerBound = getLowerBound(binIndex);
        final int upperBound = getUpperBound(binIndex);
        return (upperBound - lowerBound) / 2 + lowerBound;
    }

    @Override
    public String getName() {
        return "fs100log10";
    }
}
