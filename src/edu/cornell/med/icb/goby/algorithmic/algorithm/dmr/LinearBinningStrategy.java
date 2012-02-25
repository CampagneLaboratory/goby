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
 *         Time: 4:20 PM
 */
public class LinearBinningStrategy implements BinningStrategy {


    private static final long serialVersionUID = -5907440456960816739L;
    private int BIN_SIZE_SUM_TOTAL = 1000;

       protected final int getTheIndex(final int sumTotal) {
        final int theIndex;
        for (int index = 0; ; index++) {
            if (getLowerBound(index) <= sumTotal && sumTotal < getUpperBound(index)) {
                theIndex = index;
                break;
            }
        }
        return theIndex;
    }

    @Override
    public int getBinIndex(double covariate) {
        return covariate < 100 ? 0 : (int) (covariate / BIN_SIZE_SUM_TOTAL + 1);
    }

    @Override
    public int getLowerBound(int binindex) {
        if (binindex == 0) {
            return 0;
        }
        if (binindex == 1) {
            return 100;
        }

        return BIN_SIZE_SUM_TOTAL * (binindex - 1);

    }

    @Override
    public int getUpperBound(int binindex) {
        if (binindex == 0) {
            return 100;
        }
        return BIN_SIZE_SUM_TOTAL * binindex;

    }

    @Override
    public int getMidpoint(int binIndex) {
      final  int lowerBound = getLowerBound(binIndex);
      final  int upperBound = getUpperBound(binIndex);
        return (upperBound - lowerBound) / 2 + lowerBound;
    }
}
