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

import edu.cornell.med.icb.goby.algorithmic.algorithm.FenwickTree;
import it.unimi.dsi.fastutil.io.BinIO;

import java.io.IOException;
import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: 2/19/12
 *         Time: 3:01 PM
 */
public class DensityEstimator implements Serializable{

    private static final long serialVersionUID = -4803501043413548993L;
    private static final int MAX_ITEMS = 10000;
    FenwickTree density[];

    public DensityEstimator(int numberOfContexts) {
        density = new FenwickTree[numberOfContexts];
        for (int contextIndex = 0; contextIndex < numberOfContexts; contextIndex++) {
            density[contextIndex] = new FenwickTree(MAX_ITEMS);
        }
    }

    /**
     * ca=5 cma=10  diffA=10-5=5
     * cb=6 cmb=13  diffB=13-6=7
     * diffA,B= 5-7=2
     * <p/>
     * ca=10 cma=5  diffA=10-5=5
     * cb=6  cmb=13 diffB=6-13=-7
     * diffA,B= 5- -7=13
     *
     * @param contextIndex
     * @param cma
     * @param ca
     * @param cmb
     * @param cb
     */
    public void observe(int contextIndex, int cma, int ca, int cmb, int cb) {
        int maxA, maxB;
        int minA, minB;

        if (cma > ca) {
            maxA = cma;
            minA = ca;
            maxB = cmb;
            minB = cb;
        } else {
            maxA = ca;
            maxB = cb;
            minA = cma;
            minB = cmb;
        }
        final int diffA = maxA - minA;
        final int diffB = maxB - minB;
        final int delta = Math.abs(diffA - diffB);
        final int elementIndex = delta;
        int sumTotal = cma + ca + cmb + cb;
        System.out.printf("observing context=%d sumTotal=%d delta=%d%n", contextIndex, sumTotal, delta);
        getDensity(contextIndex, sumTotal).incrementCount(elementIndex);

    }

    private FenwickTree getDensity(final int contextIndex, final int sumTotal) {
        return density[contextIndex];
    }

    public static void store(final DensityEstimator estimator, final String filename) throws IOException {
        BinIO.storeObject(estimator,filename);
    }
    public static DensityEstimator load( final String filename) throws IOException, ClassNotFoundException {
        return (DensityEstimator) BinIO.loadObject(filename);
    }

    /**
     * Get the cumulative count for observations with delta between zero and the argument value.
     * @param  contextIndex index of the context.
     * @param  sumTotal     total sum.
     * @param  delta upper-bound on delta for counting observations.
     * @return
     */
    public long getCumulativeCount(int contextIndex, int sumTotal, int delta) {
        final FenwickTree tree = getDensity(contextIndex, sumTotal);
        return tree.getCumulativeCount(delta);
    }
}
