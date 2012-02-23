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
import edu.mssm.crover.cli.CLI;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: 2/19/12
 *         Time: 3:01 PM
 */
public class DensityEstimator implements Serializable {

    private static final long serialVersionUID = -4803501043413548993L;
    private static final int MAX_ITEMS = 10000;
    private ObjectArrayList<FenwickTree> densities;
    private int BIN_SIZE_SUM_TOTAL = 1000;


    public DensityEstimator(int numberOfContexts) {
        densities = new ObjectArrayList<FenwickTree>();

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
        final int delta = getDelta(cma, ca, cmb, cb);

        int sumTotal = cma + ca + cmb + cb;
        final int elementIndex = delta;//(int) (((((double)delta)/(1.0+sumTotal)*MAX_ITEMS*0.9)));
        //System.out.printf("observing context=%d sumTotal=%d delta=%d elementIndex=%d %n", contextIndex, sumTotal, delta, elementIndex);
        getDensity(contextIndex, sumTotal).incrementCount(elementIndex);

    }
    public void observe(int contextIndex, int delta, int sumTotal) {

       // System.out.printf("observing context=%d sumTotal=%d delta=%d elementIndex=%d %n", contextIndex, sumTotal, delta, delta);
        getDensity(contextIndex, sumTotal).incrementCount(delta);

    }
    public int getDelta(int cma, int ca, int cmb, int cb) {
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
        return Math.abs(diffA - diffB);
    }

    private FenwickTree getDensity(final int contextIndex, final int sumTotal) {


        final int index = getTheIndexFast(sumTotal);
        while (densities.size() <= index) {
            densities.add(null);
        }
        final FenwickTree tree = densities.get(index);
        if (tree != null) {
            return tree;
        } else {
            // grow the array as needed:

            final FenwickTree newTree = new FenwickTree(MAX_ITEMS);
            densities.set(index, newTree);
            return newTree;
        }
    }

    protected final int getTheIndexFast(final int sumTotal) {
        return sumTotal < 100 ? 0 : sumTotal / BIN_SIZE_SUM_TOTAL + 1;
    }

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

    public int getLowerBound(int index) {
        if (index == 0) {
            return 0;
        }
        if (index == 1) {
            return 100;
        }

        return BIN_SIZE_SUM_TOTAL * (index - 1);

    }

    public int getUpperBound(int index) {
        if (index == 0) {
            return 100;
        }
        return BIN_SIZE_SUM_TOTAL * index;

    }

    public static void store(final DensityEstimator estimator, final String filename) throws IOException {
        BinIO.storeObject(estimator, filename);
    }

    public static DensityEstimator load(final String filename) throws IOException, ClassNotFoundException {
        return (DensityEstimator) BinIO.loadObject(filename);
    }

    /**
     * Get the cumulative count for observations with delta between zero and the argument value.
     *
     * @param contextIndex index of the context.
     * @param sumTotal     total sum.
     * @param delta        upper-bound on delta for counting observations.
     * @return
     */
    public long getCumulativeCount(int contextIndex, int sumTotal, int delta) {
        final FenwickTree tree = getDensity(contextIndex, sumTotal);
        return tree.getCumulativeCount(delta);
    }

    /**
     * Get the empirical estimate of the probability that the delta could have been generated by the distribution
     * represented in the estimator.
     *
     * @param contextIndex index of the context.
     * @param sumTotal     total sum.
     * @param delta        value of delta under test.
     * @return
     */
    public double getP(final int contextIndex, final int sumTotal, final int delta) {
        final FenwickTree tree = getDensity(contextIndex, sumTotal);
        final long totalCount = tree.getTotalCount();
        final double r=totalCount- tree.getCumulativeCount(delta) ;
        final double n= totalCount;
        // estimated as per Morgan, Linda Am. J. Hum. Genet. 71:439–441, 2002
        final double p=(r+1.0d) / (n+1.0d);
        return p;
    }

    public static void main(final String[] args) throws IOException {
        boolean printDensity = CLI.isKeywordGiven(args, "--print-density");
        String filename = CLI.getOption(args, "-f", null);
        String outputFilename = CLI.getOption(args, "-o", "out.tsv");
        PrintWriter outWriter = new PrintWriter(new FileWriter(outputFilename));
        if (printDensity) {
            DensityEstimator estimated = null;
            try {
                estimated = load(filename);
                int index = 0;
                outWriter.println("midPointSumTotal\tdelta\tcount-at-delta");
                for (final FenwickTree tree : estimated.densities) {
                    if (tree != null) {
                        int low = estimated.getLowerBound(index);
                        int high = estimated.getUpperBound(index);

                        int midPointSumTotal = low + ((high - low) / 2);
                        System.out.printf("low=%d high=%d midPoint=%d %n", low, high, midPointSumTotal);
                        final long maxCumulative = tree.getCumulativeCount(tree.size() - 2);
                        for (int delta = 0; delta < tree.size() - 1; delta++) {
                            final long cumulativeCountAt = tree.getCumulativeCount(delta);
                            final long cumulativeCountAfter = tree.getCumulativeCount(delta + 1);

                            outWriter.printf("%d\t%d\t%d%n", midPointSumTotal, delta, cumulativeCountAfter - cumulativeCountAt);
                            if (cumulativeCountAfter == maxCumulative) {
                                break;
                            }
                        }

                    }
                    index++;

                }
                outWriter.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        }
    }

}
