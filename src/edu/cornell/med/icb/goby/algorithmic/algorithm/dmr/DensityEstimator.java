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
    BinningStrategy binningStrategy = new SmallAndLog10BinningStrategy();
    private StatisticAdaptor statAdaptor;

    public DensityEstimator(int numberOfContexts) {
        densities = new ObjectArrayList<FenwickTree>();
        statAdaptor = new DeltaStatisticAdaptor();
        SCALING_FACTOR = (int) Math.round(MAX_ITEMS / statAdaptor.getRange());
        binningStrategy = new Log10BinningStrategy();
    }

    public DensityEstimator(int numberOfContexts, StatisticAdaptor statAdaptor) {
        densities = new ObjectArrayList<FenwickTree>();
        this.statAdaptor = statAdaptor;
        SCALING_FACTOR = (int) Math.round(MAX_ITEMS / statAdaptor.getRange());
        binningStrategy = new Log10BinningStrategy();
    }

    /**
     * Factor by which the statistic will be scaled to use MAX_ITEMS buckets.
     */
    final int SCALING_FACTOR;

    public BinningStrategy getBinningStrategy() {
        return binningStrategy;
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
     */
    public void observe(int contextIndex, int... a) {
         int sumTotal = 0;
        for (int val : a) {
            sumTotal += val;

        }
        final int delta = (int) Math.round(statAdaptor.calculateWithCovariate(sumTotal, a) * SCALING_FACTOR);


        final int elementIndex = delta;//(int) (((((double)delta)/(1.0+sumTotal)*MAX_ITEMS*0.9)));
        //System.out.printf("observing context=%d sumTotal=%d delta=%d elementIndex=%d %n", contextIndex, sumTotal, delta, elementIndex);
        getDensity(contextIndex, sumTotal).incrementCount(elementIndex);

    }

    public void observeWithCovariate(int contextIndex, int sumTotal, int... a) {

        final int delta = (int) Math.round(statAdaptor.calculate(a) * SCALING_FACTOR);
        final int elementIndex = delta;//(int) (((((double)delta)/(1.0+sumTotal)*MAX_ITEMS*0.9)));
        //System.out.printf("observing context=%d sumTotal=%d delta=%d elementIndex=%d %n", contextIndex, sumTotal, delta, elementIndex);
        getDensity(contextIndex, sumTotal).incrementCount(elementIndex);

    }

    public void observe(int contextIndex, int delta, int sumTotal) {

        // System.out.printf("observing context=%d sumTotal=%d delta=%d elementIndex=%d %n", contextIndex, sumTotal, delta, delta);
        getDensity(contextIndex, sumTotal).incrementCount(delta);

    }


    private FenwickTree getDensity(final int contextIndex, final int sumTotal) {


        final int index = binningStrategy.getBinIndex(sumTotal);
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


    public static void store(final DensityEstimator estimator, final String filename) throws IOException {
        BinIO.storeObject(estimator, filename);
    }

    public static DensityEstimator load(final String filename) throws IOException, ClassNotFoundException {
        return (DensityEstimator) BinIO.loadObject(filename);
    }

    /**
     * Get the cumulative count for observations with delta between zero and the argument value.
     *
     * @param contextIndex    index of the context.
     * @param sumTotal        total sum.
     * @param scaledStatistic upper-bound on the scaledStatistic for counting observations.
     * @return
     */
    public long getCumulativeCount(int contextIndex, int sumTotal, int scaledStatistic) {
        final FenwickTree tree = getDensity(contextIndex, sumTotal);
        return tree.getCumulativeCount(scaledStatistic);
    }

    /**
     * Get the cumulative count for observations with delta between zero and the argument value (inclusive).
     *
     * @param contextIndex index of the context.
     * @param sumTotal     total sum.
     * @param statistic    upper-bound on the statistic for counting observations.
     * @return
     */
    public long getCumulativeCount(int contextIndex, int sumTotal, double statistic) {
        final FenwickTree tree = getDensity(contextIndex, sumTotal);
        return tree.getCumulativeCount(scale(statistic));
    }

    /**
     * Get the empirical estimate of the probability that the statistic could have been generated by the distribution
     * represented in the estimator.
     *
     * @param contextIndex index of the context.
     * @param sumTotal     total sum.
     * @param statistic    value of statistic under test.
     * @return
     */
    public double getP(final int contextIndex, final int sumTotal, final double statistic) {
        final int scaledStatistic = (int) Math.round(statistic * SCALING_FACTOR);
        final FenwickTree tree = getDensity(contextIndex, sumTotal);
        final long totalCount = tree.getTotalCount();
        final double r = totalCount - tree.getCumulativeCount(scaledStatistic);
        final double n = totalCount;
        // estimated as per Morgan, Linda Am. J. Hum. Genet. 71:439–441, 2002
        final double p = (r + 1.0d) / (n + 1.0d);
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
                String statName = estimated.getStatAdaptor().statName();
                outWriter.println("midPointSumTotal\tsumTotal range\t"+ statName +"\tcount-at-"+statName);
                final BinningStrategy binningStrategy = estimated.getBinningStrategy();
                for (final FenwickTree tree : estimated.densities) {
                    if (tree != null) {
                        int low = binningStrategy.getLowerBound(index);
                        int high = binningStrategy.getUpperBound(index);

                        int midPointSumTotal = binningStrategy.getMidpoint(index);
                        System.out.printf("low=%d high=%d midPoint=%d %n", low, high, midPointSumTotal);
                        final long maxCumulative = tree.getCumulativeCount(tree.size() - 2);
                        for (int scaledStatistic = 0; scaledStatistic < tree.size() - 1; scaledStatistic++) {
                            final long cumulativeCountAt = tree.getCumulativeCount(scaledStatistic);
                            final long cumulativeCountAfter = tree.getCumulativeCount(scaledStatistic + 1);

                            outWriter.printf("%d\t[%d-%d]\t%g\t%d%n", midPointSumTotal, binningStrategy.getLowerBound(index),
                                    binningStrategy.getUpperBound(index), estimated.unscale(scaledStatistic), cumulativeCountAfter - cumulativeCountAt);
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

    /**
     * unscale a scaled statistic.
     *
     * @param scaledStatistic Value of the statistic scaled with SCALING_FACTOR.
     * @return the raw statistic that was observed.
     */
    public double unscale(final int scaledStatistic) {
        return ((double) scaledStatistic) / (double) SCALING_FACTOR;
    }

    /**
     * unscale a scaled statistic.
     *
     * @param statistic Value of the statistic to scale.
     * @return the scaled statistic value.
     */
    public int scale(final double statistic) {
        return (int) Math.round(statistic * SCALING_FACTOR);
    }

    public StatisticAdaptor getStatAdaptor() {
        return statAdaptor;
    }

    public void setBinningStrategy(BinningStrategy theBinningStrategy) {
        binningStrategy = theBinningStrategy;
    }
}
