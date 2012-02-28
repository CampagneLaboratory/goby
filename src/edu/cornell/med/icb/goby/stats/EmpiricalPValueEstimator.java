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

import edu.cornell.med.icb.goby.algorithmic.algorithm.*;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.*;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.SamplePairEnumerator;
import edu.cornell.med.icb.goby.util.DynamicOptionClient;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.log4j.Logger;

/**
 * Helper class to derive empirical p-values from discrete observed null distributions.
 *
 * @author Fabien Campagne
 * @since Goby 1.9.8.4
 *        Date: 2/26/12
 *        Time: 10:23 PM
 */
public class EmpiricalPValueEstimator {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(EmpiricalPValueEstimator.class);

    private StatisticAdaptor statAdaptor;
    private String serializedFilename;


    private SamplePairEnumerator groupEnumerator;
    private EvidenceCombinator combinator = new MaxCombinator();
    private DensityEstimator estimator;
    private boolean densityLoadedFromDisk;

    enum combinatorNames {
        max, sum, qfast, median
    }

    enum statisticNames {
        delta, stat4, stat5, dMR, fisher
    }

    enum binningStrategyNames {
        fastslog10, log2, s100linear
    }

    static public final String[] LOCAL_DYNAMIC_OPTIONS = {

            "estimate-intra-group-differences: boolean, true indicates that pair-wise differences for sample in the same group should be tallied and written to the output. False indicates regular output.:false",
            "estimate-empirical-P: boolean, true: activates estimation of the empirical p-value.:false",
            "combinator: string, the method to combine p-values, one of qfast, average, sum, max.:max",
            "serialized-estimator-filename: string, the path to a serialized version of the density estimator populated with the empirical null-distribution.:",
            "statistic: string, the name of the statistic to evaluate between pairs of samples, one of stat4,stat5,dMR:stat5"

    };

    /**
     * Obtain dynamic options from the client of this class and configure this estimator.  For configure to work, the
     * client class must have a static registered dynamicOptionClient that was constructed with  localDynamicOptions
     *
     * @param numberOfContexts Number of discrete contexts to use. One null distribution is estimated for each context.
     * @param clientDoc        parsed dynamic options for the client of this class.
     */
    public void configure(final int numberOfContexts, final DynamicOptionClient clientDoc) {

        Boolean estimateIntraGroupDifferences = clientDoc.getBoolean("estimate-intra-group-differences");
        Boolean estimateIntraGroupP = clientDoc.getBoolean("estimate-empirical-P");
        serializedFilename = clientDoc.getString("serialized-estimator-filename");
        if (estimateIntraGroupP && serializedFilename != null) {
            try {
                estimator = DensityEstimator.load(serializedFilename);
                statAdaptor = estimator.getStatAdaptor();
                densityLoadedFromDisk = true;
            } catch (Exception e) {
                throw new RuntimeException("Unable to load serialized density with filename=" + serializedFilename);
            }
        }
        String combinatorName = clientDoc.getString("combinator");

        try {

            switch (combinatorNames.valueOf(combinatorName)) {
                case max:
                    combinator = new MaxCombinator();
                    break;
                case sum:
                    combinator = new SummedCombinator();
                    break;
                case qfast:
                    combinator = new QFast();
                    break;
                case median:
                    combinator = new MedianCombinator();
                    break;
                default:
                    new InternalError("This combinator name is not properly handled: " + combinatorName);
            }
        } catch (IllegalArgumentException e) {
            LOG.error(String.format("The combinator name %s was not recognized, using the default combinator instead (max).", combinatorName));
            combinator = new MaxCombinator();
        }
        if (statAdaptor != null) {
            LOG.info("StatAdaptor was obtained from loaded density: " + statAdaptor.statName());
        } else {
            String statisticName = clientDoc.getString("statistic");
            try {
                switch (statisticNames.valueOf(statisticName)) {
                    case delta:
                        statAdaptor = new DeltaStatisticAdaptor();
                        break;
                    case stat4:
                        statAdaptor = new Stat4StatisticAdaptor();
                        break;
                    case dMR:
                        statAdaptor = new MethylationRateDifferenceStatisticAdaptor();
                        break;
                    case fisher:
                        statAdaptor = new FisherExactTestAdaptor();
                        break;
                    default:
                    case stat5:
                        statAdaptor = new Stat5StatisticAdaptor();
                        break;
                }
            } catch (IllegalArgumentException e) {
                LOG.error(String.format("The statistic name %s was not recognized, using the default statistic instead (stat5).", combinatorName));
                statAdaptor = new Stat5StatisticAdaptor();
            }
        }
        BinningStrategy binningStrategy = null;
        if (!densityLoadedFromDisk) {

            String binningStrategyName = clientDoc.getString("binning-strategy");
            if (binningStrategyName != null) {

                try {
                    switch (binningStrategyNames.valueOf(binningStrategyName)) {
                        case log2:
                            binningStrategy = new Log2BinningStrategy();
                            break;
                        case s100linear:
                            binningStrategy = new LinearBinningStrategy();
                            break;
                        default:
                        case fastslog10:
                            binningStrategy = new FastSmallAndLog10BinningStrategy();
                            break;


                    }

                } catch (IllegalArgumentException e) {
                    LOG.error(String.format("The statistic name %s was not recognized, using the default statistic instead (stat5).", combinatorName));
                    binningStrategy = new FastSmallAndLog10BinningStrategy();
                }
            }
        } else {
            binningStrategy = estimator.getBinningStrategy();
        }
        if (estimator == null) {

            estimator = new DensityEstimator(numberOfContexts, statAdaptor);
            if (!densityLoadedFromDisk) {
                estimator.setBinningStrategy(binningStrategy);
            }
        }
    }


    void observeBetweenGroupPair(final GroupComparison comparison) {
        groupEnumerator.recordPairForGroupComparison(comparison);
    }

    void observeIntraGroupPairs(final GroupComparison comparison) {
        //   groupEnumerator = new SamplePairEnumerator(this.sampleIndexToGroupIndex, numSamples, numGroups, 0);

        groupEnumerator.recordPairForGroup(comparison.indexGroup1);
        groupEnumerator.recordPairForGroup(comparison.indexGroup2);

    }

    /**
     * Return the p-value that the difference observed between any of the pair could have been generated
     * by the distribution represented in the estimated null distribution. In this method, we compare samples across groups,
     * and use a distribution derived from pairs of samples in the same group. We therefore estimate a p-value
     * where the null-hypothesis is that the difference observed was generated by intra-group variations.
     *
     * @param contextIndex A discrete context covariate. A different null is estimated for each context.
     * @param comparison   A group comparison of interest.
     * @param dataProvider A format field counter for methylation data.
     * @return p-value.
     */
    public double estimateEmpiricalPValue(final int contextIndex,
                                          final GroupComparison comparison,
                                          final Object dataProvider) {
        combinator.reset();
        final ObjectArrayList<SamplePair> pairs = groupEnumerator.getPairs(comparison);
        for (final SamplePair pair : pairs) {
            int scaledStatistic = calculateScaledStatistic(dataProvider, pair.sampleIndexA, pair.sampleIndexB, contextIndex);
            if (!statAdaptor.ignorePair()) {
                final int[] covariates = statAdaptor.pairCovariates();
                final double p = estimator.getP(scaledStatistic, covariates);
                combinator.observe(p);

            } else {
                combinator.observe(1.0);
            }

        }
        return combinator.adjust();
    }


    /**
     * Return the p-value that the difference observed between any of the pair could have been generated
     * by the distribution represented in the estimated null distribution. In this method, we compare samples across groups,
     * and use a distribution derived from pairs of samples in the same group. We therefore estimate a p-value
     * where the null-hypothesis is that the difference observed was generated by intra-group variations.
     *
     * @param valuesA     values for sample A
     * @param valuesB     values for sample B
     * @param covariatesA covariates for sample A
     * @param covariatesB covariates for sample B
     * @return p-value.
     */
    public double estimateEmpiricalPValue(final IntArrayList[] valuesA, final IntArrayList[] valuesB, final IntArrayList[] covariatesA, final IntArrayList[] covariatesB) {
        int num = valuesA.length;
        combinator.reset();
        for (int pairIndex = 0; pairIndex < num; pairIndex++) {
            statAdaptor.reset();
            final double unscaledStatistic = statAdaptor.calculate(valuesA[pairIndex], valuesB[pairIndex], covariatesA[pairIndex], covariatesB[pairIndex]);

            if (!statAdaptor.ignorePair()) {
                final int[] covariates = statAdaptor.pairCovariates();
                final double p = estimator.getP(unscaledStatistic, covariates);
                //  System.out.println("Observing " + p);
                combinator.observe(p);

            } else {
                combinator.observe(1.0);
            }

        }
        return combinator.adjust();
    }

    /**
     * Observe differences within a group and estimate discrete null distributions.
     *
     * @param contextIndex A discrete context covariate. A different null is estimated for each context.
     * @param groupIndex   The group to derive within group differences.
     * @param dataProvider A format field counter for methylation data.
     */
    public void estimateNullDensity(final int contextIndex, final int groupIndex,
                                    final Object dataProvider) {

        // enumerate sample pairs that belong to the group of interest:
        final ObjectArrayList<SamplePair> pairs = groupEnumerator.getPairs(groupIndex);
        for (final SamplePair next : pairs) {

            observe(dataProvider, next.sampleIndexA, next.sampleIndexB, contextIndex);

        }

    }

    /**
     * Observe within group differences and estimate discrete null distributions.
     *
     * @param valuesA     values for sample A
     * @param valuesB     values for sample B
     * @param covariatesA covariates for sample A
     * @param covariatesB covariates for sample B
     */
    public void estimateNullDensity(final IntArrayList valuesA, final IntArrayList valuesB, final IntArrayList covariatesA, final IntArrayList covariatesB) {
        statAdaptor.reset();
        final double unscaledStatistic = statAdaptor.calculate(valuesA, valuesB, covariatesA, covariatesB);
        if (!statAdaptor.ignorePair()) {
            final int scaledStatistic = (int) Math.round(unscaledStatistic * estimator.getScalingFactor());

            //System.out.printf("observing context=%d sumTotal=%d scaledStatistic=%d elementIndex=%d %n", contextIndex, sumTotal, scaledStatistic, elementIndex);
            estimator.getDensity(statAdaptor.pairCovariates()).incrementCount(scaledStatistic);
        }

    }

    public final void observe(final Object sampleDataPool, final int sampleIndexA, final int sampleIndexB, final int contextIndex) {
        statAdaptor.reset();
        final double unscaledStatistic = statAdaptor.calculate(sampleDataPool, sampleIndexA, sampleIndexB, contextIndex);
        if (!statAdaptor.ignorePair()) {
            final int scaledStatistic = (int) Math.round(unscaledStatistic * estimator.getScalingFactor());

            //System.out.printf("observing context=%d sumTotal=%d scaledStatistic=%d elementIndex=%d %n", contextIndex, sumTotal, scaledStatistic, elementIndex);
            estimator.getDensity(contextIndex).incrementCount(scaledStatistic);
        }
    }

    public final int calculateScaledStatistic(final Object sampleDataPool, int sampleIndexA, int sampleIndexB, int contextIndex) {
        statAdaptor.reset();
        final double unscaledStatistic = statAdaptor.calculate(sampleDataPool, sampleIndexA, sampleIndexB, contextIndex);
        if (statAdaptor.ignorePair()) {
            return Integer.MAX_VALUE;
        } else {
            return (int) Math.round(unscaledStatistic * estimator.getScalingFactor());
        }
    }

    public final int calculateScaledStatistic(final IntArrayList valuesA, final IntArrayList valuesB,
                                              final IntArrayList covariatesA, final IntArrayList covariatesB) {
        statAdaptor.reset();
        final double unscaledStatistic = statAdaptor.calculate(valuesA, valuesB, covariatesA, covariatesB);
        if (statAdaptor.ignorePair()) {
            return Integer.MAX_VALUE;
        } else {
            return (int) Math.round(unscaledStatistic * estimator.getScalingFactor());
        }
    }

    public SamplePairEnumerator getGroupEnumerator() {
        return groupEnumerator;
    }

    public void setGroupEnumerator(SamplePairEnumerator groupEnumerator) {
        this.groupEnumerator = groupEnumerator;
    }


    public EvidenceCombinator getCombinator() {
        return combinator;
    }

    public void setCombinator(EvidenceCombinator combinator) {
        this.combinator = combinator;
    }


    public DensityEstimator getEstimator() {
        return estimator;
    }

    public void setEstimator(DensityEstimator estimator) {
        this.estimator = estimator;
    }

    public StatisticAdaptor getStatAdaptor() {
        return statAdaptor;
    }

    public void setStatAdaptor(StatisticAdaptor statAdaptor) {
        this.statAdaptor = statAdaptor;

    }


}