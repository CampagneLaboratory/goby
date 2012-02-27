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

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: 2/24/12
 *         Time: 2:09 PM
 */
public interface StatisticAdaptor extends Serializable {
    /**
     * Return the name of the statistic estimated by this adaptor.
     *
     * @return
     */
    public String statName();

    /**
     * Estimate a statistic given some integers.
     *
     * @param a values used to calculate the statistic. The meaning of the values and their order is implementation dependent.
     * @return the calculated statistic.
     */

    public double calculateNoCovariate(int... a);

    /**
     * Estimate a statistic given a covariate and some integers.
     *
     * @param covariate a relevant covariate of the statistic. Some calculator may use the covariate to normalize the
     *                  value of the statistic.
     * @param a         values used to calculate the statistic. The meaning of the values and their order is implementation dependent.
     * @return the calculated statistic.
     */

    public double calculateWithCovariate(int covariate, int... a);

    /**
     * Estimate a statistic given values and covariates for two samples A and B.
     */

    public double calculate(IntArrayList valuesSampleA, IntArrayList valuesSampleB, IntArrayList covariatesSampleA, IntArrayList covariatesSampleB);

    /**
     * Return the maximum value that the statistic can attain.
     *
     * @return a maximum bound on statistic.
     */
    double getMaximumStatistic();

    /**
     * Get the range of this statistic. The difference between the maximum statistic that can be calculated and the
     * minimum value.
     *
     * @return range of the statistic.
     */
    double getRange();

    /**
     * Calculate a statistic.
     *
     * @param dataProvider
     * @param sampleIndexA
     * @param sampleIndexB
     * @param covariate
     * @return NaN if the pair of samples should be ignored.
     */
    double calculate(Object dataProvider, int sampleIndexA, int sampleIndexB, int... covariate);

    /**
     * Return the covariates associated with the pair for which the statistic is estimated.
     *
     * @return an array of integers, where each element is a covariate associated with a pair of samples.
     */
    int[] pairCovariates();

    void reset();

    boolean ignorePair();

    /**
     * Set an optional observation writer on this statistic. Observations will be written if a writer is set.
     * @param obsWriter
     */
    public void setObservationWriter(ObservationWriter obsWriter);
}
