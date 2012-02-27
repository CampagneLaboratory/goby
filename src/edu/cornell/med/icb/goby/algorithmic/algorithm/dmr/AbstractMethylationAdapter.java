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

import edu.cornell.med.icb.goby.stats.FormatFieldCounter;
import org.apache.log4j.Logger;

/**
 * @author Fabien Campagne
 *         Date: 2/26/12
 *         Time: 12:58 PM
 */
public abstract class AbstractMethylationAdapter implements StatisticAdaptor {
    private static final long serialVersionUID = 8344532506747517703L;
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(AbstractMethylationAdapter.class);
    private boolean ignorePair;


    @Override
    public double calculate(Object dataProvider, int sampleIndexA, int sampleIndexB, int... covariate) {
        if (!(dataProvider instanceof FormatFieldCounter)) {
            throw new InternalError();
        }
        FormatFieldCounter sampleDataPool = (FormatFieldCounter) dataProvider;
        final int contextIndex = covariate[0];
        final int cma = sampleDataPool.getMethylatedCCountPerSample(contextIndex, sampleIndexA);
        final int ca = sampleDataPool.getUnmethylatedCCountPerSample(contextIndex, sampleIndexA);
        final int cmb = sampleDataPool.getMethylatedCCountPerSample(contextIndex, sampleIndexB);
        final int cb = sampleDataPool.getUnmethylatedCCountPerSample(contextIndex, sampleIndexB);

        if ((cma + ca) == 0 || (cmb + cb) == 0) {
            if (cma + ca + cmb + cb != 0) {
                if (LOG.isTraceEnabled()) {
                    LOG.trace(String.format("Zero in one intra-group sample for %d %d %d %d samplexIndexA=%d sampleIndexB=%d %n",
                            cma, ca, cmb, cb, sampleIndexA, sampleIndexB));
                }
            }
            ignorePair = true;
            return Double.NaN;
        } else {
            ignorePair = false;
        }
        // remove context index from the list of covariates. Context is irrelevant to calculate the statistic.

        final int sumTotal = cma + cmb + ca + cb;
        COVARIATES[0] = sumTotal;
        return calculateWithCovariate(sumTotal, cma, ca, cmb, cb);
    }

    private final int[] COVARIATES = new int[1];

    @Override
    public final int[] covariates() {
        if (ignorePair) {
            return null;
        }
        return COVARIATES;
    }

    @Override
    public void reset() {
        ignorePair = false;
    }

    @Override
    public boolean ignorePair() {
        return ignorePair;
    }
}
