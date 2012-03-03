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
import it.unimi.dsi.fastutil.ints.IntArrayList;
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
    public double calculate(IntArrayList valuesSampleA, IntArrayList valuesSampleB, IntArrayList covariatesSampleA, IntArrayList covariatesSampleB) {
        // we flatten the contextIndex, it is the same for a pair of samples between we always compare the same site.

        final int cma = valuesSampleA.get(0);
        final int ca = valuesSampleA.get(1);
        final int cmb = valuesSampleB.get(0);
        final int cb = valuesSampleB.get(1);

        if ((cma + ca) == 0 || (cmb + cb) == 0) {
            if (cma + ca + cmb + cb != 0) {
                if (LOG.isTraceEnabled()) {
                    LOG.trace(String.format("Zero in one intra-group sample for %d %d %d %d %n",
                            cma, ca, cmb, cb));
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

    public void setIgnorePair(boolean b) {
        ignorePair = b;

    }

    @Override
    public double calculate(final Object dataProvider, final int sampleIndexA, final int sampleIndexB, final int... covariate) {
        checkProviderType(dataProvider);
        final int contextIndex = covariate[0];
        final int cma = getCm(sampleIndexA, dataProvider, contextIndex);
        final int ca = getC(sampleIndexA, dataProvider, contextIndex);
        final int cmb = getCm(sampleIndexB, dataProvider, contextIndex);
        final int cb = getC(sampleIndexB, dataProvider, contextIndex);
        recordObservation(cma, ca, cmb, cb);
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

    protected void checkProviderType(Object dataProvider) {
        if (!(dataProvider instanceof FormatFieldCounter)) {
            throw new InternalError();
        }
    }

    protected int getC(int sampleIndexA, Object sampleDataPool, int contextIndex) {
        final FormatFieldCounter counter = (FormatFieldCounter) sampleDataPool;
        return counter.getUnmethylatedCCountPerSample(contextIndex, sampleIndexA);
    }

    protected int getCm(int sampleIndexA, Object sampleDataPool, int contextIndex) {
        final FormatFieldCounter counter = (FormatFieldCounter) sampleDataPool;
        return counter.getMethylatedCCountPerSample(contextIndex, sampleIndexA);
    }


    private IntArrayList valuesA = new IntArrayList();
    private IntArrayList valuesB = new IntArrayList();
    private IntArrayList covariatesA = new IntArrayList();
    private IntArrayList covariatesB = new IntArrayList();
    private transient ObservationWriter obsWriter;
    String[] zero = new String[0];

    public void setObservationWriter(ObservationWriter obsWriter) {
        this.obsWriter = obsWriter;
       if (obsWriter!=null) {
        obsWriter.writeHeader(new String[]{"Cma", "Ca"}, new String[]{"Cmb", "Cb"}, zero, zero);
       }
    }

    void recordObservation(int cma, int ca, int cmb, int cb) {
        if (obsWriter != null) {

            valuesA.add(cma);
            valuesA.add(ca);
            valuesB.add(cmb);
            valuesB.add(cb);
            obsWriter.observed(valuesA, valuesB, covariatesA, covariatesB);
        }
    }

    private final int[] COVARIATES = new int[1];

    @Override
    public final int[] pairCovariates() {
        if (ignorePair) {
            return null;
        }
        return COVARIATES;
    }

    @Override
    public void reset() {
        ignorePair = false;
        if (obsWriter != null) {

            valuesA.clear();
            valuesB.clear();
            covariatesA.clear();
            covariatesB.clear();
        }
    }

    @Override
    public boolean ignorePair() {
        return ignorePair;
    }
}
