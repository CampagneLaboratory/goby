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

import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import org.apache.commons.math.MathException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

/**
 * Computes -log10(chi square p-value)
 *
 * @author Fabien Campagne
 *         Date: 2/28/12
 *         Time: 8:46 AM
 */
public class ChiSquareTestAdaptor extends AbstractMethylationAdapter {
    private static final double MAXIMUM_BOUND = -Math.log10(Double.MIN_VALUE);
    private static final long serialVersionUID = -4127089751953478896L;

    boolean ignorePair = false;
    private static final Logger LOG = Logger.getLogger(ChiSquareTestAdaptor.class);

    public ChiSquareTestAdaptor() {

    }


    @Override
    public String statName() {
        return "chisquare";
    }

    double[] expectedCounts = new double[2];
    long[] observedCounts = new long[2];

    @Override
    public double calculateNoCovariate(int... a) {
        final int cma = a[0];
        final int ca = a[1];
        final int cmb = a[2];
        final int cb = a[3];

        if (cma == 0 || ca == 0 || cmb == 0 || cb == 0) {
            setIgnorePair(true);
            return -StrictMath.log10(1.0);
        }
        expectedCounts[0] = cma;
        expectedCounts[1] = ca;
        observedCounts[0] = cmb;
        observedCounts[1] = cb;
        double pValue = Double.NaN;
        final ChiSquareTest chisquare = new ChiSquareTestImpl();
        try {
            final double pValueRaw = chisquare.chiSquareTest(expectedCounts, observedCounts);
            // math commons can return negative p-values?
            pValue = Math.abs(pValueRaw);
        } catch (MaxIterationsExceededException e) {

            LOG.error("expected:" + DoubleArrayList.wrap(expectedCounts).toString());
            LOG.error("observed:" + LongArrayList.wrap(observedCounts).toString());
            LOG.error(e);
            pValue = 1.0;
            setIgnorePair(true);
        } catch (MathException e) {
            e.printStackTrace();
            setIgnorePair(true);

        }

        setIgnorePair(false);
        return -StrictMath.log10(pValue);
    }


    @Override
    /**
     *
     */
    public double calculateWithCovariate(final int covariate, final int... a) {
        return Math.min(calculateNoCovariate(a), MAXIMUM_BOUND);

    }

    @Override
    public double getMaximumStatistic() {
        return MAXIMUM_BOUND;
    }

    @Override
    public double getRange() {
        return MAXIMUM_BOUND;
    }


}
