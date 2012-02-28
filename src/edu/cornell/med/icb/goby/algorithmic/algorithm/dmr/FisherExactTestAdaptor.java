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
import org.rosuda.JRI.Rengine;

/**
 * Computes -log10(fisher exact p-value)
 * @author Fabien Campagne
 *         Date: 2/28/12
 *         Time: 8:46 AM
 */
public class FisherExactTestAdaptor extends AbstractMethylationAdapter {
    private static final double MAXIMUM_BOUND = -Math.log10(Double.MIN_VALUE);
    private static final long serialVersionUID = -4127089751953478896L;
    private boolean fisherRInstalled;
    boolean ignorePair = false;

    public FisherExactTestAdaptor() {
//activate R
        try {
            final Rengine rEngine = GobyRengine.getInstance().getRengine();
            fisherRInstalled = rEngine != null && rEngine.isAlive();
        } catch (java.lang.UnsatisfiedLinkError e) {
            System.out.println("Cannot initialize R");
            e.printStackTrace();
            throw e;
        }

    }

    @Override
    public String statName() {
        return "fisher";
    }

    @Override
    public double calculateNoCovariate(int... a) {
        final int cma = a[0];
        final int ca = a[1];
        final int cmb = a[2];
        final int cb = a[3];
        double fisherP = Double.NaN;
        if (!fisherRInstalled) {
            setIgnorePair(true);
        }

        fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                ca,
                cma,
                cb,
                cmb) : Double.NaN;
        return -Math.log10(fisherP);
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
