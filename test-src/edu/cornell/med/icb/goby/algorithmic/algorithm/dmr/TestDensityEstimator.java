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

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Fabien Campagne
 *         Date: 2/20/12
 *         Time: 4:26 PM
 */
public class TestDensityEstimator {
    class value {
        int context;
        int cma;
        int ca;
        int cmb;
        int cb;

        value(int context, int cma, int ca, int cmb, int cb) {
            this.context = context;
            this.cma = cma;
            this.ca = ca;
            this.cmb = cmb;
            this.cb = cb;
        }
    }

    @Test
    public void testObserve() throws Exception {
        final DensityEstimator estimator = new DensityEstimator(2);
        final value[] observations = {
                new value(0, 5, 13, 4, 12),
                new value(0, 0, 1, 1, 2),
                new value(0, 0, 1, 0, 3),
        };
        for (final value obs : observations) {
            estimator.observe(obs.context, obs.cma, obs.ca, obs.cmb, obs.cb);
        }
        assertEquals("", 2l, estimator.getCumulativeCount(0, 0, 0));
        assertEquals("", 3l, estimator.getCumulativeCount(0, 0, 2));
    }
}
