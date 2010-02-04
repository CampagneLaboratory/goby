/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.R;

import org.apache.commons.lang.ArrayUtils;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;
import org.rosuda.JRI.Rengine;

/**
 * Validates the functionality of the {@link edu.cornell.med.icb.goby.R.FisherExact} class.
 */
public class TestFisherExact {
    /**
     * Default epsilon for double comparsions.
     */
    private static final double EPSILON = 1.0e-12;

    /**
     * Example data from
     * <href="http://darwin.eeb.uconn.edu/eeb348/supplements-2006/chi-squared/chi-squared.html">
     * http://darwin.eeb.uconn.edu/eeb348/supplements-2006/chi-squared/chi-squared.html</a>
     */
    @Test
    public void chiSquaredUConn() {
        final FisherExact fisherExact = new FisherExact();
        final FisherExact.Result result =
                fisherExact.fexact(new int[] { 12, 17, 4, 25, 15, 4 }, 3, 2 );
        assertEquals("pValue does not match", 2.9565806126420623e-05, result.getPValue(), EPSILON);
    }

    /**
     * Validates the same input given to different convenience methods produces the same result.
     */
    @Test
    public void twoByTwoEquivalent() {
        final FisherExact fisherExact = new FisherExact();
        final FisherExact.Result result1 = fisherExact.fexact(new int[] {10, 20, 30, 40}, 2, 2);
        final FisherExact.Result result2 = fisherExact.fexact(10, 20, 30, 40);

        assertEquals("pValue does not match",
                result1.getPValue(), result2.getPValue(), EPSILON);
        assertEquals("Lower confidence interval does not match",
                result1.getConfidenceInterval()[0], result2.getConfidenceInterval()[0], EPSILON);
        assertEquals("Upper confidence interval does not match",
                result1.getConfidenceInterval()[1], result2.getConfidenceInterval()[1], EPSILON);
        assertEquals("Estimate does not match",
                result1.getEstimate(), result2.getEstimate(), EPSILON);
        assertEquals("Odds ratio does not match"
                , result1.getOddsRatio(), result2.getOddsRatio(), EPSILON);
    }

    @Test
    public void agrestiJobSatisfaction() {
        final FisherExact fisherExact = new FisherExact();
        final int[] inputTable = {
        /* income      VeryD  LittleD  ModerateS VeryS */
        /*  < 15k */    1,       3,        10,     6,
        /* 15-25k */    2,       3,        10,     7,
        /* 25-40k */    1,       6,        14,    12,
        /*  > 40k */    0,       1,         9,    11
        };

        final FisherExact.Result result = fisherExact.fexact(inputTable, 4, 4);
        assertEquals("pValue does not match", 0.7826849389656096, result.getPValue(), EPSILON);

        // everything else should be invalid since the input was not a 2x2 matrix
        assertNotNull("Confidence interval should not be null", result.getConfidenceInterval());
        assertTrue("Confidence interval should be an empty array",
                ArrayUtils.isEmpty(result.getConfidenceInterval()));
        assertTrue("Estimate should be NaN", Double.isNaN(result.getEstimate()));
        assertTrue("Odds ratio should be NaN", Double.isNaN(result.getOddsRatio()));
    }

    /**
     * If R libraries are not set up properly these tests cannot be run.
     */
    @BeforeClass
    public static void assertRAvailable() {
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        assertNotNull("R engine is not available", rengine);
        assertTrue("R is not null but is not alive either", rengine.isAlive());
    }

    /**
     * Notify {@link org.rosuda.JRI.Rengine} that the thread can be safely terminated.
     */
    @AfterClass
    public static void terminateRThread() {
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null) {
            rengine.end();
        }
    }
}
