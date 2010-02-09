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

import gominer.Fisher;
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
                fisherExact.fexact(new int[] { 12, 4, 15, 17, 25, 4 }, 3, 2 );
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

    /**
     * An an example of an R x C table from Agresti (2002, p. 57) Job Satisfaction.
     */
    @Test
    public void agrestiJobSatisfaction() {
        final FisherExact fisherExact = new FisherExact();
        final int[] inputTable = {
                /*                           income
                /* satisfaction    <15k    15-25k    25-40k   >40k */
                /*  VeryD */        1,       2,        1,      0,
                /*  LittleD */      3,       3,        6,      1,
                /*  ModerateS */    10,      10,       14,     9,
                /*  VeryS */        6,       7,        12,     11
        };

        final FisherExact.Result result = fisherExact.fexact(inputTable, 4, 4);
        assertEquals("pValue does not match", 0.7826849389656096, result.getPValue(), EPSILON);

        // everything else should be invalid since the input was not a 2x2 matrix
        assertNotNull("Confidence interval should not be null", result.getConfidenceInterval());
        assertTrue("Confidence interval should be an empty array",
                ArrayUtils.isEmpty(result.getConfidenceInterval()));
        assertTrue("Estimate should be NaN", Double.isNaN(result.getEstimate()));
        assertTrue("Odds ratio should be NaN", Double.isNaN(result.getOddsRatio()));
        assertEquals("Wrong Hypothesis for result", FisherExact.AlternativeHypothesis.twosided,
                result.getAlternativeHypothesis());

    }

    /**
     * Fisher Tea Tasting example.  A British woman claimed to be able to distinguish
     * whether milk or tea was added to the cup first.  To test, she was given 8 cups of
     * tea, in four of which milk was added first.  The null hypothesis is that there is
     * no association between the true order of pouring and the women's guess, the
     * alternative that there is a positive association (that the odds ratio
     * is greater than 1).
     */
    @Test
    public void agrestiTeaTasting() {
        final FisherExact fisherExact = new FisherExact();
        final int[] inputTable = {
                /*          Truth */
                /* Guess    Milk  Tea */
                /* Milk */    3,   1,
                /*  Tea */    1,   3
        };

        final FisherExact.Result result =
                fisherExact.fexact(inputTable, 2, 2, FisherExact.AlternativeHypothesis.greater);
        assertEquals("pValue does not match", 0.24285714285714288, result.getPValue(), EPSILON);
        assertEquals("Lower confidence interval does not match", 0.313569264110218,
                result.getConfidenceInterval()[0], EPSILON);
        assertTrue("Upper confidence interval should be infinite",
                Double.isInfinite(result.getConfidenceInterval()[1]));
        assertEquals("Estimate does not match", 6.408308867005793,  result.getEstimate(), EPSILON);
        assertEquals("Odds ratio does not match", 1.0, result.getOddsRatio(), EPSILON);
        assertEquals("Wrong Hypothesis for result", FisherExact.AlternativeHypothesis.greater,
                result.getAlternativeHypothesis());
    }

    /**
     * Fisher (1962, 1970), Criminal convictions of like-sex twins.
     */
    @Test
    public void twinConvictions() {
        final FisherExact fisherExact = new FisherExact();
        final int[] inputTable = {
                /*                   Dizygotic   Monozygotic */
                /* Convicted */         2,           10,
                /* Not Convicted */     15,          3
        };

        final FisherExact.Result result =
                fisherExact.fexact(inputTable, 2, 2, FisherExact.AlternativeHypothesis.less);
        assertEquals("pValue does not match", 0.00046518094336290525, result.getPValue(), EPSILON);
        assertEquals("Lower confidence interval does not match", 0.0,
                result.getConfidenceInterval()[0], EPSILON);
        assertEquals("Upper confidence interval does not match", 0.2849601379355694,
                result.getConfidenceInterval()[1], EPSILON);
        assertEquals("Estimate does not match", 0.04693660882769885, result.getEstimate(), EPSILON);
        assertEquals("Odds ratio does not match", 1.0, result.getOddsRatio(), EPSILON);
        assertEquals("Wrong Hypothesis for result", FisherExact.AlternativeHypothesis.less,
                result.getAlternativeHypothesis());
    }

    /**
     * Validates that the R implementation returns the same p-value for a simple
     * example set.
     */
    @Test
    public void twoTailed() {
        final Fisher gominer = new Fisher();
        final double gominerPValue = gominer.fisher(40, 10, 100, 30);

        final FisherExact fisherExact = new FisherExact();
        final double fisherExactPValue = fisherExact.twoTailed(40, 10, 100, 30);
        assertEquals("R result does not match gominer", gominerPValue, fisherExactPValue, EPSILON);
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
