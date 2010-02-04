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
import org.apache.commons.lang.BooleanUtils;
import org.apache.commons.lang.builder.ToStringBuilder;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;

import java.util.Locale;

/**
 * Java/R implementation of ACM Transactions on Mathematical Software Algorithm 643
 * FEXACT.  See <a href="http://portal.acm.org/citation.cfm?id=214326">
 * http://portal.acm.org/citation.cfm?id=214326</a>.
 * <p>
 * This requires native R libraries and rJava to be installed. From the R
 * console enter:
 * <p>&nbsp;&nbsp;<em>install.packages('rJava')</em></p>
 * <p>When running the Java code, you need to add the R and rJava libraries to the library path.
 * For example on windows, add
 * <em>-Djava.library.path="C:\Program Files (x86)\R\R-2.10.1\library\rJava\jri"</em>.  On
 * Unix, add the R and JRI paths to the <em>LD_LIBRARY_PATH</em> environment variable.
 * <p>
 * See <a href="http://www.r-project.org/">The R Project for Statistical Computing</a> and
 * <a href="http://www.rforge.net/rJava/">rJava</a> for reference.
 */
public class FisherExact {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(FisherExact.class);

    /**
     * Error string indicating R was able to execute the fisher test.
     */
    private static final String R_NOT_AVAILABLE = "R is not available.";

    /**
     * Create a new fisher exact test object.
     */
    public FisherExact() {
        super();
    }

    /**
     * Performs Fisher's exact test for testing the null of independence of rows and columns
     * in a contingency table with fixed marginals.
     * <p>
     * The order of the values in the table are "by row" so that if the array contains
     * (1,2,3,11,12,13) and nrows = 2, ncols = 3, the following matrix is created:<br/>
     * <pre>
     *      C1  C2  C3
     *  R1   1   2   3
     *  R2  11  12  13
     * </pre>
     *
     * @param vector An array of integer values used to populate the matrtix to be evaluated.
     * @param nrows The number of rows in the resulting martrix
     * @param ncols The number of columns in the resulting matrix
     * @param hybrid Wheter exact probabilities are computed (false) or a hybrid approximation
     * (true) is made
     * @return The pValue result from the fisher test
     */
    public Result fexact(final int[] vector, final int nrows, final int ncols,
                         final boolean hybrid) {
        assert vector != null : "Input vector cannot be null";
        assert vector.length > 0 : "Input vector cannot be empty";
        assert nrows >= 2 && ncols >= 2 : "Must have at least 2 rows and columns";

        final Result result;
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null) {
            final boolean vectorAssignResult = rengine.assign("vector", vector);
            if (LOG.isDebugEnabled()) {
                LOG.debug("Vector assigned: " +  vectorAssignResult);
                final REXP vectorExpression = rengine.eval("vector");
                LOG.debug("Vector: " + vectorExpression);
            }

            final StringBuilder fisherExpression = new StringBuilder("fisher.test(matrix(vector,");
            fisherExpression.append(nrows);
            fisherExpression.append(',');
            fisherExpression.append(ncols);
            fisherExpression.append(",byrow=TRUE), hybrid=");
            fisherExpression.append(
                    BooleanUtils.toStringTrueFalse(hybrid).toUpperCase(Locale.getDefault()));
            fisherExpression.append(")");
            final boolean is2x2 = nrows == 2 && ncols == 2;
            result = evaluteFisherExpression(rengine, fisherExpression.toString(), is2x2);
        } else {
            LOG.warn(R_NOT_AVAILABLE);
            result = new Result();
        }
        return result;
    }

    /**
     * Performs Fisher's exact test using two input vectors.
     * @param factor1
     * @param factor2
     * @param hybrid
     * @return
     */
    public Result fexact(final int[] factor1, final int[] factor2, final boolean hybrid) {
        assert factor1 != null && factor2 != null : "Input vector cannot be null";
        assert factor1.length == factor2.length : "Length of the two input vectors must be equal";

        final Result result;
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null) {
            final boolean xAssignResult = rengine.assign("x", factor1);
            if (LOG.isDebugEnabled()) {
                LOG.debug("X assigned: " +  xAssignResult);
                final REXP xExpression = rengine.eval("x");
                LOG.debug("X: " + xExpression);
            }

            final boolean yAssignResult = rengine.assign("y", factor2);
            if (LOG.isDebugEnabled()) {
                LOG.debug("Y assigned: " +  yAssignResult);
                final REXP yExpression = rengine.eval("y");
                LOG.debug("Y: " + yExpression);
            }

            final StringBuilder fisherExpression = new StringBuilder("fisher.test(x, y, hybrid=");
            fisherExpression.append(
                    BooleanUtils.toStringTrueFalse(hybrid).toUpperCase(Locale.getDefault()));
            fisherExpression.append(')');
            final boolean is2x2 = factor1.length == 2 && factor2.length == 2;
            result = evaluteFisherExpression(rengine, fisherExpression.toString(), is2x2);
        } else {
            LOG.warn(R_NOT_AVAILABLE);
            result = new Result();
        }
        return result;
    }

    /**
     * Pass the fisher expression to R for computation.
     * @param rengine The R engine to use to calcuate the results.
     * @param fisherExpression The string representing the expression to evaluate.
     * @param is2x2 Whether or not the data being evaluated represents a 2x2 matrix
     * @return The results of the evaluation (should never be null)
     */
    private Result evaluteFisherExpression(final Rengine rengine,
                                           final String fisherExpression,
                                           final boolean is2x2) {
        final Result result;

        if (LOG.isDebugEnabled()) {
            LOG.debug("About to evaluate: " + fisherExpression);
        }
        final REXP fisherResultExpression = rengine.eval(fisherExpression);
        LOG.debug(fisherResultExpression);

        final RVector fisherResultVector = fisherResultExpression.asVector();
        LOG.debug(fisherResultVector);

        final REXP pValueExpression = fisherResultVector.at("p.value");
        final double pValue;
        if (pValueExpression != null) {
            pValue = pValueExpression.asDouble();
        } else {
            pValue = Double.NaN;
        }

        final REXP alternativeExpression = fisherResultVector.at("alternative");
        final String alternative = alternativeExpression.asString();
        LOG.debug("alternative: " + alternative);

        if (is2x2) {
            final REXP estimateExpression = fisherResultVector.at("estimate");
            LOG.debug(estimateExpression);
            final double estimate;
            if (estimateExpression != null) {
                estimate = estimateExpression.asDouble();
            } else {
                estimate = Double.NaN;
            }
            LOG.debug("estimate: " + estimate);

            final REXP confidenceIntervalExpression = fisherResultVector.at("conf.int");
            LOG.debug(confidenceIntervalExpression);
            final double[] confidenceInterval;
            if (confidenceIntervalExpression != null) {
                confidenceInterval = confidenceIntervalExpression.asDoubleArray();
            } else {
                confidenceInterval = ArrayUtils.EMPTY_DOUBLE_ARRAY;
            }
            LOG.debug("confidenceInterval: " + ArrayUtils.toString(confidenceInterval));

            final REXP oddsRatioExpression = fisherResultVector.at("null.value");
            LOG.debug(oddsRatioExpression);
            final double oddsRatio;
            if (oddsRatioExpression != null) {
                oddsRatio = oddsRatioExpression.asDouble();
            } else {
                oddsRatio = Double.NaN;
            }
            LOG.debug("oddsRatio: " + ArrayUtils.toString(oddsRatio));

            result = new Result(pValue, confidenceInterval, estimate,
                    oddsRatio, AlternativeHypothesis.twosided);
        } else {
            result = new Result(pValue);
        }
        return result;
    }

    /**
     * Performs Fisher's exact test for testing the null of independence of rows and columns
     * in a contingency table with fixed marginals.
     * <p>
     * The order of the values in the table are "by row" so that if the array contains
     * (1,2,3,11,12,13) and nrows = 2, ncols = 3, the following matrix is created:<br/>
     * <pre>
     *      C1  C2  C3
     *  R1   1   2   3
     *  R2  11  12  13
     * </pre>
     *
     * @param vector An array of integer values used to populate the matrtix to be evaluated.
     * @param nrows The number of rows in the resulting martrix
     * @param ncols The number of columns in the resulting matrix
     * @return The pValue result from the fisher test
     */
    public Result fexact(final int[] vector, final int nrows, final int ncols) {
        return fexact(vector, nrows, ncols, false);
    }

    /**
     * Performs Fisher's exact test for testing the null of independence of rows and columns
     * in a contingency table with fixed marginals.
     * <p>
     * The paramters to this method create a matrix of the following form.
     * <pre>
     *        C1    C2
     *  R1  r1c1  r1c2
     *  R2  r2c1  r2c2
     * </pre>
     *
     * @param r1c1 Value for row1/column1 of the contingency matrix
     * @param r1c2 Value for row1/column2 of the contingency matrix
     * @param r2c1 Value for row1/column1 of the contingency matrix
     * @param r2c2 Value for row2/column2 of the contingency matrix
     * @return The pValue result from the fisher test
     */
    public Result fexact(final int r1c1, final int r1c2, final int r2c1, final int r2c2) {
        return fexact(new int[] {r1c1, r1c2, r2c1, r2c2}, 2, 2);
    }

    /**
     *
     */
    public enum AlternativeHypothesis {
        /**
         * Two-sided tests are based on the probabilities of the tables, and take as "more extreme"
         * all tables with probabilities less than or equal to that of the observed table, the
         * p-value being the sum of such probabilities.
         */
        twosided,
        greater,
        less
    }

    /**
     * Result from the fisher exact test.  Note that some fields may only considered
     * valid when the input was a 2 by 2 contingency matrix.
     */
    public class Result {
        /**
         * Indicates that R did in fact return the results contained in this object.
         */
        private final boolean valid;
        /**
         * The p-value returned by the test.
         */
        private final double pValue;
        /**
         * A confidence interval for the odds ratio. Only present in the 2 by 2 case.
         */
        private final double[] confidenceInterval;
        /**
         * An estimate of the odds ratio. Note that the conditional Maximum Likelihood
         * Estimate (MLE) rather than the unconditional MLE (the sample odds ratio) is used.
         * Only present in the 2 by 2 case.
         */
        private final double estimate;
        /**
         * The odds ratio. Only present in the 2 by 2 case
         */
        private final double oddsRatio;
        /**
         * Indicates the alternative hypothesis used to compute the result.
         * Only used in the 2 by 2 case.
         */
        private final AlternativeHypothesis alternativeHypothesis;

        /**
         * Used to create an invalid result.
         */
        private Result() {
            super();
            this.pValue = Double.NaN;
            this.confidenceInterval = ArrayUtils.EMPTY_DOUBLE_ARRAY;
            this.estimate = Double.NaN;
            this.oddsRatio = Double.NaN;
            this.alternativeHypothesis = AlternativeHypothesis.twosided;

            valid = false;
        }

        /**
         * Used to create a result containing only pValues.  All other fields are
         * set to default values.
         * @param pValue The pValue of the result.
         */
        public Result(final double pValue) {
            this(pValue, ArrayUtils.EMPTY_DOUBLE_ARRAY, Double.NaN, Double.NaN,
                    AlternativeHypothesis.twosided);
        }

        /**
         * Used to create a result containing values returned by the R fisher exact method.
         * @param pValue The pValue of the result.
         * @param confidenceInterval The confidence interval for the odds ratio.
         * @param estimate The estimate of the odds ratio
         * @param oddsRatio The odds ratio
         * @param alternativeHypothesis The alternative hypothesis used to compute the result
         */
        public Result(final double pValue, final double[] confidenceInterval,
                      final double estimate, final double oddsRatio,
                      final AlternativeHypothesis alternativeHypothesis) {
            super();
            this.pValue = pValue;
            this.confidenceInterval = confidenceInterval;
            this.estimate = estimate;
            this.oddsRatio = oddsRatio;
            this.alternativeHypothesis = alternativeHypothesis;

            valid = true;
        }

        /**
         * Indicates that R did in fact return the results contained in this object.
         * @return whether or not values were actually computed or an error occurred
         */
        public boolean isValid() {
            return valid;
        }

        /**
         * Get the p-value returned by the test or {@link Double#NaN} if R was not avaialble.
         * @return The p-value
         */
        public double getPValue() {
            return pValue;
        }

        /**
         * Get a confidence interval for the odds ratio. Will be an empty array if the input
         * was not a 2 by 2 matrix or R was not avaialble.
         * @return A confidence interval for the odds ratio
         */
        public double[] getConfidenceInterval() {
            return confidenceInterval;
        }

        /**
         * Get an estimate of the odds ratio.  Will be {@link Double#NaN} if the input
         * was not a 2 by 2 matrix or R was not avaialble.
         * @return An estimate of the odds ratio
         */
        public double getEstimate() {
            return estimate;
        }

        /**
         * Get the odds ratio. Will be {@link Double#NaN} if the input
         * was not a 2 by 2 matrix or R was not avaialble.
         * @return The odds ratio
         */
        public double getOddsRatio() {
            return oddsRatio;
        }

        /**
         * Get the alternative hypothesis used to compute the result.  Will be
         * {@link edu.cornell.med.icb.goby.R.FisherExact.AlternativeHypothesis#twosided} if
         * the input was not a 2 by 2 matrix or R was not avaialble.
         * @return the alternative hypothesis used to compute the result
         */
        public AlternativeHypothesis getAlternativeHypothesis() {
            return alternativeHypothesis;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public String toString() {
            return new ToStringBuilder(this)
                    .append("pValue", pValue)
                    .append("confidenceInterval", confidenceInterval)
                    .append("estimate", estimate)
                    .append("oddsRatio", oddsRatio)
                    .append("alternativeHypothesis", alternativeHypothesis)
                    .append("valid", valid)
                    .toString();
        }
    }

    public static void main(final String[] args) {
        // TODO - convert to unit tests
        final FisherExact fisherExact = new FisherExact();

        Result result;
        result = fisherExact.fexact(new int[] { 12, 17, 4, 25, 15, 4 }, 3, 2, true);
        System.out.println("pValue: " + result);

        result = fisherExact.fexact(new int[] { 1, 2, 3, 11, 12, 13}, 2, 3);
        System.out.println("pValue: " + result);


        result = fisherExact.fexact(new int[] { 22, 13, 5,  4,  5,  3,  2,  1, 7,  1,  4,  3,  1,  2,  3,  4 }, 8, 2);
        System.out.println("pValue: " + result);

        result = fisherExact.fexact(new int[] { 5, 1, 1, 3 }, 2, 2);
        System.out.println("pValue: " + result);

        result = fisherExact.fexact(625, 256, 81, 16);
        System.out.println("pValue: " + result);

        System.out.println(" ========= gominer ========");

        final Fisher fisher = new Fisher();
        final double pValue = fisher.calculateFisherFromMatrix(625, 256, 81, 16);
        System.out.println("pValue: " + pValue);
        System.out.println("  left: " + fisher.getLeft());
        System.out.println(" right: " + fisher.getRight());
        System.out.println("2-tail: " + fisher.getTwotail());

        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null) {
            rengine.end();
        }
    }
}
