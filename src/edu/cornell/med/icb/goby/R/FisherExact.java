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
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;

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
    public double fexact(final int[] vector, final int nrows, final int ncols, final boolean hybrid) {
        final double pValue;
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null) {
            final boolean vectorAssignResult = rengine.assign("vector", vector);
            LOG.debug("Vector assigned: " +  vectorAssignResult);

            final REXP vectorExpression = rengine.eval("vector");
            LOG.debug("Vector: " + vectorExpression);

            final StringBuilder fisherExpression = new StringBuilder("fisher.test(matrix(vector,");
            fisherExpression.append(nrows);
            fisherExpression.append(',');
            fisherExpression.append(ncols);
            fisherExpression.append(",byrow=TRUE), hybrid=");
            fisherExpression.append(BooleanUtils.toStringTrueFalse(hybrid).toUpperCase());
            fisherExpression.append(")");

            final REXP fisherResultExpression = rengine.eval(fisherExpression.toString());
            LOG.debug(fisherResultExpression);

            final RVector fisherResultVector = fisherResultExpression.asVector();
            LOG.debug(fisherResultVector);

            final REXP pValueExpression = fisherResultVector.at("p.value");
            if (pValueExpression != null) {
                pValue = pValueExpression.asDouble();
            } else {
                pValue = Double.NaN;
            }

            if (nrows == 2 && ncols == 2) {
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
                    confidenceInterval = null;
                }
                LOG.debug("confidenceInterval: " + ArrayUtils.toString(confidenceInterval));

                final REXP nullValueExpression = fisherResultVector.at("null.value");
                LOG.debug(nullValueExpression);
                final double[] nullValue;
                if (nullValueExpression != null) {
                    nullValue = nullValueExpression.asDoubleArray();
                } else {
                    nullValue = null;
                }
               LOG.debug("nullValue: " + ArrayUtils.toString(nullValue));

            }
        } else {
            LOG.warn("R is not available.  Returning invalid pValue");
            pValue = Double.NaN;
        }
        return pValue;
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
    public double fexact(final int[] vector, final int nrows, final int ncols) {
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
    public double fexact(final int r1c1, final int r1c2, final int r2c1, final int r2c2) {
        return fexact(new int[] {r1c1, r1c2, r2c1, r2c2}, 2, 2);
    }

    public static void main(final String[] args) {
        // TODO - convert to unit tests
        final FisherExact fisherExact = new FisherExact();

        // example comes from http://darwin.eeb.uconn.edu/eeb348/supplements-2006/chi-squared/chi-squared.html
        double pValue = fisherExact.fexact(new int[] { 12, 17, 4, 25, 15, 4 }, 3, 2);
        System.out.println("pValue: " + pValue);

        pValue = fisherExact.fexact(new int[] { 12, 17, 4, 25, 15, 4 }, 3, 2, true);
        System.out.println("pValue: " + pValue);

        pValue = fisherExact.fexact(new int[] { 1, 2, 3, 11, 12, 13}, 2, 3);
        System.out.println("pValue: " + pValue);

        pValue = fisherExact.fexact(10 , 20 , 30 , 40);
        System.out.println("pValue: " + pValue);

        pValue = fisherExact.fexact(new int[] { 22, 13, 5,  4,  5,  3,  2,  1, 7,  1,  4,  3,  1,  2,  3,  4 }, 8, 2);
        System.out.println("pValue: " + pValue);

        pValue = fisherExact.fexact(new int[] { 5, 1, 1, 3 }, 2, 2);
        System.out.println("pValue: " + pValue);

        pValue = fisherExact.fexact(625, 256, 81, 16);
        System.out.println("pValue: " + pValue);

        System.out.println(" ========= gominer ========");

        final Fisher fisher = new Fisher();
        pValue = fisher.calculateFisherFromMatrix(625, 256, 81, 16);
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
