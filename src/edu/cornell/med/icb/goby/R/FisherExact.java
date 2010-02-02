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
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;

/**
 * Java/R implementation of ACM Transactions on Mathematical Software Algorithm 643
 * FEXACT.  See <a href="http://portal.acm.org/citation.cfm?id=214326">
 * http://portal.acm.org/citation.cfm?id=214326</a>.  This requires native R libraries
 * and rJava to be installed. From the R console <em>install.packages('rJava')</em>.
 *
 * See <a href="http://www.r-project.org/">The R Project for Statistical Computing</a> and
 * <a href="http://www.rforge.net/rJava/">rJava</a> for reference.
 */
public class FisherExact {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(FisherExact.class);

    public FisherExact() {
        super();
    }

    public double fexact(final int[] rows, final int[] cols) {
        return fexact(ArrayUtils.addAll(rows, cols), rows.length, cols.length);
    }

    public double fexact(final int[] vector, final int nrows, final int ncols) {
        final double pValue;
        final Rengine rengine = GobyRengine.getInstance().getRengine();
        if (rengine != null) {
            final boolean vectorAssignResult = rengine.assign("vector", vector);
            LOG.debug("assign: " +  vectorAssignResult);

            final REXP vectorExpression = rengine.eval("vector");
            LOG.debug(vectorExpression);

            final StringBuilder fisherExpression = new StringBuilder("fisher.test(matrix(vector,");
            fisherExpression.append(nrows);
            fisherExpression.append(',');
            fisherExpression.append(ncols);
            fisherExpression.append("))");

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

    public double fexact(final int n11, final int n12, final int n21, final int n22) {
        return fexact(new int[] {n11, n12, n21, n22}, 2, 2);
    }

    public static void main(final String[] args) throws InterruptedException {
        final FisherExact fisherExact = new FisherExact();

        // example comes from http://darwin.eeb.uconn.edu/eeb348/supplements-2006/chi-squared/chi-squared.html
        double pValue = fisherExact.fexact(new int[] { 12, 4, 15, 17, 25, 4 }, 3, 2);
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
