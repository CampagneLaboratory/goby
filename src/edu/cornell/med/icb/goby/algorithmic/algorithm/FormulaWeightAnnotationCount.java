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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;

/**
 * @author Fabien Campagne
 *         Date: May 25, 2010
 *         Time: 5:43:34 PM
 */
public class FormulaWeightAnnotationCount implements AnnotationCountInterface {
    AnnotationWeightCount weightCounter;
    AnnotationCount regularCounter;

    public FormulaWeightAnnotationCount(final WeightsInfo weights) {
        this.weightCounter = new AnnotationWeightCount(weights);
        this.regularCounter = new AnnotationCount();
    }


    public ComputeCountInterface getBaseCounter() {
        throw new UnsupportedOperationException("This operation is currently not implemented.");
    }

    public void startPopulating() {
        weightCounter.startPopulating();
        regularCounter.startPopulating();
    }

    public void populate(final int startPosition, final int endPosition, final int queryIndex) {
        weightCounter.populate(startPosition, endPosition, queryIndex);
        regularCounter.populate(startPosition, endPosition, queryIndex);
    }

    public void sortReads() {
        weightCounter.sortReads();
        regularCounter.sortReads();
    }

    public float averageReadsPerPosition(final int geneStart, final int geneEnd) {
        return (float) evaluateFormula(formulaChoice, weightCounter.averageReadsPerPosition(geneStart, geneEnd),
                regularCounter.averageReadsPerPosition(geneStart, geneEnd));
    }


    public double countReadsPartiallyOverlappingWithInterval(final int geneStart, final int geneEnd) {
        return evaluateFormula(formulaChoice, weightCounter.countReadsPartiallyOverlappingWithInterval(geneStart, geneEnd),
                regularCounter.countReadsPartiallyOverlappingWithInterval(geneStart, geneEnd));
    }

    public double countReadsStriclyWithinInterval(final int geneStart, final int geneEnd) {
        return evaluateFormula(formulaChoice, weightCounter.countReadsStriclyWithinInterval(geneStart, geneEnd),
                regularCounter.countReadsStriclyWithinInterval(geneStart, geneEnd));
    }

    public enum FormulaChoice {

        FORMULA1,
        FORMULA2,
        FORMULA3,
        /**
         * Fit against both Helicos and SOLID on HBR datasets.
         */
        FORMULA4

    }

    public void setFormulaChoice(final FormulaChoice formulaChoice) {
        this.formulaChoice = formulaChoice;
    }

    private FormulaChoice formulaChoice = FormulaChoice.FORMULA2;

    public static double evaluateFormula(final FormulaChoice formulaChoice, final double sumGamma, final double rawCount) {
        if (rawCount == 0) {
            return 0;
        }
        if (sumGamma == 0) {
            return rawCount;
        }
        double value;
        switch (formulaChoice) {
            case FORMULA1:
                value = (float) Math.exp(-0.898699452975139d +
                        4.32171188401506d * Math.log(rawCount)
                        - 3.42375631012767d * Math.log(sumGamma));
                value = Math.max(value, 0);
                return value;
            case FORMULA2: {
                final double logGC_a = Math.log(sumGamma) - Math.log(rawCount);
                value = (float) Math.exp(-1.57361718748031d - 3.62887641327823 * logGC_a + Math.log(rawCount));
                value = Math.max(value, 0);
                return value;
            }
            case FORMULA3: {
                final double logGC_a_plus1 = Math.log(sumGamma + 1) - Math.log(rawCount + 1);
                value = Math.exp(-0.843924877396631d + 0.887303234304011d * Math.log(rawCount + 1) -
                        3.45874660923795d * logGC_a_plus1);
                value = Math.max(value, 0);
                return value;
            }
            case FORMULA4: {
                // These estimates were obtained by comparing the Bullard Illumina HBR dataset to the SEQC Helicos
                // and SOLID datasets. A covariate was used to represent the Helicos or SOLID target platform.

                final double logGC_a = Math.log(sumGamma) - Math.log(rawCount);
                value = (float) Math.exp(-1.4050204825287 - 3.5820783386146 * logGC_a + Math.log(rawCount));
                value = Math.max(value, 0);
                return value;
            }
            default:
                return Double.NaN;
        }

    }

    public double geneExpressionCount(final Annotation annot) {
        final double weightExpression = weightCounter.geneExpressionCount(annot);
        final double regularExpression = regularCounter.geneExpressionCount(annot);
        return evaluateFormula(formulaChoice, weightExpression, regularExpression);
    }

    public void accumulate() {
        weightCounter.accumulate();
        regularCounter.accumulate();
    }

    public void baseCount() {
        weightCounter.baseCount();
        regularCounter.baseCount();
    }

}
