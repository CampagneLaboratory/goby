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

package edu.cornell.med.icb.goby.algorithmic.data;

import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCountInterface;
import edu.cornell.med.icb.goby.algorithmic.algorithm.ComputeCountInterface;
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationWeightCount;
import edu.cornell.med.icb.goby.algorithmic.algorithm.AnnotationCount;

/**
 * @author Fabien Campagne
 *         Date: May 25, 2010
 *         Time: 5:43:34 PM
 */
public class FormulaWeightCount implements AnnotationCountInterface {
    AnnotationWeightCount weightCounter;
    AnnotationCount regularCounter;

    public FormulaWeightCount(WeightsInfo weights) {
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

    public void populate(int startPosition, int endPosition, int queryIndex) {
        weightCounter.populate(startPosition, endPosition, queryIndex);
        regularCounter.populate(startPosition, endPosition, queryIndex);
    }

    public void sortReads() {
        weightCounter.sortReads();
        regularCounter.sortReads();
    }

    public float averageReadsPerPosition(int geneStart, int geneEnd) {
        return (float) evaluateFormula(weightCounter.averageReadsPerPosition(geneStart, geneEnd),
                regularCounter.averageReadsPerPosition(geneStart, geneEnd));
    }


    public double countReadsPartiallyOverlappingWithInterval(int geneStart, int geneEnd) {
        return evaluateFormula(weightCounter.countReadsPartiallyOverlappingWithInterval(geneStart, geneEnd),
                regularCounter.countReadsPartiallyOverlappingWithInterval(geneStart, geneEnd));
    }

    public double countReadsStriclyWithinInterval(int geneStart, int geneEnd) {
        return evaluateFormula(weightCounter.countReadsStriclyWithinInterval(geneStart, geneEnd),
                regularCounter.countReadsStriclyWithinInterval(geneStart, geneEnd));
    }

    enum FormulaChoice {

        FORMULA1,
        FORMULA2,
        FORMULA3,


    }

    private static FormulaChoice formulaChoice = FormulaChoice.FORMULA2;

    private double evaluateFormula(double sumGamma, double rawCount) {
        double value;
        switch (formulaChoice) {
            case FORMULA1:
                value = (float) Math.exp(-0.898699452975139d +
                        4.32171188401506d * Math.log(rawCount)
                        - 3.42375631012767d * Math.log(sumGamma));
                return value;
            case FORMULA2: {
                double logGC_a = Math.log(sumGamma) - Math.log(rawCount);
                value = (float) Math.exp(-1.57361718748031d - 3.62887641327823 * logGC_a + Math.log(rawCount));
                return value;
            }
            case FORMULA3: {
                double logGC_a_plus1 = Math.log(sumGamma + 1) - Math.log(rawCount + 1);
                value = Math.exp(-0.843924877396631d + 0.887303234304011d * Math.log(rawCount + 1) -
                        3.45874660923795d * logGC_a_plus1);
                return value;
            }
            default:
                return Double.NaN;
        }

    }

    public double geneExpressionCount(Annotation annot) {
        double weightExpression = weightCounter.geneExpressionCount(annot);
        double regularExpression = regularCounter.geneExpressionCount(annot);
        return evaluateFormula(weightExpression, regularExpression);
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
