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

import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Map;

/**
 * Calculates reweighted counts for .counts output.
 * @author Fabien Campagne
 *         Date: May 25, 2010
 *         Time: 5:43:34 PM
 */
public class FormulaWeightCount implements ComputeCountInterface {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(FormulaWeightCount.class);

    /**
     * Counter that sums read weights.
     */
    ComputeWeightCount weightCounter;
    /**
     * counter that produces raw counts.
     */
    ComputeCount regularCounter;

    public FormulaWeightCount(final WeightsInfo weights) {
        this.weightCounter = new ComputeWeightCount(weights);
        this.regularCounter = new ComputeCount();
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
        regularCounter.populate(startPosition, endPosition);
    }


    public void setFormulaChoice(final FormulaWeightAnnotationCount.FormulaChoice formulaChoice) {
        this.formulaChoice = formulaChoice;
    }

    private FormulaWeightAnnotationCount.FormulaChoice formulaChoice = FormulaWeightAnnotationCount.FormulaChoice.FORMULA2;


    public void accumulate() {
        weightCounter.accumulate();
        regularCounter.accumulate();
    }

    public void baseCount(final CountsWriter writer) throws IOException {
        final IntSortedSet joints = new IntAVLTreeSet();
        joints.addAll(regularCounter.starts.keySet());
        joints.addAll(regularCounter.starts.keySet());
        joints.addAll(weightCounter.starts.keySet());
        joints.addAll(weightCounter.ends.keySet());

        final int[] jointsArray = joints.toArray(new int[joints.size()]);
        double prevCount = 0;
        int lengthConstant = 0;
        double count;
        int line = 0;
        LOG.debug("counting");
        int rawStartValue = regularCounter.starts.get(0);
        int rawEndValue = regularCounter.ends.get(0);
        double weightStartValue = weightCounter.starts.get(0);
        double weightEndValue = weightCounter.ends.get(0);

        for (int i = 1; i < jointsArray.length; i++) {
            if (line % 1000000 == 0) {
                LOG.debug("line " + line);
            }
            line++;
            final int curKey = jointsArray[i];
            final int prevKey = jointsArray[i - 1];

            if (regularCounter.starts.containsKey(curKey)) {
                rawStartValue = regularCounter.starts.get(curKey);
            }
            if (weightCounter.starts.containsKey(curKey)) {
                weightStartValue = weightCounter.starts.get(curKey);
            }
            if (regularCounter.ends.containsKey(curKey)) {
                rawEndValue = regularCounter.ends.get(curKey);
            }
            if (weightCounter.ends.containsKey(curKey)) {
                weightEndValue = weightCounter.ends.get(curKey);
            }
            final double sumGamma = weightStartValue - weightEndValue;
            final int rawCount = rawStartValue - rawEndValue;

            count = FormulaWeightAnnotationCount.evaluateFormula(formulaChoice,
                    sumGamma,
                    rawCount);
            //     System.out.printf("sumGamma: %g rawCount %d reweighted: %g%n", sumGamma, rawCount, count);

            lengthConstant += curKey - prevKey;

            if (Math.abs(count - prevCount) > 1) {

                writer.appendCount((int) prevCount, lengthConstant);
                prevCount = count;
                lengthConstant = 0;
            }
        }
        writer.close();

    }

    public void baseCount() {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    public Map getCountPerBase() {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    public IntList getCountKeys() {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    public void populate(final int startPosition, final int endPosition, final boolean forwardStrand, final int queryIndex) {
        weightCounter.populate(startPosition, endPosition, forwardStrand, queryIndex);
        regularCounter.populate(startPosition, endPosition, forwardStrand, queryIndex);
    }

    public void populate(final int startPosition, final int endPosition, final boolean forwardStrand) {
        weightCounter.populate(startPosition, endPosition, forwardStrand);
        regularCounter.populate(startPosition, endPosition, forwardStrand);
    }

}
