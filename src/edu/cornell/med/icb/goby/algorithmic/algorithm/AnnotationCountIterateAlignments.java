/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.modes.WeightParameters;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.io.IOException;

/**
 * Iterate through an alignment to populate an array of AnnotationCountInterface with read matches.
 *
 * @author Fabien Campagne
 *         Date: Jun 23, 2010
 *         Time: 2:21:19 PM
 */
public class AnnotationCountIterateAlignments extends IterateAlignments {
    private WeightParameters weightParams;
    private WeightsInfo weights;

    /**
     * Retrieves the populated instances of AnnotationCountInterface.
     *
     * @return
     */
    public AnnotationCountInterface[] getAlgs() {
        return algs;
    }

    private int numAlignedReadsInSample;

    public void setWeightInfo(final WeightParameters weightParams, final WeightsInfo weights) {
        this.weightParams = weightParams;
        this.weights = weights;

    }

    @Override
    public void processNumberOfReferences(final String basename, final int numberOfReferences) throws IOException {
        algs = new AnnotationCountInterface[numberOfReferences];
    }

    private AnnotationCountInterface[] algs;

    @Override
    public void processAlignmentEntry(final AlignmentReaderImpl alignmentReader,
                                      final Alignments.AlignmentEntry alignmentEntry) {
        final int referenceIndex = alignmentEntry.getTargetIndex();
        final int startPosition = alignmentEntry.getPosition();

        final int alignmentLength = alignmentEntry.getQueryAlignedLength();
        //shifted the ends populating by 1
        for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {
            algs[referenceIndex].populate(startPosition, startPosition + alignmentLength, alignmentEntry.getQueryIndex());
            ++numAlignedReadsInSample;

        }
    }


    @Override
    public void prepareDataStructuresForReference(final AlignmentReaderImpl alignmentReader, final int referenceIndex) {
        AnnotationCountInterface algo = new AnnotationCount();

        algo = chooseAlgorithm(weightParams, weights, algo);
        algs[referenceIndex] = algo;
        algs[referenceIndex].startPopulating();
        referencesSelected.add(referenceIndex);
    }


    public int getNumAlignedReadsInSample() {
        return numAlignedReadsInSample;
    }
    // determine which algorithm to use based on weight parameters.
    private AnnotationCountInterface chooseAlgorithm(final WeightParameters params,
                                                     final WeightsInfo weights,
                                                     AnnotationCountInterface algo) {
        if (params.useWeights) {
            if (!params.adjustGcBias) {

                algo = new AnnotationWeightCount(weights);
            } else {

                final FormulaWeightAnnotationCount algo1 = new FormulaWeightAnnotationCount(weights);

                algo1.setFormulaChoice(FormulaWeightAnnotationCount.FormulaChoice.valueOf(params.formulaChoice));
                algo = algo1;
            }
        }
        return algo;
    }
    private IntSet referencesSelected = new IntOpenHashSet();

     public IntSet getReferencesSelected() {
         return referencesSelected;
     }
}
