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

package edu.cornell.med.icb.goby.stats;

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2DoubleArrayMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;

import java.util.Map;
import java.util.Set;

/**
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 3:35:24 PM
 */
public class DifferentialExpressionCalculator {
    private final Set<String> groups;
    private final Map<String, String> sampleToGroupMap;
    private final IndexedIdentifier elementLabels;
    private int elementsPerSample;
    private int numberOfSamples;
    private final Map<String, IntArrayList> sampleToCounts;
    private Object2DoubleMap<String> sampleProportions;
    private final IntArrayList lengths;
    private final Map<MutableString, MutableString> elementLabelToElementType;
    /**
     * The number of alignment entries observed in each sample.
     */
    private final Object2IntMap<String> numAlignedInSample;

    public DifferentialExpressionCalculator() {
        super();
        groups = new ObjectArraySet<String>();
        elementLabels = new IndexedIdentifier();
        sampleToGroupMap = new Object2ObjectOpenHashMap<String, String>();
        numAlignedInSample = new Object2IntOpenHashMap<String>();
        sampleToCounts = new Object2ObjectOpenHashMap<String, IntArrayList>();
        lengths = new IntArrayList();
        elementLabelToElementType = new Object2ObjectOpenHashMap<MutableString, MutableString>();
    }

    public double calculateNormalized(final int readCountInt, final int annotLength, final double normalizationFactor) {
        final double readCount = readCountInt;
        final double length = annotLength; // in bases
        final double sampleReadCount = normalizationFactor; // in reads
        return readCount / (length / 1000.0d) / (normalizationFactor / 1E6d);
    }

    public synchronized void defineGroup(final String label) {
        groups.add(label);
    }

    /**
     * Define an element but don't specify the type. The version with the type is the preferred version.
     * @param label element-id
     * @return the index of the element
     */
    public synchronized int defineElement(final String label) {
        return defineElement(label, "");
    }

    /**
     * Define an element and it's type
     * @param label element-id
     * @param type the type (gene, exon, ...)
     * @return the index of the element
     */
    public synchronized int defineElement(final String label, final String type) {
        final MutableString elementLabel = new MutableString(label);
        if (elementLabelToElementType.get(elementLabel) == null) {
            elementLabelToElementType.put(elementLabel, new MutableString(type));
        }
        return elementLabels.registerIdentifier(elementLabel);
    }

    public Map<MutableString, MutableString> getElementLabelToElementTypeMap() {
        return elementLabelToElementType;
    }

    /**
     * Define the number of sequence bases that this element spans. The bases do not need to be contiguous (i.e., multiple
     * exons of a transcript).
     *
     * @param elementIndex index of the element associated with length.
     * @param length       length of the element.
     */
    public synchronized void defineElementLength(final int elementIndex, final int length) {
        if (lengths.size() <= elementIndex) {
            lengths.size(elementIndex + 1);
        }
        lengths.set(elementIndex, length);
    }

    public synchronized void associateSampleToGroup(final String sample, final String group) {
        sampleToGroupMap.put(sample, group);
    }

    /**
     * Return the length of an element.
     *
     * @param elementId
     * @return
     */
    public int getElementLength(final MutableString elementId) {
        final int elementIndex = elementLabels.getInt(elementId);
        return lengths.get(elementIndex);
    }

    /**
     * Observe counts and RPKM for a sample.
     *
     * @param sample    sample id.
     * @param elementId element id.
     * @param count     Number of reads that can be assigned to the element.
     */
    public void observe(final String sample, final String elementId, final int count) {


        IntArrayList counts = sampleToCounts.get(sample);
        // the following looks a bit complicated. We are trying to avoid synchronizing every time observe is called.
        // This would slow the whole process too much. Instead, we synchronize only when we need to create a counts
        // data structure for a new sample, which should not happen too often.
        if (counts == null) {
            synchronized (this) {
                counts = sampleToCounts.get(sample);
                if (counts == null) {
                    counts = new IntArrayList(elementsPerSample);
                    counts.size(elementsPerSample);
                    sampleToCounts.put(sample, counts);
                }
            }
        }


        final int elementIndex = elementLabels.get(new MutableString(elementId));
        counts.set(elementIndex, count);
    }

    /**
     * Define the number of alignment entries found in each sample.
     *
     * @param sampleId            The sample
     * @param numAlignedInSamples The number of alignment entries observed in the sample.
     */
    public synchronized void setNumAlignedInSample(final String sampleId, final int numAlignedInSamples) {
        numAlignedInSample.put(sampleId, numAlignedInSamples);
    }

    /**
     * Return the number of alignment entries found in the sample.
     *
     * @param sampleId Identifier of the sample.
     * @return the number of alignment entries found in the sample.
     */
    public int getNumAlignedInSample(final String sampleId) {
        return numAlignedInSample.get(sampleId);
    }


    public DifferentialExpressionResults compare(DifferentialExpressionResults results,
                                                 final NormalizationMethod method,
                                                 final StatisticCalculator tester,
                                                 final String... group) {
        if (results == null) {
            results = new DifferentialExpressionResults();
        }
        if (tester.installed()) {
            assert !tester.canDo(group) : "The number of groups to compare is not supported by the specified calculator.";
            tester.setResults(results);
            return tester.evaluate(this, method, results, group);
        } else {
            return results;
        }
    }

    public DifferentialExpressionResults compare(final StatisticCalculator tester,
                                                 final NormalizationMethod method,
                                                 final String... group) {
        final DifferentialExpressionResults results = new DifferentialExpressionResults();
        if (tester.installed()) {
            assert !tester.canDo(group) : "The number of groups to compare is not supported by the specified calculator.";
            tester.setResults(results);
            return tester.evaluate(this, method, results, group);
        } else {
            return results;
        }
    }

    /**
     * Reserve storage for specified number of elements and samples. One DE test will be conducted for each element.
     *
     * @param elementsPerSample
     * @param numberOfSamples
     */
    public void reserve(final int elementsPerSample, final int numberOfSamples) {
        this.elementsPerSample = elementsPerSample;
        this.numberOfSamples = numberOfSamples;
        lengths.size(elementsPerSample);
    }

    public ObjectArraySet<String> getSamples(final String groupA) {
        final ObjectArraySet<String> samples = new ObjectArraySet<String>();
        for (final String sampleId : sampleToGroupMap.keySet()) {
            if (sampleToGroupMap.get(sampleId).equals(groupA)) {
                samples.add(sampleId);
            }
        }
        return samples;
    }

    public ObjectSet<MutableString> getElementIds() {
        return elementLabels.keySet();
    }

    /**
     * Get the normalized expression value an element in a given sample.  Equivalent to calling the normalizationMethod
     * getNormalizedExpressionValue method.
     *
     * @param sampleId            sample identifier
     * @param normalizationMethod Normalization method (adjusts the denominator of the RPKM value).
     * @param elementId           the element for which a normalized expression value is sought.
     * @return normalized expression value scaled by length and global normalization method.
     */
    public double getNormalizedExpressionValue(final String sampleId, final NormalizationMethod normalizationMethod, final MutableString elementId) {
        return normalizationMethod.getNormalizedExpressionValue(this, sampleId, elementId);

    }

    /**
     * Get the stored overlap count for an element in a given sample.
     *
     * @param sample
     * @param elementId
     * @return
     */
    public int getOverlapCount(final String sample, final MutableString elementId) {
        final IntArrayList counts = sampleToCounts.get(sample);
        if (counts == null) {
            return 0;
        }
        final int elementIndex = elementLabels.get(elementId);

        return counts.get(elementIndex);
    }

    /**
     * Returns the sum of counts in a given sample.
     */
    public int getSumOverlapCounts(final String sample) {
        int sumCounts = 0;
        final IntArrayList counts = sampleToCounts.get(sample);
        if (counts == null) {
            return 0;
        }
        for (final int count : counts) {
            sumCounts += count;
        }
        return sumCounts;
    }

    public String[] samples() {
        final Set<String> sampleSet = sampleToGroupMap.keySet();
        return sampleSet.toArray(new String[sampleSet.size()]);
    }

    /**
     * Returns the proportion of counts that originate from a certain sample. The number of counts in the sample
     * is divided by the sum of counts over all the samples in the experiment.
     *
     * @param sample
     * @return
     */
    public double getSampleProportion(final String sample) {
        if (sampleProportions == null) {
            sampleProportions = new Object2DoubleArrayMap<String>();
            sampleProportions.defaultReturnValue(-1);
        }
        final double proportion;

        if (!sampleProportions.containsKey(sample)) {
            int sumOverSamples = 0;
            for (final String s : samples()) {
                sumOverSamples += getSumOverlapCounts(s);
            }
            for (final String s : samples()) {
                final double sampleProportion =
                        ((double) getSumOverlapCounts(s)) / (double) sumOverSamples;
                this.sampleProportions.put(s, sampleProportion);
            }
        }
        proportion = sampleProportions.get(sample);
        assert proportion != -1 : " Proportion must be defined for sample " + sample;
        return proportion;
    }


}
