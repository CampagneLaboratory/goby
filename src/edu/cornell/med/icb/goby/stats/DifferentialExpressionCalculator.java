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
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2DoubleArrayMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 3:35:24 PM
 */
public class DifferentialExpressionCalculator {
    private final ObjectSet<String> groups;
    private final Object2ObjectMap<String, String> sampleToGroupMap;
    private final IndexedIdentifier elementLabels;
    private int elementsPerSample;
    private int numberOfSamples;
    private final Object2ObjectMap<String, IntArrayList> sampleToCounts =
            new Object2ObjectOpenHashMap<String, IntArrayList>();
    private final Object2ObjectMap<String, DoubleArrayList> sampleToRPKMs =
            new Object2ObjectOpenHashMap<String, DoubleArrayList>();
    private Object2DoubleMap<String> sampleProportions;

    public DifferentialExpressionCalculator() {
        super();
        groups = new ObjectArraySet<String>();
        elementLabels = new IndexedIdentifier();
        sampleToGroupMap = new Object2ObjectOpenHashMap<String, String>();
    }

    public void defineGroup(final String label) {
        groups.add(label);
    }

    public int defineElement(final String label) {
        return elementLabels.registerIdentifier(new MutableString(label));
    }

    public void associateSampleToGroup(final String sample, final String group) {
        sampleToGroupMap.put(sample, group);
    }

    public void observe(final String sample, final String elementId, final int count, final double RPKM) {
        IntArrayList counts = sampleToCounts.get(sample);
        DoubleArrayList rpkms = sampleToRPKMs.get(sample);

        if (counts == null) {
            counts = new IntArrayList(elementsPerSample);
            counts.size(elementsPerSample);
            sampleToCounts.put(sample, counts);
        }
        if (rpkms == null) {
            rpkms = new DoubleArrayList(elementsPerSample);
            rpkms.size(elementsPerSample);
            sampleToRPKMs.put(sample, rpkms);
        }

        final int elementIndex = elementLabels.get(new MutableString(elementId));
        counts.set(elementIndex, count);
        rpkms.add(elementIndex, RPKM);
    }

    public DifferentialExpressionResults compare(final DifferentialExpressionResults results,
                                                 final StatisticCalculator tester,
                                                 final String... group) {
        assert !tester.canDo(group) : "The number of groups to compare is not supported by the specified calculator.";
        tester.setResults(results);
        return tester.evaluate(this, results, group);
    }

    public DifferentialExpressionResults compare(final StatisticCalculator tester, final String... group) {
        final DifferentialExpressionResults results = new DifferentialExpressionResults();

        assert !tester.canDo(group) : "The number of groups to compare is not supported by the specified calculator.";
        tester.setResults(results);
        return tester.evaluate(this, results, group);
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
     * Get the stored RPKM for an element in a given sample.
     *
     * @param sample
     * @param elementId
     * @return
     */
    public double getRPKM(final String sample, final MutableString elementId) {
        final DoubleArrayList rpkms = sampleToRPKMs.get(sample);
        if (rpkms == null) {
            return 0;
        }
        final int elementIndex = elementLabels.get(new MutableString(elementId));

        return rpkms.get(elementIndex);
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
        final int elementIndex = elementLabels.get(new MutableString(elementId));

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
        final ObjectSet<String> sampleSet = sampleToGroupMap.keySet();
        return sampleSet.toArray(new String[sampleSet.size()]);
    }

    /**
     * Returns the proportion of counts that originate from a certain sample. The number of counts in the sample
     * is divided by the sum of counts over all the samples in the experiment.
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
