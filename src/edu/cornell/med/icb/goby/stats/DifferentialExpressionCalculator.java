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

import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.lang.MutableString;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

/**
 * @author Fabien Campagne
 *         Date: Jan 11, 2010
 *         Time: 3:35:24 PM
 */
public class DifferentialExpressionCalculator {

    ObjectSet<String> groups;
    Object2ObjectMap<String, String> sampleToGroupMap;
    IndexedIdentifier elementLabels;
    private int elementsPerSample;
    private int numberOfSamples;
    private Object2ObjectMap<String, IntArrayList> sampleToCounts = new Object2ObjectOpenHashMap<String, IntArrayList>();
    private Object2ObjectMap<String, DoubleArrayList> sampleToRPKMs = new Object2ObjectOpenHashMap<String, DoubleArrayList>();


    @Override
    public String toString() {
        return super.toString();
    }

    public DifferentialExpressionCalculator() {
        groups = new ObjectArraySet<String>();
        elementLabels = new IndexedIdentifier();
        sampleToGroupMap = new Object2ObjectOpenHashMap<String, String>();
    }

    public void defineGroup(String label) {
        groups.add(label);
    }

    public int defineElement(String label) {
        return elementLabels.registerIdentifier(new MutableString(label));
    }

    public void associateSampleToGroup(String sample, String group) {
        sampleToGroupMap.put(sample, group);
    }

    public void observe(String sample, String elementId, int count, double RPKM) {
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

    public DifferentialExpressionResults compare(DifferentialExpressionResults results, StatisticCalculator tester, String... group) {


        assert !tester.canDo(group) : "The number of groups to compare is not supported by the specified calculator.";
        tester.setResults(results);
        return tester.evaluate(this, results, group);
    }

    public DifferentialExpressionResults compare(StatisticCalculator tester, String... group) {

        DifferentialExpressionResults results = new DifferentialExpressionResults();

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
    public void reserve(int elementsPerSample, int numberOfSamples) {
        this.elementsPerSample = elementsPerSample;
        this.numberOfSamples = numberOfSamples;
    }

    public ObjectArraySet<String> getSamples(String groupA) {
        ObjectArraySet<String> samples = new ObjectArraySet<String>();
        for (String sampleId : sampleToGroupMap.keySet()) {
            if (sampleToGroupMap.get(sampleId).equals(groupA)) samples.add(sampleId);
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
    public double getRPKM(String sample, MutableString elementId) {
        DoubleArrayList rpkms = sampleToRPKMs.get(sample);
        if (rpkms == null) return 0;
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
    public int getOverlapCount(String sample, MutableString elementId) {
        IntArrayList counts = sampleToCounts.get(sample);
        if (counts == null) return 0;
        final int elementIndex = elementLabels.get(new MutableString(elementId));

        return counts.get(elementIndex);
    }

    /**
     * Returns the sum of counts in a given sample.
     */
    public int getSumOverlapCounts(String sample) {
        int sumCounts = 0;
        IntArrayList counts = sampleToCounts.get(sample);
        if (counts == null) return 0;
        for (int count : counts) {
            sumCounts += count;
        }
        return sumCounts;
    }

    public String[] samples() {
        final ObjectSet<String> sampleSet = sampleToGroupMap.keySet();
        return sampleSet.toArray(new String[sampleSet.size()]);
    }
}
