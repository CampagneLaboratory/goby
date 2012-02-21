/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

import edu.cornell.med.icb.goby.stats.SamplePair;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

/**
 * @author Fabien Campagne
 *         Date: 2/19/12
 *         Time: 1:24 PM
 */
public class IntraGroupEnumerator {
    private int[] sampleIndexToGroupIndex;
    private int numSamples;
    private int numGroups;

    public IntraGroupEnumerator(final int[] sampleIndexToGroupIndex, int numSamples, int numGroups) {
        samplePairsForGroup = new ObjectArrayList[numGroups];
        this.sampleIndexToGroupIndex = sampleIndexToGroupIndex;
        this.numSamples = numSamples;
        this.numGroups = numGroups;

    }

    public ObjectArrayList<SamplePair> getPairs(int groupIndex) {
        final ObjectArrayList<SamplePair> samplePairs = samplePairsForGroup[groupIndex];
        assert samplePairs != null : "sample pairs must have been defined for group " + groupIndex;
        return samplePairs;
    }

    public void recordPairForGroup(int groupIndex) {
        samplePairsForGroup[groupIndex] = new ObjectArrayList<SamplePair>();
        IntSet sampleIndicesInGroup = new IntAVLTreeSet();
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            if (sampleIndexToGroupIndex[sampleIndex] == groupIndex) {
                sampleIndicesInGroup.add(sampleIndex);
            }
        }
        for (int sampleIndexA : sampleIndicesInGroup) {
            for (int sampleIndexB : sampleIndicesInGroup) {
                if (sampleIndexA < sampleIndexB) {
                    samplePairsForGroup[groupIndex].add(new SamplePair(sampleIndexA, sampleIndexB));
                }
            }
        }
    }

    private ObjectArrayList<SamplePair>[] samplePairsForGroup;
}
