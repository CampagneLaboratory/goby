/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.util.Arrays;

/**
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:16:18 AM
 */
public class QualityScoreFilter extends GenotypeFilter {
    private byte scoreThreshold = 30;

    public String describe() {
        return "q<" + scoreThreshold;
    }

    int[] removed = new int[5];

    @Override
    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredList) {
        resetCounters();
        initStorage(sampleCounts.length);
        Arrays.fill(removed, 0);
        int removedVarCount = 0;
        for (PositionBaseInfo info : list) {
            numScreened++;
            if (!info.matchesReference && info.qualityScore < scoreThreshold) {
                if (!filteredList.contains(info)) {
                    filteredList.add(info);
                    final SampleCountInfo countInfo = sampleCounts[info.readerIndex];
                    final int baseIndex = countInfo.baseIndex(info.to);
                    countInfo.counts[baseIndex]--;
                    removed[baseIndex]++;
                    this.varCountRemovedPerSample[info.readerIndex]++;
                    numFiltered++;
                }
            }
        }
        // adjust refCount and varCount:
        adjustRefVarCounts(sampleCounts);
    }


}
