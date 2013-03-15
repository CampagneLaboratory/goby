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

import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.util.Arrays;

/**
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:16:18 AM
 */
public class QualityScoreFilter extends GenotypeFilter {
    private byte scoreThreshold = 30;
    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(QualityScoreFilter.class, "scoreThreshold:Phred score threshold to keep bases.:30");


    public static DynamicOptionClient doc() {
        return doc;
    }

    public QualityScoreFilter() {
        scoreThreshold = doc.getByte("scoreThreshold");
    }

    public String describe() {
        return "q<" + scoreThreshold;
    }

    int thresholdPerSample[];

    @Override
    public int getThresholdForSample(int sampleIndex) {
        return thresholdPerSample[sampleIndex];
    }

    int[] removed = new int[5];


    @Override
    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredList) {
        resetCounters();
        initStorage(sampleCounts.length);
        if (thresholdPerSample == null) {
            thresholdPerSample = new int[sampleCounts.length];

        } else {
            Arrays.fill(thresholdPerSample, 0);
        }
        for (final PositionBaseInfo info : list) {
            numScreened++;
            if (!info.matchesReference && info.qualityScore < scoreThreshold) {
                if (!filteredList.contains(info)) {
                    if (info.to != '-' && info.from != '-') {
                        // indels have a quality score  of zero but should not be removed at this stage.
                        final SampleCountInfo sampleCountInfo = sampleCounts[info.readerIndex];
                        final int baseIndex = sampleCountInfo.baseIndex(info.to);

                        sampleCountInfo.suggestRemovingGenotype(baseIndex);
                        removeGenotype(info, filteredList);
                        thresholdPerSample[info.readerIndex]++;
                    }
                }
            }
        }
        filterIndels(list, sampleCounts);
        /*
       TODO: enable this when we store quality score for context of indels:
       if (list.hasCandidateIndels()) {
            // remove candidate indels if they don't make the base quality threshold (threshold determined by bases observed
            // at that position):
            for (final EquivalentIndelRegion indel : list.getIndels()) {
                for (final byte baseQuality : indel.getQualityScoresInContext()) {
                    if (baseQuality < scoreThreshold) {
                        list.failIndel(indel);
                        if (indel.matchesReference()) {
                            refCountRemovedPerSample[indel.sampleIndex]++;
                        } else {
                            varCountRemovedPerSample[indel.sampleIndex]++;
                        }
                        break;
                    }
                }
            }
        }

        */
        adjustGenotypes(list, filteredList, sampleCounts);

        // adjust refCount and varCount:
        adjustRefVarCounts(sampleCounts);
    }


}
