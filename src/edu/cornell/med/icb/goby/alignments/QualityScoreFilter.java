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

import edu.cornell.med.icb.goby.util.DynamicOptionClient;
import it.unimi.dsi.fastutil.objects.ObjectSet;

/**
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:16:18 AM
 */
public class QualityScoreFilter extends GenotypeFilter  {
    private byte scoreThreshold = 30;
    public static DynamicOptionClient doc=new DynamicOptionClient("scoreThreshold:Phred score threshold to keep bases.:30");

    public QualityScoreFilter() {
       scoreThreshold=doc.getByte("scoreThreshold");
    }

    public String describe() {
        return "q<" + scoreThreshold;
    }

    @Override
    public int getThresholdForSample(int sampleIndex) {
        throw new UnsupportedOperationException("This filter does not support method getThresholdForSample()");
    }

    int[] removed = new int[5];



    @Override
    public void filterGenotypes(DiscoverVariantPositionData list,
                                SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredList) {
        resetCounters();
        initStorage(sampleCounts.length);
        /*ByteArrayList bytes = new ByteArrayList();
        int position = 0;
        for (final PositionBaseInfo info : list) {
            position = info.position;
            bytes.add(info.qualityScore);
        }
        System.out.println(position + ": " + bytes.toString());*/

        for (final PositionBaseInfo info : list) {
            numScreened++;
            if (!info.matchesReference && info.qualityScore < scoreThreshold) {
                if (!filteredList.contains(info)) {
                    if (info.to != '-') {
                        // deleted bases have a quality score  of zero but should not be removed at this stage.
                        filteredList.add(info);
                        final SampleCountInfo countInfo = sampleCounts[info.readerIndex];
                        final int baseIndex = countInfo.baseIndex(info.to);
                        countInfo.counts[baseIndex]--;
                        varCountRemovedPerSample[info.readerIndex]++;
                        numFiltered++;
                    }
                }
            }
        }
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
        // adjust refCount and varCount:
        adjustRefVarCounts(sampleCounts);
    }


}
