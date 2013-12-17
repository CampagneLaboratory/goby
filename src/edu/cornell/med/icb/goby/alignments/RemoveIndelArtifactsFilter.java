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

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import it.unimi.dsi.fastutil.objects.ObjectSet;

/**
 * Remove base observations that are a direct consequence of indels spanning the location.
 *
 * @author campagne
 *         Date: 7/1/11
 *         Time: 4:00 PM
 */
public class RemoveIndelArtifactsFilter extends GenotypeFilter {
    @Override
    public void filterGenotypes(DiscoverVariantPositionData list, SampleCountInfo[] sampleCounts,
                                ObjectSet<PositionBaseInfo> filteredList) {
        resetCounters();
        initStorage(sampleCounts.length);

        if (list.hasCandidateIndels()) {
            for (EquivalentIndelRegion indel : list.getIndels()) {
                for (final PositionBaseInfo info : list) {


                    if (info.readerIndex == indel.sampleIndex && info.to == '-' ) {
                        //     list.remove(info);

                        final SampleCountInfo sampleCountInfo = sampleCounts[info.readerIndex];
                        final char base = info.to == '-'? info.from : info.to;
                        final int baseIndex = sampleCountInfo.baseIndex(base);
                        if (!filteredList.contains(info)) {
                            sampleCountInfo.suggestRemovingGenotype(baseIndex);
                            removeGenotype(info, filteredList);
                        }
                        /*  if (varCountRemovedPerSample[info.readerIndex] > sampleCount.varCount) {

                            System.out.println("STOP3");
                        }*/
                    }
                }
            }

        }
        numScreened += list.size();
     //   filterIndels(list, sampleCounts);
     // do not adjust for thresholding effects.
        adjustGenotypes(list, filteredList, sampleCounts);
        adjustRefVarCounts(sampleCounts);
    }

    @Override
    public int getThresholdForSample(int sampleIndex) {
        throw new UnsupportedOperationException("This filter does not support method getThresholdForSample()");
    }

    private void checkCountPositive(SampleCountInfo sampleCount, int baseIndex) {
       /* if (sampleCount.counts[baseIndex] < 0) {
            System.out.println("STOP");
        } */
    }
}
