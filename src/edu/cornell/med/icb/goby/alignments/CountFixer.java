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

import java.util.Arrays;

/**
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:08:56 AM
 */
public class CountFixer implements CountFixerInterface {
    /**
     * The default implementation removes the likely errors from list. Override this method to implement alternative
     * treatment of errors (e.g., patching of the base to one of the types already detected without errors).
     *
     * @param list
     * @param sampleCounts
     * @param likelyErrors
     */
    @Override
    public void fix(final DiscoverVariantPositionData list,
                    final SampleCountInfo[] sampleCounts,
                    final ObjectSet<PositionBaseInfo> likelyErrors) {

        // do not decrement counts again. The Filters have done this already..

        for (SampleCountInfo sci : sampleCounts) {
            for (int i = 0; i < sci.counts.length; i++) {
                assert sci.counts[i] >= 0 : ("Counts must never be negative. This would happen if a GenotypeFilter removed counts directly. Value was " + sci.counts[i]);
            }
        }
        // calculate failed Count in each sample:
        for (PositionBaseInfo failed : likelyErrors) {
            ++(sampleCounts[failed.readerIndex].failedCount);
        }
        for (final EquivalentIndelRegion failedIndel : list.getFailedIndels()) {
            sampleCounts[failedIndel.sampleIndex].removeIndel(failedIndel);
            ++(sampleCounts[failedIndel.sampleIndex].failedCount);
        }
        list.removeAll(likelyErrors);
        for (SampleCountInfo sci : sampleCounts) {
            sci.clearFiltered();
        }
    }

}
