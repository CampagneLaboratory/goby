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

import java.util.Arrays;

/**
 * Holds information for counts of methylated and un-methylated bases at a position. Info is organized per sample and per group.
 *
 * @author Fabien Campagne
 *         Date: 1/28/12
 *         Time: 2:45 PM
 */
public class MethylCountInfo {
    /**
     * One count per sample for the number of bases that are not methylated at this position.
     */
    public int[] unmethylatedCCountPerSample;
    /**
     * One count per sample for the number of bases that are methylated at this position.
     */
    public int[] methylatedCCountPerSample;
    /**
     * One count per group for the number of bases that are not methylated at this position.
     */
    public int[] unmethylatedCCountPerGroup;
    /**
     * One count per group for the number of bases that are methylated at this position.
     */
    public int[] methylatedCCountPerGroup;
    /**
     * Number of conversion/non-conversion events observed at this site.
     */
    public int eventCountAtSite = 0;
    /**
     * Strand '+' or '-' at the given site.
     */
    public char strandAtSite;

    /**
     * Construct a new MethylCountInfo instance.
     * @param numberOfSamples The number of samples analyzed.
     * @param numberOfGroups  The number of groups analyzed.
     */
    public MethylCountInfo(final int numberOfSamples, final int numberOfGroups) {
        unmethylatedCCountPerGroup = new int[numberOfGroups];
        methylatedCCountPerGroup = new int[numberOfGroups];
        methylatedCCountPerSample = new int[numberOfSamples];
        unmethylatedCCountPerSample = new int[numberOfSamples];
    }

    /**
     * Resets all counts to zero and strand to unknown.
     */
    public void reset() {
        Arrays.fill(methylatedCCountPerGroup, 0);
        Arrays.fill(unmethylatedCCountPerGroup, 0);
        Arrays.fill(methylatedCCountPerSample, 0);
        Arrays.fill(unmethylatedCCountPerSample, 0);
        eventCountAtSite = 0;

        strandAtSite = '?';
    }
}
