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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

/**
 * Provide a strategy for filtering bases and reduce the impact of sequencing errors on downstream
 * statistics.
 *
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:07:11 AM
 */
public abstract class BaseFilter {

    /**
     * Adjust list and sampleCounts to remove/reduce the effect of sequencing errors.
     * @param list Variation or reference bases at position
     * @param sampleCounts Counts for alleles at position each each sample under study.
     */
    public abstract void filterBases(ObjectArrayList<PositionBaseInfo> list,
                      SampleCountInfo[] sampleCounts,
                      ObjectArrayList<PositionBaseInfo> filteredList);

    /**
     * Returns a short description of the fitlering criteria.
      * @return a short description of the fitlering criteria.
     */
    public String describe() {
        return this.getClass().getSimpleName();
    }
    int numFiltered=0;
    int numScreened=0;

    public double getPercentFilteredOut() {
        double rate=numFiltered;
        rate/=numScreened;
        return rate*100d;
    }

    public String getName() {
        return "filter ("+describe()+")";
    }

     void resetCounters() {
        numScreened=0;
        numFiltered=0;
    }

    protected void adjustRefVarCounts(SampleCountInfo[] sampleCounts, int[] removed) {
       for (SampleCountInfo sci : sampleCounts) {
            final int refBaseIndex = sci.baseIndex(sci.referenceBase);

            for (int otherBaseIndex = 0; otherBaseIndex < removed.length; otherBaseIndex++) {
                if (refBaseIndex == otherBaseIndex) {
                    sci.refCount -= removed[refBaseIndex];
                } else {
                    sci.varCount -= removed[otherBaseIndex];
                }
            }
            assert sci.refCount >= 0;
            assert sci.varCount >= 0;
        }
    }
}
