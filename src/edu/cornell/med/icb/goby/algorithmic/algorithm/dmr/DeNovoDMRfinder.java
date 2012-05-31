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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.List;


/**
 * Find differentially methylated regions(DMRs) across groups of samples without prior
 * annotation information
 *
 * @Author: Nyasha Chambwe
 * @Date: 3/22/12
 * @Time: 3:18 PM
 */
public class DeNovoDMRfinder {

    /**
     * Threshold to consider a differentially methylated region significant
     */
    private double pThreshold;

    private ObjectArrayList<WindowRange> dmrResultList = new ObjectArrayList<WindowRange>();
    private WindowRange searchSpace;
    private int windowLength;
    private WindowRange fullRange;
    private PVAlueProvider provider;

    public DeNovoDMRfinder(final double pThreshold, int windowLength,
                           final PVAlueProvider provider) {
        this.pThreshold = pThreshold;
        this.windowLength = windowLength;
        this.provider = provider;
        fullRange = new WindowRange(0, windowLength);
    }


    public ObjectArrayList<WindowRange> search(final GroupComparison groupComp) {
        // assert that methylatedCounts length >= 1
        searchSpace = fullRange;
        return search(groupComp, searchSpace, true, true, false);
    }

    public List<WindowRange> getDMRs() {
        return dmrResultList;
    }

    private ObjectArrayList<WindowRange> search(final GroupComparison groupComp,
                                                final WindowRange currentWindow,
                                                final boolean searchFirst,
                                                final boolean searchMiddle,
                                                final boolean searchLast) {

        if (currentWindow.length >= 1) {

            if (searchFirst) {
                searchFirst(groupComp, currentWindow);
            }
            if (searchMiddle) {
                searchMiddle(groupComp, currentWindow);
            }
            if (searchLast) {
                searchLast(groupComp, currentWindow);
            }
        }
        return dmrResultList;
    }

    private void searchFirst(final GroupComparison groupComp, final WindowRange intervalToSearch) {
        final int queryWindowSize = Math.round(intervalToSearch.length / 3);
        final WindowRange first = new WindowRange(intervalToSearch.start, intervalToSearch.start + queryWindowSize);
        final double pFirst = provider.getPValue(first.start, first.end, groupComp);
        first.setPforRange(pFirst);
        search(groupComp, first);
    }

    private void searchMiddle(final GroupComparison groupComp, final WindowRange intervalToSearch) {
        final int queryWindowSize = Math.round(intervalToSearch.length / 3);
        final int newStart = intervalToSearch.start + queryWindowSize;
        final int newEnd = newStart + queryWindowSize;
        final WindowRange middle = new WindowRange(newStart, newEnd);
        final double pMiddle = provider.getPValue(middle.start, middle.end, groupComp);
        middle.setPforRange(pMiddle);
        search(groupComp, middle);
    }

    private void searchLast(final GroupComparison groupComp, WindowRange intervalToSearch) {
        final int queryWindowSize = Math.round(intervalToSearch.length / 3);
        final int newStart = intervalToSearch.start + 2 * intervalToSearch.start;
        final int newEnd = newStart + queryWindowSize;
        final WindowRange last = new WindowRange(newStart, newEnd);
        final double pEnd = provider.getPValue(last.start, last.end, groupComp);
        last.setPforRange(pEnd);
        search(groupComp, intervalToSearch);
    }

    private void search(final GroupComparison groupComp, final WindowRange intervalToSearch) {
        if (intervalToSearch.pForRange < pThreshold) {
            dmrResultList.add(intervalToSearch);
        } else {
            search(groupComp, intervalToSearch, false, true, false);
        }
    }


}
