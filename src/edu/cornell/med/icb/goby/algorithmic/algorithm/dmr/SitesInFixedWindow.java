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

import it.unimi.dsi.fastutil.ints.IntArrayList;

/**
 * Class to count how many sites are observed in a moving window of fixed length.
 *
 * @author Fabien Campagne
 *         Date: 3/1/12
 *         Time: 11:16 AM
 */
public class SitesInFixedWindow {

    private final int D;
    private IntArrayList refIndices = new IntArrayList();
    private IntArrayList positions = new IntArrayList();

    /**
     * Construct an instance with fixed window of length D.
     *
     * @param D length of the fixed window.
     */
    public SitesInFixedWindow(int D) {
        this.D = D;
    }

    /**
     * Add a site to the window.
     *
     * @param referenceIndex index of the reference sequence where the site is located.
     * @param position       position of the site on the reference sequence.
     */
    public void add(final int referenceIndex, final int position) {
        prune(referenceIndex, position);
        refIndices.add(referenceIndex);
        positions.add(position);
    }

    /**
     * Move the window to a new forward location and prune sites that are outside the window.
     *
     * @param referenceIndex index of the reference sequence where the window should be moved.
     * @param position       position of the site where the window should be moved.
     */
    public void prune(final int referenceIndex, final int position) {
       int size=n();
        if (size == 0) {
            return;
        }
        int pruneIndex = 0;
        int i=0;

        while (i<size && outsideOfWindow(refIndices.getInt(i), positions.getInt(i), referenceIndex, position) ) {
            pruneIndex++;
            i++ ;
        }
        if (pruneIndex>0) {
        refIndices.removeElements(0, pruneIndex);
        positions.removeElements(0, pruneIndex);
        }
    }

    /**
     * @param refIndexA
     * @param positionA
     * @param referenceIndex
     * @param position       position must be larger than positionA when both are on the same reference sequence.
     * @return
     */
    private boolean outsideOfWindow(int refIndexA, int positionA, int referenceIndex, int position) {
        if (refIndexA != referenceIndex) {
            return true;
        }
        assert position > positionA : "query position must be larger than head position.";
        return position - positionA > D;

    }

    /**
     * Return the number of sites in the window.
     *
     * @return n the number of sites currently in the window.
     */
    public int n() {
        int size = refIndices.size();
        assert size == positions.size() : "refIndices and positions must have the same size!";
        return size;
    }

}
