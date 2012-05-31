/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.methylation;

import it.unimi.dsi.lang.MutableString;

/**
 * @author Fabien Campagne
 *         Date: Oct 24, 2010
 *         Time: 12:02:47 PM
 */
public class MethylationSiteIterator {

    private MethylationData data;
    private int currentSiteIndex = -1;
    private char strand;
    private MethylationSite nextSite;
    private int chromosomeIndex;

    public void skipTo(String chromosomeSelected) {
        final MutableString chromosomeMutable = new MutableString(chromosomeSelected);
        skipTo(chromosomeMutable);

    }

    public void skipTo(MutableString chromosomeMutable) {
        chromosomeIndex = data.chromosomes.get(chromosomeMutable);
        if (chromosomeIndex == -1) {
            throw new IllegalArgumentException(String.format("Chromosome %d does not exist in dataset. ", chromosomeMutable));
        }
        if (nextSite != null && nextSite.chromosome == chromosomeIndex) return;
        else {
            nextSite = null;
            currentSiteIndex = 0;
            // advance current site index until the site it points to matches the selected chromosome:
            while (data.sites.get(currentSiteIndex).chromosome != chromosomeIndex) {
                currentSiteIndex++;
            }
        }
    }

    public MethylationSiteIterator(MethylationData data, char strand) {
        this.data = data;
        this.strand = strand;
    }


    public boolean hasNextSite() {
        if (nextSite != null) {
            return true;
        } else {
            do {
                currentSiteIndex++;
                if (currentSiteIndex >= data.sites.size()) {

                    // no more sites to iterate through
                    return false;
                }
                nextSite = data.sites.get(currentSiteIndex);
                if (nextSite.chromosome != chromosomeIndex) {
                    // no more sites on the selected iterator chromosome.
                    return false;
                }
            } while (nextSite.strand != strand);
            return true;
        }
    }

    public MethylationSite nextSite() {
        try {
            if (hasNextSite()) {
                return nextSite;

            } else {
                return null;
            }
        } finally {
            nextSite = null;

        }
    }

    public void skipToPosition(int positionStart) {
        if (nextSite != null && nextSite.position == positionStart) return;
        nextSite = null;
        while (hasNextSite()) {

            if (nextSite.position >= positionStart) return;
            nextSite = null;
        }
    }
}
