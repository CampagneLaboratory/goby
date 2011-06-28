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
import edu.rit.mp.buf.ObjectArrayBuf;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import javax.swing.*;

/**
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 11:37:42 AM
 */
public class SampleCountInfo {
    public static final int BASE_A_INDEX = 0;
    public static final int BASE_T_INDEX = 1;
    public static final int BASE_C_INDEX = 2;
    public static final int BASE_G_INDEX = 3;
    public static final int BASE_OTHER_INDEX = 4;
    public static final int BASE_MAX_INDEX = BASE_OTHER_INDEX + 1;

    public char referenceBase;
    public IntSet distinctReadIndices = new IntArraySet();
    public int sampleIndex;
    public int varCount;
    public int refCount;
    /**
     * Number of bases that failed base filters, for any reason.
     */
    public int failedCount;
    public int[] counts = new int[5];
    /**
     * List of indel eirs that start at the position in this sample.
     */
    private ObjectArrayList<EquivalentIndelRegion> indels;

    /**
     * Add an indel to this sample info. This method lazily creates the indels collection when the first
     * indel is added. If the same indel was already added (equality by range of eir start-end), the frequency
     * is incremented.
     *
     * @param indel eir observed starting at the position this sampleCountInfo is associated with.
     */
    public void addIndel(final EquivalentIndelRegion indel) {

        if (indels == null) {
            indels = new ObjectArrayList<EquivalentIndelRegion>();
        }
        for (final EquivalentIndelRegion prevIndel : indels) {
            if (prevIndel.equals(indel)) {
                prevIndel.frequency+=1;
            }
        }
        indels.add(indel);
    }

    /**
     * Return true if this sample has indel observations.
     * @return  True or false.
     */
    public boolean hasIndels() {
        return !(indels==null);
    }

    public final char base(final int baseIndex) {
        switch (baseIndex) {
            case BASE_A_INDEX:
                return 'A';
            case BASE_C_INDEX:
                return 'C';
            case BASE_T_INDEX:
                return 'T';
            case BASE_G_INDEX:
                return 'G';
            default:
                return 'N';
        }
    }

    public final int baseIndex(final char to) {
        switch (to) {
            case 'A':
                return BASE_A_INDEX;
            case 'C':
                return BASE_C_INDEX;
            case 'T':
                return BASE_T_INDEX;
            case 'G':
                return BASE_G_INDEX;
            default:
                return BASE_OTHER_INDEX;
        }
    }

    public ObjectArrayList<EquivalentIndelRegion> getEquivalentIndelRegions() {
        return indels;
    }
}
