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

import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;

/**
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 11:37:42 AM
 */
public class SampleCountInfo {
    static final public int BASE_A_INDEX = 0;
    static final public int BASE_T_INDEX = 1;
    static final public int BASE_C_INDEX = 2;
    static final public int BASE_G_INDEX = 3;
    static final public int BASE_OTHER_INDEX = 4;
    static final public int BASE_MAX_INDEX = BASE_OTHER_INDEX;

    public char referenceBase;
    public IntSet distinctReadIndices = new IntArraySet();
    public int sampleIndex;
    public int varCount;
    public int refCount;
    public int[] counts = new int[5];

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
}
