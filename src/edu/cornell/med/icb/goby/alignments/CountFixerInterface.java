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
 * Implementations of this interface provide strategies to eliminate or correct likely errors.
 *
 * @author Fabien Campagne
 *         Date: Mar 23, 2011
 *         Time: 11:09:25 AM
 */
public interface CountFixerInterface {
    /**
     * Implementations of this method decide how to fix list to eliminate or correct likely errors.
     *
     * @param list         List of variations or reference bases.
     * @param sampleCounts allele frequencies per sample
     * @param likelyErrors List of suspicious variations or reference bases.
     */
    void fix(ObjectArrayList<PositionBaseInfo> list, SampleCountInfo[] sampleCounts,
             ObjectList<PositionBaseInfo> likelyErrors);
}
