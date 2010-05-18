/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntList;

import java.util.Map;
import java.io.IOException;

import edu.cornell.med.icb.goby.counts.CountsWriter;

/**
 * @author Fabien Campagne
 *         Date: May 16, 2010
 *         Time: 12:55:53 PM
 */
public interface ComputeCountInterface {
    void startPopulating();

    void accumulate();

    void baseCount(final CountsWriter writer) throws IOException;

    void baseCount();

    Map getCountPerBase();

    IntList getCountKeys();

    void populate(int startPosition, int endPosition, boolean forwardStrand, int queryIndex);
}
