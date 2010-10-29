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

package edu.cornell.med.icb.goby.reads;

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

/**
 * @author Fabien Campagne
 *         Date: Jun 20, 2009
 *         Time: 10:51:11 AM
 */
public class MultiReads {
    private int readIndex = -1;
    private int currentByteIndex;
    private final ByteArrayList bytes = new ByteArrayList();
    private final IntArrayList ends = new IntArrayList();
    private final IntArrayList starts = new IntArrayList();

    public void newRead() {
        // previous read ends at previous index.
        ends.set(readIndex, currentByteIndex - 1);
        readIndex += 1;
        starts.set(readIndex, currentByteIndex);
    }

    public void addByte(final byte base) {
        currentByteIndex += 1;
        bytes.add(base);
    }
}
