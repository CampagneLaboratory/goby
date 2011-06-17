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

package edu.cornell.med.icb.goby.counts;


import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.io.InputBitStream;

import java.io.ByteArrayInputStream;
import java.io.DataOutput;
import java.io.IOException;
import java.io.InputStream;

/**
 * Build an index for Counts information.
 * *
 * @author Fabien Campagne
 *         Date: 6/15/11
 *         Time: 5:56 PM
 */
public class CountIndexBuilder {
    IntArrayList positions = new IntArrayList();
    IntArrayList offsets = new IntArrayList();
    IntArrayList counts = new IntArrayList();
    private int numIndexEntries;
    private int transitionsPerIndex = 10000;

    public CountIndexBuilder(int transitionsPerIndex) {
        this.transitionsPerIndex = transitionsPerIndex;
    }

    public CountIndexBuilder() {

    }

    /**
     * Construct an index for counts data given as a byte array.
     *
     * @param countBytes The compressed counts data to read with CountsReader
     * @param indexPart  The data output where to write the index.
     * @throws IOException If an error occurs writing the index or reading the counts.
     */
    public void buildIndex(byte[] countBytes, DataOutput indexPart) throws IOException {

        final InputStream stream = new ByteArrayInputStream(countBytes);
        InputBitStream inputBitStream = new InputBitStream(stream);
        final CountsReaderI reader = new CountsReader(inputBitStream);
        int transitionNum = 0;
        long bitsWritten=0;
        while (reader.hasNextTransition()) {

            reader.nextTransition();

            int position = reader.getPosition();
            int count = reader.getCount();
            ++transitionNum;
            if (transitionNum % transitionsPerIndex == 0) {

                offsets.add((int) bitsWritten);
                positions.add(position);
                counts.add(count);
                // we cast the offset to an int because we do not expect to see offsets larger than a few million bits.

                ++numIndexEntries;
            }
            bitsWritten = inputBitStream.readBits();
        }
        indexPart.writeInt(numIndexEntries);
        BinIO.storeInts(positions.toIntArray(), indexPart);
        BinIO.storeInts(offsets.toIntArray(), indexPart);
        BinIO.storeInts(counts.toIntArray(), indexPart);

    }

}
