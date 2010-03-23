/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

/**
 * A set of read/query indices and associated multiplicity. Multiplicity is an int associated
 * with each read. It can represent anything, but is often used to stored the number of times
 * a read is repeated in a sample file or dataset.  The read indices are stored with delta
 * coding, and multiplicity values with gamma coding, keeping the file representation
 * of the set small even when there are million of reads.
 *
 * @author Fabien Campagne
 *         Date: Jun 3, 2009
 *         Time: 3:22:49 PM
 */
public class ReadSet {
    private final IntSet filter = new IntOpenHashSet();
    private final Int2IntMap multiplicityMap;
    private int smallestStoredMultiplicity;

    public ReadSet() {
        super();
        multiplicityMap = new Int2IntOpenHashMap();
    }

    public void clear() {
        filter.clear();
        multiplicityMap.clear();
    }

    public void add(final int readIndex, final int multiplicity) {
        assert readIndex >= 0 : "read indices must be positive";
        filter.add(readIndex);
     //   System.out.println("adding " + readIndex + " " + multiplicity);
        if (multiplicity >= smallestStoredMultiplicity) {

            multiplicityMap.put(readIndex, multiplicity);
        }
    }

    public boolean contains(final int readIndex) {
        return filter.contains(readIndex);
    }

    /**
     * Save the filter to disk.  The name of the actual file writen will be of the form
     * "${basename}-${suffix}.filter"
     *
     * @param basename basename of the filter file to load
     * @param suffix   suffix of the filter file to load
     * @throws IOException if the file cannot be found or parsed.
     */
    public void save(final String basename, final String suffix) throws IOException {
        final String filename = basename + "-" + suffix + ".filter";
        final OutputBitStream out = new OutputBitStream(filename);

        // sort in increasing order:
        final IntList sorted = new IntArrayList();
        sorted.addAll(filter);
        Collections.sort(sorted);

        // append smallest stored multiplicity
        out.writeGamma(smallestStoredMultiplicity);

        // append to output:
        out.writeGamma(sorted.size());
        int previous = -1;
        for (final int readIndex : sorted) {
            out.writeDelta(readIndex - previous);
            final int multiplicity = multiplicityMap.get(readIndex);
      //      System.out.println(readIndex + " " + multiplicity);
            out.writeGamma(multiplicity);
            previous = readIndex;

        }
        out.close();
    }

    /**
     * Load the filter from disk.  The name of the actual file loaded will be of the form
     * "${basename}-${suffix}.filter"
     *
     * @param basename basename of the filter file to load
     * @param suffix   suffix of the filter file to load
     * @throws IOException if the file cannot be found or parsed.
     */
    public void load(final String basename, final String suffix) throws IOException {
        final String filename = basename + "-" + suffix + ".filter";
        load(new File(filename));
    }

    /**
     * Load a read set from disk.
     *
     * @param file File to load.
     * @throws IOException If an error occurs reading the read set file.
     */
    public void load(final File file) throws IOException {
        InputBitStream in = null;
        try {
            in = new InputBitStream(file);
            filter.clear();
            smallestStoredMultiplicity = in.readGamma();

            final int numDeltas = in.readGamma();
            int previous = -1;
            for (int i = 0; i < numDeltas; i++) {
                final int delta = in.readDelta();
                final int multiplicity = in.readGamma();
                final int readIndex = previous + delta;
                filter.add(readIndex);
                multiplicityMap.put(readIndex, multiplicity);
                previous = readIndex;
            }
        } finally {
            if (in != null) {
                in.close();
            }
        }
    }

    public void add(final int readIndex) {
        add(readIndex, 0);
    }

    public int size() {
        return filter.size();
    }

    /**
     * Returns the multiplicity of the read, or zero is the read is not in this read set.
     *
     * @param readIndex
     * @return
     */
    public int getMultiplicity(final int readIndex) {
        if (filter.contains(readIndex)) {
            return multiplicityMap.get(readIndex);
        } else {
            return 0;
        }
    }

    /**
     * Add all the read indices with the given multiplicity.
     *
     * @param otherReadIndices Set of read indices to add.
     * @param multiplicity
     */
    public void addAll(final IntSet otherReadIndices, final int multiplicity) {
        for (final int readIndex : otherReadIndices) {
            add(readIndex, multiplicity);
        }
    }

    /**
     * Returns the maximum query index stored in this filter.
     * @return
     */
    public int getMaxReadIndex() {
        int maxQueryIndex = -1;
        final IntIterator intIterator = filter.iterator();
        while (intIterator.hasNext()) {
            int queryIndex = intIterator.next();
            maxQueryIndex = Math.max(maxQueryIndex, queryIndex);

        }
        return maxQueryIndex;
    }

    /**
     * Set the default multiplicity value.
     *
     * @param value
     */
    public void smallestStoredMultiplicity(final int value) {
        multiplicityMap.defaultReturnValue(value);
        smallestStoredMultiplicity = value;
    }


}
