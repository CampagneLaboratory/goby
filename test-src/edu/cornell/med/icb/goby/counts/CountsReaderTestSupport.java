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

import java.io.IOException;
import java.util.NoSuchElementException;

/**
 * An implementation of CountsReader to facilitate writting JUnit tests.
 *
 * @author Fabien Campagne
 *         Date: May 21, 2011
 *         Time: 11:09:20 AM
 */
public class CountsReaderTestSupport implements CountsReaderI {
    private int[] lengths;
    private int[] counts;
    private int index;
    private int size;
    private int currentPosition = 0;

    /**
     * Creates a CountsReader with arrays of counts and lengths. Paired elements of the array
     * provide count transaction information.
     *
     * @param format
     */
    public CountsReaderTestSupport(String format) {
        String[] tokens = format.split("[() ]+");
        IntArrayList lengths = new IntArrayList();
        IntArrayList counts = new IntArrayList();
        for (int i = 0; i < tokens.length; i++) {
            String token = tokens[i];
            if (token.length() > 0) {
                String[] t = token.split(",");
                lengths.add(Integer.parseInt(t[0]));
                counts.add(Integer.parseInt(t[1]));
            }
        }
        // System.out.printf("lengths: %s counts: %s", lengths, counts);
        init(lengths.toIntArray(), counts.toIntArray());
    }

    /**
     * Creates a CountsReader with arrays of counts and lengths. Paired elements of the array
     * provide count transaction information.
     *
     * @param lengths
     * @param counts
     */
    public CountsReaderTestSupport(int[] lengths, int[] counts) {
        init(lengths, counts);
    }

    private void init(int[] lengths, int[] counts) {
        this.lengths = lengths;
        this.counts = counts;
        this.index = -1;
        this.size = lengths.length;
        if (lengths.length != counts.length)
            throw new IllegalArgumentException("lengths and size must have same dimension.");
    }

    public int getPosition() {
        //      System.out.println("returning position="+currentPosition);
        return currentPosition;
    }

    public boolean hasNextTransition() throws IOException {

        return index + 1 < size;
    }

    public void nextTransition() throws IOException {
        if (!hasNextTransition()) {
            throw new NoSuchElementException("No such element.");
        }
        currentPosition +=  index>=0?lengths[index]:0;

        index++;


    }

    public int getCount() {
        //    System.out.println("returning count="+counts[index]);
        return counts[index];
    }

    public void skipTo(final int position) throws IOException {
        // skip to the specified position
        while (hasNextTransition()) {
            nextTransition();
            if (getPosition() >= position) {
                break;
            }
        }
    }

    public int getLength() {
        return lengths[index];
    }

    public void close() throws IOException {
        counts = null;
        lengths = null;
    }
}
