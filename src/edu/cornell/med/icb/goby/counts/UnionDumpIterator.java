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

/**
 * A naive implementation of the union iterator. Beware, this is slow and uses a lot of memory (proportional to the
 * length of the sequence).
 *
 * @author Fabien Campagne
 *         Date: May 23, 2011
 *         Time: 4:25:02 PM
 */
public class UnionDumpIterator implements CountsAggregatorI {
    private CountsReaderI[] readers;
    IntArrayList positions[];
    IntArrayList counts[];
    IntArrayList lengths[];
    int sizes[];
    int numReaders;
    private int position = 0;
    private boolean hasNextTransition;
    private int index;
    private int count[];
    private int length;
    private int maxSize;
    private int beforePosition;
    private boolean doneOnNext;


    public UnionDumpIterator(final CountsReaderI... countReader) throws IOException {
        readers = countReader;
        positions = new IntArrayList[readers.length];
        counts = new IntArrayList[readers.length];
        sizes = new int[readers.length];
        int i = 0;
        numReaders = readers.length;
        count = new int[numReaders];
        for (CountsReaderI reader : readers) {
            counts[i] = new IntArrayList();
            positions[i] = new IntArrayList();
            int position = -1;
            while (reader.hasNextTransition()) {

                reader.nextTransition();
                final int count = reader.getCount();
                final int length = reader.getLength();
                int pos = reader.getPosition();
                //    for (int j = position; j < pos; ++j) {
                //        counts[i].add(0);
                //positions[i].add(j);
                //  }
                for (int j = 0; j < length; j++) {

                    counts[i].add(count);
                    position++;
                    // positions[i].add(position);
                }

            }
            sizes[i] = counts[i].size();
            maxSize = Math.max(sizes[i], maxSize);
            i++;
        }

    }

    public int getPosition() {
        return beforePosition;
    }

    public boolean hasNextTransition() throws IOException {
        if (hasNextTransition) {
            return true;
        }
        if (position > maxSize) {
            beforePosition = position - 1;
            return false;
        }

        beforePosition = position;
        while (sameCounts(position, position + 1) && position < maxSize) {
            ++position;
        }

        length = position - beforePosition + 1;
        // System.out.printf("Position %d length=%d%n",position, length);

        for (int i = 0; i < numReaders; i++) {
            final int count = beforePosition >= sizes[i] ? 0 : counts[i].getInt(beforePosition);
            this.count[i] = count;
            //   System.out.printf("Position %d count[%d]=%d%n",position, i,count);
        }

        position++;
        if (position >= maxSize) {
            length -= 1;
        }
        // position+=length-1;

        hasNextTransition = true;
        return true;
    }

    private boolean sameCounts(int position, int nextPos) {
        if (position == -1) return true;

        for (int i = 0; i < numReaders; i++) {


            if (nextPos < sizes[i] && counts[i].getInt(position) != counts[i].getInt(nextPos)) return false;
            if (nextPos==sizes[i]
                    && counts[i].getInt(position) !=0 && nextPos!=maxSize) {
                // the reader ends in a count different from zero, we have to detect that transition.
                return false;
           }
        }

        return true;
    }

    public final int getCount(final int readerIndex) {

        return count[readerIndex];
        //if (readerIndex == 1) {
        //    System.out.printf("position %d annotation count=%d%n", position, count);
        // }


    }

    public void nextTransition() throws IOException {
        if (!hasNextTransition()) throw new IllegalStateException("No such element");
        hasNextTransition = false;
    }

    /**
     * Return the sum of counts over the readers that have non zero counts at the current position.
     */
    public int getCount() {
        int count = 0;
        for (int i = 0; i < numReaders; i++) {
            count += getCount(i);
        }
        return count;
    }

    public void skipTo(int position) throws IOException {
        for (CountsReaderI reader : readers) {
            reader.skipTo(position);
        }
    }

    @Override
    public void reposition(int position) throws IOException {
         throw new UnsupportedOperationException("this implementation does not support this method.");
    }

    public int getLength() {
        return length;
    }

    public void close() throws IOException {
        int index = 0;
        for (CountsReaderI reader : readers) {
            counts[index].clear();
            positions[index].clear();
            reader.close();
            index++;
        }
    }
}
