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
import it.unimi.dsi.fastutil.ints.IntList;

import java.io.IOException;

/**
 * A helper class to facilitate writing counts with a CountsWriter.
 *
 * @author Fabien Campagne
 *         Date: Jun 12, 2009
 *         Time: 4:44:06 PM
 */
public class CountWriterHelper2 implements CountsWriterHelperI {
    private final CountsWriterI delegate;
    private int previousCount;
    private int previousPosition = -1;


    public CountWriterHelper2(final CountsWriterI delegate) {
        this.delegate = delegate;
        previousCount = delegate.getInitialCount();
        counts = new IntArrayList();
        positions = new IntArrayList();
        positions.add(-1);
        counts.add(0);
    }


    private final IntList counts;
    private final IntList positions;


    @Override
    public void appendCountAtPosition(final int count, final int position) throws IOException {
        /*   System.out.printf(" count=%d position=%d previousCount=%d%n",
                 count, position, previousCount);
        */

        if (position != previousPosition + 1 && count != 0) {
            // add transition that goes back to zero:
            counts.add(0);
            positions.add(position - 1);
        }
        counts.add(count);
        positions.add(position);


        if (positions.size() >= 2) {

            final int i = counts.size() - 2;
            final int j = counts.size() - 1;

            if (counts.getInt(i) == counts.getInt(j)) {

                positions.removeElements(i, j);
                counts.removeElements(i, j);


            }

        }
        while (positions.size() > 2) {
            final int diffPos = positions.getInt(1) - (positions.getInt(0) + 1);

            int aCount = counts.getInt(1);
            if (aCount < 0) {
                System.out.printf("Count can never be negative (found value=%d at position %d). Setting count to zero", aCount, positions.getInt(0));
                aCount = 0;
            }
            delegate.appendCount(aCount, diffPos + 1);
            positions.removeElements(0, 1);
            counts.removeElements(0, 1);
            previousCount = count;

        }
        previousPosition = position;

    }

    @Override
    public void close() throws IOException {
        int n = counts.size() - 1;
        if (counts.get(n) != 0) {
            // return the count to zero:
            counts.add(0);
            positions.add(positions.get(n) + 1);
        }
        final int count = counts.get(0);
        // write the last counts:
        while (positions.size() >= 2) {
            final int diffPos = positions.getInt(1) - (positions.getInt(0) + 1);
            int aCount = counts.getInt(1);
            if (aCount < 0) {
                System.out.printf("Count can never be negative (found value=%d at position %d). Setting count to zero", aCount, positions.getInt(0));
                aCount = 0;
            }
            delegate.appendCount(aCount, diffPos + 1);
            positions.removeElements(0, 1);
            counts.removeElements(0, 1);
            previousCount = count;
        }

        delegate.close();
    }
}
