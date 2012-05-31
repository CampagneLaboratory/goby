/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.counts;

import java.io.IOException;

/**
 * A helper class to facilitate writing counts with a CountsWriter.
 *
 * @author Fabien Campagne
 *         Date: Jun 12, 2009
 *         Time: 4:44:06 PM
 */
public class CountWriterHelper implements CountsWriterHelperI {
    private final CountsWriterI delegate;
    private int previousCount;
    private int accumulatedLength = 0;
    private int previousPosition = 0;
    private boolean firstAppend;

    public CountWriterHelper(final CountsWriterI delegate) {
        this.delegate = delegate;

        firstAppend = true;
    }

    @Override
    public void appendCountAtPosition(final int count, final int position) throws IOException {
        System.out.printf(" count=%d position=%d previousCount=%d length-constant=%d %n",
                count, position, previousCount, accumulatedLength);
        if (firstAppend && count != previousCount) {
            accumulatedLength = position;
            if (accumulatedLength != 0) {
                delegate.appendCount(0, accumulatedLength);
                accumulatedLength = 1;
                firstAppend = false;
                previousPosition = position;
                previousCount = count;
                return;
            }
        }
        if (firstAppend || (count != previousCount && previousPosition == position - 1)) {
            if (accumulatedLength > 0) {
                delegate.appendCount(previousCount, accumulatedLength);
                previousCount = count;
                previousPosition = position;
                accumulatedLength = 1;
            }
        } else {

            if (previousPosition == position - 1) {

                accumulatedLength++;
            } else {
                if (accumulatedLength > 0) {
                    //return to zero:
                    delegate.appendCount(previousCount, accumulatedLength);
                    accumulatedLength = position - previousPosition - 1;
                    previousCount = 0;
                    delegate.appendCount(0, accumulatedLength);
                    previousCount = count;
                    accumulatedLength = 1;
                }
            }
        }

        previousPosition = position;
        firstAppend = false;
    }

    @Override
    public void close() throws IOException {
        delegate.appendCount(previousCount, accumulatedLength);
        if (previousCount != 0) {
            delegate.appendCount(0, 1);
        }
        delegate.close();
    }
}
