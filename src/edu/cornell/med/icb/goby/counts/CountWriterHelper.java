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

package edu.cornell.med.icb.goby.counts;

import java.io.Closeable;
import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Jun 12, 2009
 *         Time: 4:44:06 PM
 */
public class CountWriterHelper implements Closeable {
    private final CountsWriter delegate;
    private int previousPosition = -1;
    private int previousCount;
    private int lengthConstant = 1;
    private int previousPositionNotWritten;

    public CountWriterHelper(final CountsWriter delegate) {
        this.delegate = delegate;
    }

    public void appendCountAtPosition(final int count, final int position) throws IOException {
        System.out.printf("// count=%d position=%d previousCount=%d %n",
                count, position, previousCount);
        lengthConstant++;
        if (count == previousCount) {
            lengthConstant += position - previousPositionNotWritten;
            previousPositionNotWritten = position;
        } else {
            delegate.appendCount(previousCount, lengthConstant);
            previousCount = count;
            lengthConstant = 0;
            previousPosition = position;
            previousPositionNotWritten = position;
        }
    }

    public void close() throws IOException {
        //  if (previousCount != 0) delegate.appendCount(0, 1);
        delegate.close();
    }
}
