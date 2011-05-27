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

import java.io.IOException;

/**
 * Skeleton for algorithm development discussion.
 * @author Fabien Campagne
 *         Date: 5/26/11
 *         Time: 10:14 PM
 */
public class UnionAlgorithmSkeleton  implements CountsReaderI {
    private int numReaders;
    private CountsReaderI[] readers;

    public UnionAlgorithmSkeleton(CountsReaderI ... readers) {
        this.numReaders=readers.length;
        this.readers=readers;
    }

    public int getPosition() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean hasNextTransition() throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void nextTransition() throws IOException {
        //To change body of implemented methods use File | Settings | File Templates.
    }


    public void skipTo(int position) throws IOException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getLength() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void close() throws IOException {
        //To change body of implemented methods use File | Settings | File Templates.
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

    public final CountsReaderI[] getReaders() {
        return readers;
    }

    /**
     * Return the count for a specific reader.
     * @param readerIndex  Index ((zero-based) of the reader when provided as parameter to the constructor
     * @return count for the reader identified by readerIndex.
     */
    public final int getCount(final int readerIndex) {
      return 0;
      //To change body of implemented methods use File | Settings | File Templates.
    }
}
