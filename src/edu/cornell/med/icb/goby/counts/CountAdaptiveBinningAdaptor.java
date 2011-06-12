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
import java.util.NoSuchElementException;

/**
 * An adapter over a CountReaderI which returns average counts over bins. Bins are determined adaptively. A bin is
 * defined over count peaks whose value is above zero. This adapter is used in IGV to create lower
 * resolution views of the count data. The bin size can be chosen to return just a few bins over an entire chromosome.
 *
 * @author Fabien Campagne
 *         Date: 6/12/11
 *         Time: 12:31 PM
 */
public class CountAdaptiveBinningAdaptor implements CountBinningAdapterI {
    final CountsReaderI delegate;
    int position = -1;
    int length;
    /**
     * The average count over the bin, counting only site where count>0.
     */
    double average;
    /**
     * The maximum count seen in the bin.
     */
    int max;
    long sumBasesOverBin;
    private boolean binLoaded;
    /**
     * The next transition if available was encountered at the end of the last bin, just outside its right limit.
     */
    private int nextCount;
    private int nextLength;
    private int nextPosition;
    private boolean haveCachedNextTransition;


    public CountAdaptiveBinningAdaptor(final CountsReaderI delegate) {
        this.delegate = delegate;

    }

    @Override
    public void close() throws IOException {
        delegate.close();
    }

    @Override
    public int getPosition() {
        return position;
    }

    boolean previousCountWasZero = true;

    @Override
    public boolean hasNextTransition() throws IOException {
        if (binLoaded) {
            return true;
        } else {

            if (!haveCachedNextTransition && !delegate.hasNextTransition()) {
                // no more transitions in the source, we are done already.
                return false;
            }
            length = 0;
            sumBasesOverBin = 0;
            max = 0;
            position = Integer.MAX_VALUE;


            boolean binCompleted = false;

            if (haveCachedNextTransition) {
                length = nextLength;
                sumBasesOverBin = (long) nextCount * nextLength;
                max = nextCount;
                haveCachedNextTransition = false;
                position = Math.min(nextPosition, position);
                previousCountWasZero=nextCount==0;
            }

            while (delegate.hasNextTransition() && !binCompleted) {
                delegate.nextTransition();
                // the length of the bin transition is the position of the current position minus the position of the
                // very first non-zero count observed in the bin, plus the length of the current transition:
                final int delegateLength = delegate.getLength();
                final int count = delegate.getCount();
                if (count == 0 && previousCountWasZero || count != 0 && !previousCountWasZero) {
                    // the count is still in the current bin (count!0 peak or count=0 stretch):
                    length += delegateLength;

                   sumBasesOverBin += (long) count * delegateLength;
                    max = Math.max(count, max);
                    //set position to the first non-zero count encountered in a bin:

                    position = Math.min(delegate.getPosition(), position);
                    previousCountWasZero = count == 0;
                } else {
                    binCompleted = true;
                    nextCount = count;
                    nextLength = delegateLength;
                    nextPosition = delegate.getPosition();
                    haveCachedNextTransition = true;

                }
            }
            // set count to average count over all counts with count>0 observed within bin size:
            average = ((double) sumBasesOverBin) / (double) length;
        }
        binLoaded = true;
        return true;
    }

    @Override
    public void nextTransition() throws IOException {
        if (!hasNextTransition()) {
            throw new NoSuchElementException("No such element.");
        }
        binLoaded = false;


    }

    @Override
    public int getCount() {
        return (int) average;
    }


    @Override
    public double getAverage() {
        return average;
    }

    /**
     * Advance up to or past the specified position. The reader is advanced until the position returned by getPosition()
     * is at least equal, or greater to the specified position.
     *
     * @param position position to skip to.
     * @throws java.io.IOException
     */
    @Override
    public void skipTo(final int position) throws IOException {
        // skip to the specified position
        while (hasNextTransition() && getPosition() < position) {
            nextTransition();

        }
    }

    @Override
    public int getLength() {
        return length;
    }

    /**
     * Get the maximum count observed over the bin.
     *
     * @return return the maximum count in the current bin.
     */
    @Override
    public int getMax() {
        return max;
    }
}
