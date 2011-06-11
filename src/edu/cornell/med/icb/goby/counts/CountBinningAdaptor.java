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
 * An adapter over a CountReaderI which returns average counts over bins. This adapter is used in IGV to create lower
 * resolution views of the count data. The bin size can be chosen to return just a few bins over an entire chromosome.
 *
 * @author Fabien Campagne
 *         Date: 6/11/11
 *         Time: 12:31 PM
 */
public class CountBinningAdaptor implements CountsReaderI {
    final CountsReaderI delegate;
    private final int binSize;
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


    public CountBinningAdaptor(final CountsReaderI delegate, final int binSize) {
        this.delegate = delegate;
        this.binSize = binSize;
    }

    @Override
    public void close() throws IOException {
        delegate.close();
    }

    @Override
    public int getPosition() {
        return position;
    }

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
            int numSitesObserved = 0;
            max = 0;
            position = Integer.MAX_VALUE;

            if (haveCachedNextTransition) {
                length = nextLength;
                sumBasesOverBin = (long) nextCount * nextLength;
                max = nextCount;
                haveCachedNextTransition = false;
                if (nextCount != 0) {
                    position = Math.min(nextPosition, position);
                    numSitesObserved += nextLength;
                }
            }
            while (delegate.hasNextTransition() && numSitesObserved < binSize) {
                delegate.nextTransition();
                // the length of the bin transition is the position of the current position minus the position of the
                // very first non-zero count observed in the bin, plus the length of the current transition:
                final int delegateLength = delegate.getLength();
                final int newLength = Math.max(delegateLength, delegate.getPosition() - position + delegateLength);
                final int count = delegate.getCount();
        //        if (count==0 || newLength < binSize) {
                    // the count is still in the current bin:
                    length = newLength;

                    if (count != 0) {
                        numSitesObserved += delegateLength;
                    }
                    sumBasesOverBin += (long) count * delegateLength;
                    max = Math.max(count, max);
                    //set position to the first non-zero count encountered in a bin:
                    if (count != 0) {
                        position = Math.min(delegate.getPosition(), position);
                    }

            /*    } else {
                    // we need to store the next transition for latter use:
                    nextCount = delegate.getCount();
                    nextLength = delegate.getLength();
                    nextPosition = delegate.getPosition();
                    haveCachedNextTransition = true;
                    return true;
                }    */
            }
            // set count to average count over all counts with count>0 observed within bin size:
            average = ((double) sumBasesOverBin) / (double) numSitesObserved;
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


    public double getAverage() {
        return average;
    }

    /**
     * Advance up to or past the specified position. The reader is advanced until the position returned by getPosition()
     * is at least equal, or greater to the specified position.
     *
     * @param position position to skip to.
     * @throws IOException
     */
    @Override
    public void skipTo(final int position) throws IOException {
        // skip to the specified position
        while (hasNextTransition() && getPosition()< position) {
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
    public int getMax() {
        return max;
    }
}
