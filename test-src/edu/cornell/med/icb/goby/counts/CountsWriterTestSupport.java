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

import it.unimi.dsi.lang.MutableString;

import java.io.IOException;

/**
 * CountsWriter useful in JUnit tests.
 *
 * @author Fabien Campagne
 *         Date: 10/29/11
 *         Time: 1:45 PM
 */
public class CountsWriterTestSupport implements CountsWriterI {
    private int numTransitions;
    private final MutableString transitionsAsText = new MutableString();
    private int previousCount;
    private int initialCount = 0;

    /**
     * Initial count is set to the value specified.
     *
     * @param initialCount count at position zero.
     */
    public CountsWriterTestSupport(final int initialCount) {
        this.initialCount = initialCount;
    }

    /**
     * Initial count is set to zero.
     */
    public CountsWriterTestSupport() {
        this.initialCount = 0;
        previousCount = initialCount;
    }

    @Override
    public long getNumberOfBitsWritten() {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    @Override
    public int getNumberOfTransitions() {
        return numTransitions;
    }

    @Override
    public void appendCount(final int count, final int lengthConstant) throws IOException {
        assert count != previousCount : " count must be different to start a new transition.";
        previousCount = count;
        numTransitions++;
        System.out.println(String.format(" appending (count=%d,length=%d)", count, lengthConstant));
        transitionsAsText.append(String.format("(c=%d,l=%d)", count, lengthConstant));


    }

    /**
     * Returns transitions in the format [(length,count)]+
     *
     * @return transitions in the format [(length,count)]+
     */
    public String countsAsText() {
        return "initial-count=" + initialCount + " " + transitionsAsText.toString();
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public long getNumberOfBasesSeen() {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    @Override
    public long getNumberOfSitesSeen() {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    @Override
    public int getInitialCount() {

        return initialCount;
    }
}
