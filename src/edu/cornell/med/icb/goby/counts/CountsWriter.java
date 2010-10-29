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

import it.unimi.dsi.io.OutputBitStream;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Writes counts in a compressed format.
 *
 * @author Fabien Campagne
 *         Date: May 6, 2009
 *         Time: 2:43:45 PM
 */
public class CountsWriter implements Closeable {
    private static final Log LOG = LogFactory.getLog(CountsWriter.class);
    private OutputBitStream out;
    private int previousCount;
    private boolean dataAlreadyWritten;
    private int numberOfCountsWritten;
    private int bitsWritten;
    private int position;

    @Deprecated
    public CountsWriter(final OutputStream output, final int initialCount) throws IOException {
        out = new OutputBitStream(output);
        setInitialCount(initialCount);
    }

    @SuppressWarnings("deprecation")
    public CountsWriter(final OutputStream output) throws IOException {
        this(output, 1);
    }

    /**
     * Set the count at position zero.  Always set initial count before you append counts
     * (see appendCount() method).
     *
     * @param count The count to use for the initial value.
     * @throws IOException if there is an problem writing the counts
     */
    private void setInitialCount(final int count) throws IOException {
        if (dataAlreadyWritten) {
            throw new IllegalStateException("Data must not have been written before setting initial count.");
        }

        previousCount = count;
        bitsWritten += out.writeDelta(previousCount + 1);  // Delta cannot be zero, so add 1.
    }

    public void appendCount(final int count, final int lengthConstant) throws IOException {
        assert lengthConstant > 0 : "length must be greater than zero.";
        //  System.out.printf("appending %d (%d, %d)  %n", position, count, lengthConstant);
        final int deltaCount = count - previousCount;
        final int deltaCountEncoded;
        deltaCountEncoded = encodeDeltaCount(deltaCount);

        assert deltaCountEncoded > 0 : " delta count integer must not be zero";
        bitsWritten += out.writeGamma(deltaCountEncoded);
        bitsWritten += out.writeGamma(lengthConstant);
        dataAlreadyWritten = true;
        //  if (LOG.isTraceEnabled()) {
        //     LOG.trace(deltaCountEncoded + ", " + lengthConstant + " written "
        //            + bitsWritten + " bits");
        // }
        previousCount = count;
        ++numberOfCountsWritten;
        position += lengthConstant;
    }

    protected static int encodeDeltaCount(final int deltaCount) {
        final int deltaCountEncoded;
        if (deltaCount < 0) {
            // odd numbers encode negative integers:
            deltaCountEncoded = -deltaCount * 2 + 1;
        } else {
            // even numbers encode positive integers:
            deltaCountEncoded = deltaCount * 2;
        }
        return deltaCountEncoded;
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        if (out != null) {
            bitsWritten += out.writeGamma(CountsReader.END_OF_DATA_MARKER);

            out.flush();
            out.close();
            out = null;
            if (LOG.isInfoEnabled()) {
                LOG.info("bits written: " + bitsWritten);
                LOG.info("bytes written: " + bitsWritten / 8);
                LOG.info("bits/count_transition (average): " + (float) bitsWritten / (float) numberOfCountsWritten);
            }
        }
    }
}
