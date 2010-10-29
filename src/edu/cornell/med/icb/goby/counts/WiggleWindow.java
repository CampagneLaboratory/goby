/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.PrintWriter;

/**
 * Class to assist with writing wiggles at a specific resolution (window size).
 * Averages the data within the specified window.
 *
 * @author Kevin Dorff
 */
public class WiggleWindow {
    /** Used to log debug and informational messages. */
    private static final Log LOG = LogFactory.getLog(WiggleWindow.class);

    /** The start position for the current window. */
    private int startPosition = -1;

    /** The position in the current window. */
    private int curWindowPosition;

    /** The writer to output wiggle data to. */
    private final PrintWriter writer;

    /** The window size to use. */
    private final int windowSize;

    /** The maximum data size, will not create a window with this value or larger. */
    private int maxDataSize;

    /** The running total for the current window. */
    private int windowTotal;

    /**
     * Create a new WiggleWindow with the given parameters.
     * @param writer The writer to output wiggle data to
     * @param windowSize The window size to use
     * @param maxDataSize The maximum data size, will not create a window with this value or larger
     */
    public WiggleWindow(final PrintWriter writer, final int windowSize, final int maxDataSize) {
        this.writer = writer;
        this.windowSize = windowSize;
        this.maxDataSize = maxDataSize;
    }

    /**
     * Add data to the output.
     * @param position the starting position for the data
     * @param length the length of the data with the given value
     * @param value the value of the data
     */
    public void addData(final int position, final int length, final int value) {
        if (value == 0) {
            return;
        }
        if (startPosition == -1) {
            // New window
            startPosition = position;
        } else  if (position < startPosition + windowSize) {
            // Fits in the current window
            curWindowPosition = position - startPosition;
        } else  {
            // Didn't fit in current window
            outputCurrentWindow();
            startPosition = position;
        }
        writeDataToWindow(length, value);
    }

    /**
     * Write the data to the window.
     * @param length the length of the data with the given value
     * @param value the value of the data
     */
    private void writeDataToWindow(final int length, final int value) {
        for (int i = 0; i < length; i++) {
            windowTotal += value;
            if (LOG.isTraceEnabled()) {
                LOG.trace(String.format("pos[%d]=%d", curWindowPosition + startPosition, value));
            }
            curWindowPosition++;
            if (curWindowPosition == windowSize) {
                outputCurrentWindow();
                startPosition += windowSize;
            }
        }
    }

    /**
     * Finish writing the data. You should call this after calling {@link #addData(int, int, int)}
     * for all of the data to finish writing the output.
     */
    public void finish() {
        outputCurrentWindow();
        if (writer != null) {
            writer.flush();
        }
    }

    /**
     * Reset for whole new run, ie, a new chromosome.
     */
    public void reset() {
        startPosition = -1;
        curWindowPosition = 0;
        windowTotal = 0;
    }

    /**
     * Get the max data size.
     * @return the max data size
     */
    public int getMaxDataSize() {
        return maxDataSize;
    }

    /**
     * Set the max data size.
     * @param maxDataSize the max data size
     */
    public void setMaxDataSize(final int maxDataSize) {
        this.maxDataSize = maxDataSize;
    }

    /**
     * Write a completed window of data.
     */
    private void outputCurrentWindow() {
        if (curWindowPosition == 0) {
            return;
        }
        final int windowAverage = Math.round((float) windowTotal / (float) curWindowPosition);
        if (LOG.isTraceEnabled()) {
            LOG.trace(String.format("... average=%d (%d/%d)",
                    windowAverage, windowTotal, curWindowPosition));
        }
        if ((startPosition + windowSize) <= maxDataSize) {
            // positions start at 1 for UCSC genome browser
            if (writer != null) {
                writer.printf("%d %d%n", startPosition + 1, windowAverage);
            }
        } else {
            LOG.warn(String.format("Not writing %d %d", startPosition + 1, windowAverage));
        }
        curWindowPosition = 0;
        windowTotal = 0;
    }
}
