/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.counts;

import java.io.IOException;
import java.util.Iterator;

/**
 * Aggregate counts for peaks. Peaks are defined as a contiguous set of bases such that no base has zero count.
 * This class defines the peaks and aggregate counts for the peak.
 *
 * @author Fabien Campagne
 *         Date: May 27, 2009

 *         Time: 6:41:18 PM
 */
public class PeakAggregator implements Iterator<Peak>, Iterable<Peak> {
    CountsReaderI reader;
    private boolean nextLoaded;
    private int peakDetectionThreshold;

    /**
     * Will detect peaks in the counts information provided by reader.
     * @param reader
     */
    public PeakAggregator(final CountsReaderI reader) {
        this.reader = reader;
        currentPeak = new Peak();
    }


    private Peak currentPeak;

    /**
     * Returns true if another peak is detected above the detection threshold.
     *
     * @return
     * @throws java.io.IOException When an error occurs reading counts.
     */
    public final boolean hasNextPeak() throws IOException {
        if (nextLoaded) {
            return true;
        }
        currentPeak.start = -1;
        currentPeak.length = 0;
        currentPeak.count = 0;
        // find start of peak:

        while (reader.hasNextTransition()) {
            reader.nextTransition();

            final int baseCount = reader.getCount();
            if (baseCount > peakDetectionThreshold) {
                nextLoaded = true;
                //start of a new peak.
                currentPeak.count += baseCount;
                currentPeak.length += reader.getLength();
                break;
            }

        }
        currentPeak.start = reader.getPosition();

        // find end of peak:
        while (reader.hasNextTransition()) {
            reader.nextTransition();

            final int baseCount = reader.getCount();
            if (baseCount <= peakDetectionThreshold) {
                // past end of the peak.
                break;
            }
            currentPeak.length += reader.getLength();
            currentPeak.count += baseCount;

        }


        return nextLoaded;
    }

    /**
     * Return the next detected peak. The same instance of Peak is reused to provide information
     * about the peak. Make a copy if you need to save the peak for some purpose.
     * @return
     */
    public final Peak nextPeak() {

        if (hasNext()) {
            nextLoaded = false;
            return currentPeak;
        } else {
            throw new IllegalStateException();
        }
    }

    public boolean hasNext() {
        try {
            return hasNextPeak();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public Peak next() {

        return nextPeak();

    }

    public void remove() {
        throw new UnsupportedOperationException("This operation is not supported. Peaks are read-only.");
    }

    public Iterator<Peak> iterator() {
        return this;
    }

    /**
     * Set the threshold for peak detection. Zero by default.
     * @param peakDetectionThreshold The value of count that triggers the detection of a peak.
     */
    public void setPeakDetectionThreshold(final int peakDetectionThreshold) {
        this.peakDetectionThreshold = peakDetectionThreshold;
    }
}
