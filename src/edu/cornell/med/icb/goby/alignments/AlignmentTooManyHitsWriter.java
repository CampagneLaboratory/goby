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

package edu.cornell.med.icb.goby.alignments;

import java.io.Closeable;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * This class write the too many hit datastructure Protocol Buffer format.
 * See Alignements.proto for the specification of this format.
 *
 * @author Fabien Campagne
 *         Date: May 5, 2009
 *         Time: 2:03 PM
 */
public class AlignmentTooManyHitsWriter implements Closeable {
    private boolean tooManyHitsWritten;
    private FileOutputStream tooManyHitsOutput;
    private Alignments.AlignmentTooManyHits.Builder tooManyHits;
    private Alignments.AmbiguousLocation.Builder newAmbiguousLocation;


    public AlignmentTooManyHitsWriter(final String outputBasename, final int alignerThreshjold) throws FileNotFoundException {
        tooManyHitsOutput = new FileOutputStream(outputBasename + ".tmh");
        newAmbiguousLocation = Alignments.AmbiguousLocation.newBuilder();
        tooManyHits = Alignments.AlignmentTooManyHits.newBuilder();
        tooManyHits.setAlignerThreshold(alignerThreshjold);

    }

    /**
     * Update the aligner threshold.
     * @param alignerThreshjold the new threshold to write in the too many hits file.
     */
    public void setAlignerThreshold(final int alignerThreshjold) {
        tooManyHits.setAlignerThreshold(alignerThreshjold);
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        write();
    }

    public void write() throws IOException {
        if (!tooManyHitsWritten) {
            tooManyHits.build().writeTo(tooManyHitsOutput);
            tooManyHitsOutput.close();
            tooManyHitsWritten = true;
        }
    }

    /**
     * Obtain the too many hits set that is being prepared. Set values on the set, then call appendTooManyHits()
     *
     * @return the current alignment entry.
     */
    public Alignments.AmbiguousLocation.Builder getNewAmbiguousLocation() {
        return newAmbiguousLocation;
    }

    /**
     * Append record defined by 3 arguments
     * previously called appendTooManyHits()
     */
    public void append(final int queryIndex, final int howManyHits, final int lengthOfMatch) {
        newAmbiguousLocation.setQueryIndex(queryIndex);
        newAmbiguousLocation.setAtLeastNumberOfHits(howManyHits);
        newAmbiguousLocation.setLengthOfMatch(lengthOfMatch);
        append();
    }

    /**
     * Append the current too many hits record.
     * - previously appendTooManyHits()
     */
    public void append() {
        assert (tooManyHits.hasAlignerThreshold()) : "append> writer missing aligner threshold";
        assert (newAmbiguousLocation.hasAtLeastNumberOfHits()) : "append> new record missing atLeastNumberOfHits";
        if (newAmbiguousLocation.getAtLeastNumberOfHits() > tooManyHits.getAlignerThreshold()) {
            tooManyHits.addHits(newAmbiguousLocation.build());
        }
        //tooManyHits.addHits(newAmbiguousLocation.build());
        // whether or not the hit was added, reset/create a new one
        newAmbiguousLocation = Alignments.AmbiguousLocation.newBuilder();
    }
}
