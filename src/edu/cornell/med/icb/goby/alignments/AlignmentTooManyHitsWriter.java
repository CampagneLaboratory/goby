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

package edu.cornell.med.icb.goby.alignments;

import java.io.Closeable;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;

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
    private final OutputStream tooManyHitsOutput;
    private final Alignments.AlignmentTooManyHits.Builder tooManyHits;
    private Alignments.AmbiguousLocation.Builder newAmbiguousLocation;


    public AlignmentTooManyHitsWriter(final String outputBasename, final int alignerThreshold) throws IOException {
        tooManyHitsOutput = new GZIPOutputStream(new FileOutputStream(outputBasename + ".tmh"));
        newAmbiguousLocation = Alignments.AmbiguousLocation.newBuilder();
        tooManyHits = Alignments.AlignmentTooManyHits.newBuilder();
        tooManyHits.setAlignerThreshold(alignerThreshold);
    }

    /**
     * Update the aligner threshold.
     * @param alignerThreshold the new threshold to write in the too many hits file.
     */
    public void setAlignerThreshold(final int alignerThreshold) {
        tooManyHits.setAlignerThreshold(alignerThreshold);
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

    public Alignments.AmbiguousLocation.Builder getNewAmbiguousLocation() {
        return newAmbiguousLocation;
    }

    /**
     * Append record defined by 3 arguments.
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
     */
    public void append() {
        assert tooManyHits.hasAlignerThreshold() : "append> writer missing aligner threshold";
        assert newAmbiguousLocation.hasAtLeastNumberOfHits() : "append> new record missing atLeastNumberOfHits";
        if (newAmbiguousLocation.getAtLeastNumberOfHits() > tooManyHits.getAlignerThreshold()) {
            tooManyHits.addHits(newAmbiguousLocation.build());
        }
        //tooManyHits.addHits(newAmbiguousLocation.build());
        // whether or not the hit was added, reset/create a new one
        newAmbiguousLocation = Alignments.AmbiguousLocation.newBuilder();
    }
}
