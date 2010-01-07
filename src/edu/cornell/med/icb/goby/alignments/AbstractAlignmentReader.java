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

import edu.cornell.med.icb.identifier.IndexedIdentifier;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: May 20, 2009
 *         Time: 6:23:26 PM
 */
public abstract class AbstractAlignmentReader {
    protected int numberOfQueries;
    protected IndexedIdentifier queryIdentifiers;
    protected IndexedIdentifier targetIdentifiers;
    protected int numberOfTargets;
    protected int[] queryLengths;

    public int getNumberOfQueries() {
        assert isHeaderLoaded() : "Header must be loaded to access number of queries";
        return Math.max(numberOfQueries,0);
    }


    private boolean headerLoaded;

    protected boolean isHeaderLoaded() {
        return this.headerLoaded;
    }

    protected void setHeaderLoaded(final boolean headerLoaded) {
        this.headerLoaded = headerLoaded;
    }


    /**
     * Read the header of this alignment.
     *
     * @throws java.io.IOException If an error occurs.
     */
    public abstract void readHeader() throws IOException;

    /**
     * Return the query identifiers, if the header has been read. Null otherwise.
     *
     * @return query identifiers or null.
     * @see #readHeader()
     */
    public final IndexedIdentifier getQueryIdentifiers() {
        assert isHeaderLoaded() : "Header must be loaded to access query identifiers";

        return queryIdentifiers;
    }

    /**
     * Return the target identifiers, if the header has been read. Null otherwise.
     *
     * @return target identifiers or null.
     * @see #readHeader()
     */
    public final IndexedIdentifier getTargetIdentifiers() {
        assert isHeaderLoaded() : "Header must be loaded to access target identifiers";

        return targetIdentifiers;
    }

    public int getNumberOfTargets() {
        assert isHeaderLoaded() : "Header must be loaded to access number of targets";
        return Math.max(0, numberOfTargets);
    }

    /*
   Returns query lengths. An array of size the number of query sequences, where each element
   indicates the length of the query sequence.
    */
    public final int[] getQueryLength() {
        assert isHeaderLoaded() : "Header must be loaded to access query lengths";

        return queryLengths;
    }

    /**
     * Returns true if the input has more sequences.
     *
     * @return true if the input has more sequences, false otherwise.
     */
    public abstract boolean hasNextAligmentEntry();

    /**
     * Returns the next read entry from the input stream.
     *
     * @return the next read entry from the input stream.
     */
    public abstract Alignments.AlignmentEntry nextAlignmentEntry();
}
