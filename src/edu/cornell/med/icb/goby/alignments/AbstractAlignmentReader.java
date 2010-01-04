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
        return Math.max(0,numberOfTargets);
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
