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

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;


/**
 * Abstract class for reading Goby compact alignments.
 *
 * @author Fabien Campagne
 *         Date: May 20, 2009
 *         Time: 6:23:26 PM
 */
public abstract class AbstractAlignmentReader implements Closeable,
        Iterator<Alignments.AlignmentEntry>, Iterable<Alignments.AlignmentEntry> {
    /**
     * Mapping from query indicies to query identifier strings.  Note that it is
     * not necessary for every query to have an associated identifer.
     *
     * @see #numberOfQueries
     */
    protected IndexedIdentifier queryIdentifiers;

    /**
     * The number of query sequences represented in this alignment.
     */
    protected int numberOfQueries;

    /**
     * Length of each query sequence, or null when constantQueryLengths is true.
     */
    protected int[] queryLengths;

    /**
     * Mapping from target indicies to target identifier strings.  Note that it is
     * not necessary for every target to have an associated identifer.
     *
     * @see #numberOfTargets
     */
    protected IndexedIdentifier targetIdentifiers;

    /**
     * The number of target sequences represented in this alignment.
     */
    protected int numberOfTargets;

    /**
     * Length of each target sequence, or null when constantQueryLengths==true.
     */
    protected int[] targetLengths;

    /**
     * Indicates that the alignment header has been processed.
     */
    private boolean headerLoaded;
    /**
     * True if all the query sequences have the same lengths.
     */
    protected boolean constantQueryLengths;
    /**
     * The length of all the query sequences. Only valid if constantQueryLengths==true.
     */
    protected int constantLength;

    /**
     * Get the number of query sequences represented in this alignment.
     * @return The number of query sequences in this alignment
     */
    public int getNumberOfQueries() {
        assert isHeaderLoaded() : "Header must be loaded to access number of queries";
        return Math.max(numberOfQueries, 0);
    }

    /**
     * Determine whether or not the alignment header informaton has been processed.
     * @return true if the header has been loaded
     */
    protected boolean isHeaderLoaded() {
        return this.headerLoaded;
    }

    /**
     * Indicate whether or not the alignment header informaton has been processed.
     * @param headerLoaded whether or not the header has been loaded
     */
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

    /**
     * Get the number of target sequences used in this alignment.
     * @return The number of target sequences in this alignment
     */
    public int getNumberOfTargets() {
        assert isHeaderLoaded() : "Header must be loaded to access number of targets";
        return Math.max(0, numberOfTargets);
    }

    /**
     * Returns query lengths. An array of size the number of query sequences, where each element
     * indicates the length of the query sequence.
     * @return an array of target lengths
     */
    public final int[] getQueryLengths() {
        assert isHeaderLoaded() : "Header must be loaded to access query lengths";
        if (constantQueryLengths) {
            final int[] localQueryLengths = new int[numberOfQueries];
            for (int i = 0; i < localQueryLengths.length; ++i) {
                localQueryLengths[i] = constantLength;

            }
            return localQueryLengths;
        } else {
            return queryLengths;
        }
    }

    /**
     * Returns the length of a query.
     * @param queryIndex Index of the query sequence.
     * @return Length of the specified query sequence.
     */
    public final int getQueryLength(final int queryIndex) {
        assert isHeaderLoaded() : "Header must be loaded to access query lengths";
        if (constantQueryLengths) {
            return constantLength;
        } else {
            assert queryLengths != null : "Query lengths must exist in the header.";
            return queryLengths[queryIndex];
        }
    }

    /**
     * Returns the length of a target.
     * @param targetIndex Index of the target sequence.
     * @return Length of the specified target sequence.
     */
    public final int getTargetLength(final int targetIndex) {
        assert isHeaderLoaded() : "Header must be loaded to access target lengths";

        assert targetLengths != null : "Target lengths must exist in the header.";
        return targetLengths[targetIndex];
    }



    /**
     * Returns target lengths. An array of size the number of target sequences, where each
     * element indicates the length of the target sequence.
     * @return an array of target lengths
     */
    public final int[] getTargetLength() {
        assert isHeaderLoaded() : "Header must be loaded to access target lengths";
        return targetLengths;
    }

    /**
     * @return True if the alignment stores a constant query length.
     */
    public boolean isConstantQueryLengths() {
        return constantQueryLengths;
    }
}
