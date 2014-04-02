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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.modes.UpgradeMode;
import edu.cornell.med.icb.goby.util.FileExtensionHelper;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.lang.StringUtils;

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
     * Length of each target sequence.
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
     * The length of all the query sequences. Only valid if {@link #constantQueryLengths} is true.
     */
    protected int constantLength;
    /**
     * The smallest possible query index in this alignment. Data stored as an array where
     * queryIndex is the array index will be stored with only the elements in the inclusive
     * range [smallestSplitQueryIndex largestSplitQueryIndex]
     * Such data structures include queryLength and some arrays in the TooManyHits data
     * structure.
     */
    protected int smallestQueryIndex;

    /**
     * This constructor will upgrade the alignment to the latest version of the Goby data structure.
     *
     * @param upgrade
     * @param basename Alignment basename.
     */
    protected AbstractAlignmentReader(boolean upgrade, String basename) {
        // do not upgrade when the upgrade flag is false;
        if (basename != null && upgrade) {
            // we can only upgrade alignments for which we have a basename
            UpgradeMode upgrader = new UpgradeMode();
            upgrader.setSilent(true);
            upgrader.upgrade(basename);
        }
    }

    /**
     * Return the basename corresponding to the input alignment filename.  Note
     * that if the filename does have the extension known to be a compact alignment
     * the returned value is the original filename
     *
     * @param filename The name of the file to get the basename for
     * @return basename for the alignment file
     */
    public static String getBasename(final String filename) {
        for (final String ext : FileExtensionHelper.COMPACT_ALIGNMENT_FILE_EXTS) {
            if (StringUtils.endsWith(filename, ext)) {
                return StringUtils.removeEnd(filename, ext);
            }
        }

        // perhaps the input was a basename already.
        return filename;
    }

    /**
     * Return the basenames corresponding to the input filenames. Less basename than filenames
     * may be returned (if several filenames reduce to the same baseline after removing
     * the extension). This method returns the unique set of basenames in the same order they are
     * provided as argument.
     *
     * @param filenames The names of the files to get the basnames for
     * @return An array of basenames
     */
    public static String[] getBasenames(final String... filenames) {
        final ObjectSet<String> result = new ObjectArraySet<String>();
        final ObjectList<String> unique = new ObjectArrayList<String>();
        if (filenames != null) {
            for (final String filename : filenames) {


                final String newBasename = getBasename(filename);
                if (!result.contains(newBasename)) {
                    unique.add(newBasename);
                    result.add(getBasename(filename));
                }
            }
        }
        return unique.toArray(new String[unique.size()]);
    }

    /**
     * The smallest possible query index in this alignment.
     *
     * @return The smallest possible query index in this alignment.
     */
    public int getSmallestSplitQueryIndex() {
        return smallestQueryIndex;
    }

    /**
     * The largest possible query index in this alignment.
     *
     * @return The largest possible query index in this alignment.
     */
    public int getLargestSplitQueryIndex() {
        return largestQueryIndex;
    }

    /**
     * The largest possible query index in this alignment. Data stored as an array where
     * queryIndex is the array index will be stored with only the elements in the inclusive
     * range [smallestSplitQueryIndex largestSplitQueryIndex]
     * Such data structures include queryLength and some arrays in the TooManyHits data
     * structure.
     */
    protected int largestQueryIndex;

    /**
     * Get the number of query sequences represented in this alignment.
     *
     * @return the number of query sequences represented in this alignment.
     */
    public int getNumberOfQueries() {
        assert isHeaderLoaded() : "Header must be loaded to access number of queries";
        return Math.max(0,Math.max(numberOfQueries, largestQueryIndex-smallestQueryIndex));
    }

    /**
     * Has the alignment header has been processed?
     *
     * @return true if the reader has loaded the header
     */
    protected final boolean  isHeaderLoaded() {
        return this.headerLoaded;
    }

    /**
     * Indicate that the alignment header has been processed.
     *
     * @param headerLoaded whether or not the reader has loaded the header
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
     * Get the number of target sequences represented in this alignment.
     *
     * @return the number of target sequences represented in this alignment.
     */
    public int getNumberOfTargets() {
        assert isHeaderLoaded() : "Header must be loaded to access number of targets";
        return Math.max(0, numberOfTargets);
    }


    /**
     * Returns the length of a target.
     *
     * @param targetIndex Index of the target sequence.
     * @return Length of the specified target sequence.
     */
    public final int getTargetLength(final int targetIndex) {
        assert isHeaderLoaded() : "Header must be loaded to access target lengths";

        assert targetLengths != null : "Target lengths must exist in the header.";
        return targetLengths[targetIndex];
    }

    /**
     * Returns target lengths. An array of size the number of target sequences, where each element
     * indicates the length of the target sequence.
     *
     * @return an array containing the lengths of all the targets represented in the alignment
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


    protected String alignerName;

    public String getAlignerName() {
        return alignerName;
    }

    public String getAlignerVersion() {
        return alignerVersion;
    }

    /**
     * the list of aligner versions in the order of names, with duplicates removed.
     */
    protected String alignerVersion;

    protected ObjectArrayList<String> sampleBasenames=new ObjectArrayList<String>();

    /**
     * Return the sample basename associated with the given sampleIndex.
     * @param sampleIndex index of the sample, from AlignmentEntry.
     * @return The sample basename/identifier for the index.
     */
    public String getSampleBasename(final int sampleIndex) {
        return sampleBasenames.get(sampleIndex);
    }
    /**
     * Return the minimum location in this alignment.
     *
     * @return The minimum location in this alignment.
     * @throws IOException
     */
    public abstract  ReferenceLocation getMinLocation() throws IOException;

    /**
     * Return the maximum location in this alignment.
     *
     * @return The maximum location in this alignment.
     * @throws IOException
     */
    public abstract ReferenceLocation getMaxLocation() throws IOException;
}
