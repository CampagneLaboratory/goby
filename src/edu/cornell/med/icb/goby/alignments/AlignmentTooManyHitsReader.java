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

import com.google.protobuf.CodedInputStream;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

/**
 * Reads alignments too many hits data structure written with
 * {@link edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsWriter}.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 6:36:04 PM
 */
public class AlignmentTooManyHitsReader {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignmentTooManyHitsReader.class);

    /**
     * The underlying protocol buffer stream.
     */
    private final InputStream tooManyHitsStream;

    /**
     * A map from query index to the number of hits at that index.
     */
    private Int2IntMap queryIndex2NumHits = new Int2IntOpenHashMap();

    /**
     * A map from query index to depth/length of match.
     */
    private final Int2IntMap queryIndex2Depth = new Int2IntOpenHashMap();

    /**
     * The threshold used by the aligner to determine that a query is ambiguous and
     * should be dropped.
     */
    private int alignerThreshold;

    public AlignmentTooManyHitsReader(final String basename) throws IOException {
        final String filename = basename + ".tmh";
        final File optionalFile = new File(filename);
        this.queryIndex2NumHits = new Int2IntOpenHashMap();
        this.queryIndex2NumHits.defaultReturnValue(-1);
        this.queryIndex2Depth.defaultReturnValue(-1);

        InputStream tmhStream = null;
        if (optionalFile.exists()) {
            try {
                tmhStream = new GZIPInputStream(new FileInputStream(optionalFile));
            } catch (IOException e) {
                // try not compressed for compatibility with 1.6-:
                LOG.trace("falling back to legacy 1.6- uncompressed TMH.");

                tmhStream = new FileInputStream(optionalFile);
            }
            tooManyHitsStream = tmhStream;
            // accept very large too many hits messages, since these may describe more than 60 million reads:
            final CodedInputStream codedInput = CodedInputStream.newInstance(tooManyHitsStream);
            codedInput.setSizeLimit(Integer.MAX_VALUE);

            final Alignments.AlignmentTooManyHits tmh = Alignments.AlignmentTooManyHits.parseFrom(codedInput);
            for (final Alignments.AmbiguousLocation hit : tmh.getHitsList()) {
                queryIndex2NumHits.put(hit.getQueryIndex(), hit.getAtLeastNumberOfHits());
                if (hit.hasLengthOfMatch()) {
                    queryIndex2Depth.put(hit.getQueryIndex(), hit.getLengthOfMatch());
                }

            }
            this.alignerThreshold = tmh.getAlignerThreshold();
        } else {
            tooManyHitsStream = null;
            // the file does not exist. Log this fact, and act as if no query had too many hits.
            LOG.info("basename " + optionalFile + " has no 'too many hits' information ("
                    + basename + ".tmh does not exist)."
                    + " Assuming no queries have too many hits.");
        }
    }

    public AlignmentTooManyHitsReader(final InputStream tooManyHitsStream) {
        this.tooManyHitsStream = tooManyHitsStream;
    }

    /**
     * The number of hits against the reference that the aligner considered was too many to report.
     *
     * @return alignerThreshold.
     */
    public int getAlignerThreshold() {
        return alignerThreshold;
    }

    /**
     * Returns the number 'at least number of hits' reported by the alignment tool against
     * the reference.
     *
     * @param queryIndex The index of the query sequence.
     * @return The number of hits that triggered membership in the too many hits list.
     *         The query may hit more locations than reported here, since some alignment
     *         tools will just drop queries that match above a threshold and stop counting.
     *         This number can be >=k.
     */
    public final int getNumberOfHits(final int queryIndex) {
        return queryIndex2NumHits.get(queryIndex);
    }

    /**
     * Returns the length of match that resulted in the number of occurence againt the specific
     * reference sequence(s) this alignment was produced against.
     *
     * @param queryIndex The index of the query sequence.
     * @return The length of the longest match between the query and the reference sequence(s)
     *         that yielded the number of hits.
     */
    public final int getLengthOfMatch(final int queryIndex) {
        return queryIndex2Depth.get(queryIndex);
    }

    public final IntSet getQueryIndices() {
        return queryIndex2NumHits.keySet();
    }

    /**
     * Returns true if the query was considered ambiguous by the alignment tool.
     *
     * @param queryIndex The index of the query sequence.
     * @return True or false.
     */
    public boolean isQueryAmbiguous(final int queryIndex) {
        return queryIndex2NumHits.containsKey(queryIndex);
    }

    /**
     * Returns true if the query matched at least k number of locations in the reference.
     * However, if k >= alignerThreshold, then true is returned regardless.
     * <p/>
     * TODO : provide rationale for logic of k >= alignerThreshold behaviour
     * <p/>
     * TODO : discuss removal of this method as result is logically equivalent to isQueryAmbiguous(queryIndex)
     * TODO : given that tmhWriter.append() guarantees numHits > threshold
     *
     * @param queryIndex The index of the query sequence.
     * @param k          The parameter k.
     * @return True or false.
     */
    public final boolean isQueryAmbiguous(final int queryIndex, final int k) {
        final int atLeastNumberOfHits = queryIndex2NumHits.get(queryIndex);
        if (atLeastNumberOfHits == -1) {
            return false;
        }

        if (k >= alignerThreshold) {
            // since k is larger than the aligner threshold, we have to assume the query is
            // ambiguous at k, this is the safe choice.
            return true;
        } else {
            return (atLeastNumberOfHits >= k);
        }
    }

    /**
     * Returns true if the query matched at least k number of locations in the reference, at the
     * specified match length or less.
     * <p/>
     * TODO : discuss removal of 2nd argument from this method as result does not depend on its value
     * TODO : see comments at isQueryAmbiguous above
     *
     * @param queryIndex  The index of the query sequence.
     * @param k           The parameter k.
     * @param matchLength The match length.
     * @return True or false.
     */
    public final boolean isQueryAmbiguous(final int queryIndex, final int k, final int matchLength) {
        if (matchLength < getLengthOfMatch(queryIndex)) {
            return true;
        } else {
            return isQueryAmbiguous(queryIndex, k);
        }
    }
}
