/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.*;
import org.junit.Test;
import org.junit.BeforeClass;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNotNull;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.io.File;
import java.util.Map;
import java.util.Arrays;

import edu.cornell.med.icb.io.TsvToFromMap;
import edu.cornell.med.icb.maps.LinkedHashToMultiTypeMap;
import edu.cornell.med.icb.iterators.TsvLineIterator;


/**
 * Test iterated sorted alignments with a set of sequence variations
 * that were hand coded for alignment with a synthetic, known reference.
 * This will use the .compact-reads file and the alignment as input.
 * The compact-reads file is used to determine the queryLength of
 * the input reads.
 *
 * The file "seq-var-reads-gsnap.display-seq-var-tsv-base.tsv" is output
 * from DisplaySequenceVariations in the TAB_SINGLE_BASE format. It
 * has been hand checked for accuracy.
 *
 * The the MyIterateSortedAlignments class below will be used to read
 * the sequence variations to compare against the variations found
 * in "seq-var-reads-gsnap.display-seq-var-tsv-base.tsv".
 * User: kdorff
 */
public class TestIteratedSortedAlignment2 {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(TestIteratedSortedAlignment2.class);

    @BeforeClass
    public static void beforeClass() throws IOException {
    }

    @AfterClass
    public static void afterClass() throws IOException {

    }

    // TODO: test the refPositions in an insertion that the refPositions within the insertion
    // TODO: are constant AND that it matches the refPosition just before the insertion, if
    // TODO: there is one.
    //
    // TODO: test the readIndex in an deletion that the readIndexes within the deltion
    // TODO: are constant AND that it matches the readIndex just before (forward) or
    // TODO: after (reverse) the deletion.
    //
    // TODO: Test an insertion followed immediately by a deletion and the opposite?
    // TODO: Test an insertion and deletion at the head or tail? Although I suspect GSNAP would
    // TODO: just pad these out, might be the case for a synthetically created Alignment+Seqvars.
    //
    // TODO: Convert these to a purely synthetically created Alignment to make it easier to add
    // TODO: other tests.

    @Test
    public void testIterateSorted() throws IOException {
        MyIterateSortedAlignments alignmentIterator = new MyIterateSortedAlignments();
        final String[] singleBasename = {"test-data/seq-var-test/sorted-seq-var-reads-gsnap"};
        alignmentIterator.iterate(singleBasename);
        alignmentIterator.removeWithoutSeqvars();
        Int2ObjectMap<PerQueryAlignmentData> alignmentDataMap = alignmentIterator.queryIndexToAlignmentDataMap;
        int[] alignmentQueryIndexes = alignmentDataMap.keySet().toIntArray();
        Arrays.sort(alignmentQueryIndexes);
        for (int queryIndex : alignmentQueryIndexes) {
            PerQueryAlignmentData align = alignmentDataMap.get(queryIndex);

            if (LOG.isDebugEnabled()) {
                LOG.debug(String.format("%n------------------------------------%nqueryIndex=%d, data=%s",
                        queryIndex, align.toString()));
                LOG.debug(String.format("leftPadding(%d)=queryPosition(%d)",
                        align.leftPadding, align.queryPosition));
                LOG.debug(String.format("rightPadding(%d)=(queryLength(%d) + numDeletions(%d)) - " +
                        "(targetAlignedLength(%d) + numInsertions(%d)) - leftPadding(%d)",
                        align.rightPadding, align.queryLength, align.numDeletions, 
                        align.targetAlignedLength, align.numInsertions, align.leftPadding));
            }

            assertEquals(String.format("queryIndex=%d readLength(%d) should equal queryLength(%d) + numDeletions(%d)",
                    queryIndex,
                    align.readLength, align.queryLength, align.numDeletions),
                    align.queryLength + align.numDeletions, align.readLength);
            assertEquals(String.format("queryIndex=%d readLength(%d) should equal refLength(%d)", queryIndex,
                    align.readLength, align.refLength),
                    align.readLength, align.refLength);
            assertEquals(String.format("queryIndex=%d  minRefPosition(%d) should equal" +
                    "targetPosition(%d) - queryPosition(%d) + 1", queryIndex,
                    align.minRefPosition, align.targetPosition, align.queryPosition),
                    align.targetPosition - align.queryPosition + 1, align.minRefPosition);
            assertEquals(String.format("queryIndex=%d maxRefPosition(%d) should be equal to " +
                    "targetPosition(%d) + targetAlignedLength(%d) + rightPadding(%d)",
                    queryIndex,
                    align.maxRefPosition, align.targetPosition, align.targetAlignedLength, align.rightPadding),
                    align.targetPosition + align.targetAlignedLength + align.rightPadding, align.maxRefPosition);
            assertEquals(String.format("queryIndex=%d minReadIndex(%d) should be equal to 1",
                    queryIndex,
                    align.minReadIndex),
                    align.minReadIndex, 1);
            assertEquals(String.format("queryIndex=%d maxReadIndex(%d) should be equal to queryLength(%d)",
                    queryIndex,
                    align.maxReadIndex, align.queryLength),
                    align.maxReadIndex, align.queryLength);
            if (align.reverseStrand) {
                assertEquals(String.format("queryIndex=%d (reverse strand) firstReadIndex(%d) should be equal to " +
                        "queryLength(%d)", queryIndex,
                        align.firstReadIndex, align.queryLength),
                        align.firstReadIndex, align.queryLength);
                assertEquals(String.format("queryIndex=%d (reverse strand) firstReadIndex(%d) should be equal to " +
                        "maxReadIndex(%d)", queryIndex,
                        align.firstReadIndex, align.maxReadIndex),
                        align.firstReadIndex, align.maxReadIndex);
            } else {
                assertEquals(String.format("queryIndex=%d (forward strand) firstReadIndex(%d) should be equal to 1",
                        queryIndex, align.firstReadIndex),
                        align.firstReadIndex, 1);
                assertEquals(String.format("queryIndex=%d (forward strand) firstReadIndex(%d) should be equal to" +
                        "minReadIndex(%d)",
                        queryIndex,
                        align.firstReadIndex, align.minReadIndex),
                        align.firstReadIndex, align.minReadIndex);
            }
        }

        // Read the per-base sequence variations from here
        Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap = readSeqVarFile(
                "test-data/seq-var-test/seq-var-reads-gsnap.display-seq-var-tsv-base.tsv");
        // ... Update the data from the sequence variations file to add the queries from the compact-reads file

        int[] seqvarQueryIndexes = seqvarDataMap.keySet().toIntArray();
        Arrays.sort(seqvarQueryIndexes);
        assertArrayEquals("alignmentDataMap indexes doesn't match seqvarDataMap indexes",
                seqvarQueryIndexes, alignmentQueryIndexes);
        for (int queryIndex : seqvarQueryIndexes) {
            PerQueryAlignmentData align = alignmentDataMap.get(queryIndex);
            PerQueryAlignmentData var = seqvarDataMap.get(queryIndex);
            if (LOG.isDebugEnabled()) {
                LOG.debug(String.format("%n------------------------------------%nqueryIndex=%d, data=%s",
                        queryIndex, var.toString()));
            }
            Map<String, String> alignSeqVarsMap = align.refPositionReadIndexToBaseMap;
            Map<String, String> varSeqVarsMap = var.refPositionReadIndexToBaseMap;
            assertEquals(String.format("queryIndex=%d alignSeqVarsMap.size()(%d) should equal varSeqVarsMap.size()(%d)",
                    queryIndex,
                    alignSeqVarsMap.size(), varSeqVarsMap.size()),
                    alignSeqVarsMap.size(), varSeqVarsMap.size());
            for (Map.Entry<String, String> varEntry : varSeqVarsMap.entrySet()) {
                // Make sure the sequence variations match
                final String varEntryBases = varEntry.getValue();
                final String alignEntryBases = alignSeqVarsMap.get(varEntry.getKey());
                assertNotNull(String.format("queryIndex=%d Could not find alignSeqVarsMap entry for %s",
                        queryIndex, varEntry.getKey()),
                        alignEntryBases);
                assertEquals(String.format("queryIndex=%d alignEntryBases(%s) should equal varEntryBases(%s)",
                        queryIndex, alignEntryBases, varEntryBases),
                        alignEntryBases, varEntryBases);
            }
        }
    }

    private Int2ObjectMap<PerQueryAlignmentData> readSeqVarFile(final String filename) throws IOException {
        Int2ObjectMap<PerQueryAlignmentData> seqvarDataMap = new Int2ObjectOpenHashMap<PerQueryAlignmentData>();
        final File datafile = new File(filename);
        final TsvToFromMap tsvReader = TsvToFromMap.createFromTsvFile(datafile);
        for (final LinkedHashToMultiTypeMap<String> dataline : new TsvLineIterator(datafile, tsvReader)) {
            int queryIndex = dataline.getInt("query-index");
            int refPosition = dataline.getInt("position-on-reference");
            int readIndex = dataline.getInt("read-index");
            char fromBase = dataline.getString("var-from").charAt(0);
            char toBase = dataline.getString("var-to").charAt(0);

            PerQueryAlignmentData var = seqvarDataMap.get(queryIndex);
            if (var == null) {
                var = new PerQueryAlignmentData();
                seqvarDataMap.put(queryIndex, var);
            }
            var.observe(refPosition, readIndex, fromBase, toBase);
        }
        return seqvarDataMap;
    }

    private class PerQueryAlignmentData {
        String query;           // From compact-reads file
        int queryLength;        // From compact-reads file, verified against alignment
        int minRefPosition;     // From seq-var
        int maxRefPosition;     // From seq-var
        int minReadIndex;       // From seq-var
        int maxReadIndex;       // From seq-var
        int numDeletions;       // From seq-var, from display-seq-var TSV file
        int numInsertions;      // From seq-var, from display-seq-var TSV file
        int numMismatches;      // From seq-var, from display-seq-var TSV file
        int refLength;          // Calculated, should equal queryLength
        int readLength;         // Calculated, should equal queryLength
        int targetPosition;     // From alignment
        int queryPosition;      // From alignment
        int targetAlignedLength; // From alignment
        int queryAlignedLength;   // From alignment
        boolean reverseStrand;  // From alignment
        int firstReadIndex;     // From alignment
        int leftPadding;        // Calculated
        int rightPadding;       // Calculated
        Object2ObjectMap<String, String> refPositionReadIndexToBaseMap;
        public PerQueryAlignmentData() {
            query = "";
            queryLength = 0;
            minRefPosition = Integer.MAX_VALUE;
            maxRefPosition = Integer.MIN_VALUE;
            minReadIndex = Integer.MAX_VALUE;
            maxReadIndex = Integer.MIN_VALUE;
            numDeletions = 0;
            numInsertions = 0;
            numMismatches = 0;
            refLength = 0;
            readLength = 0;
            targetPosition = 0;
            queryPosition = 0;
            targetAlignedLength = 0;
            queryAlignedLength = 0;
            reverseStrand = false;
            firstReadIndex = -1;
            queryPosition = -1;
            rightPadding = -1;
            refPositionReadIndexToBaseMap = new Object2ObjectOpenHashMap<String, String>();
        }
        public void observe(int refPosition, int readIndex) {
            minRefPosition = Math.min(refPosition, minRefPosition);
            maxRefPosition = Math.max(refPosition, maxRefPosition);
            minReadIndex = Math.min(readIndex, minReadIndex);
            maxReadIndex = Math.max(readIndex, maxReadIndex);
            refLength = maxRefPosition - minRefPosition + numInsertions + 1;
            readLength = maxReadIndex - minReadIndex + numDeletions + 1;
        }
        public void observe(int refPosition, int readIndex, char fromBase, char toBase) {
            if (toBase == '-') {
                numDeletions++;
            } if (fromBase == '-') {
                numInsertions++;
            } else {
                numMismatches++;
            }
            observe(refPosition, readIndex);
            refPositionReadIndexToBaseMap.put(
                    String.format("%d:%d", refPosition, readIndex),
                    String.format("%c->%c", fromBase, toBase));
            leftPadding = queryPosition;
            rightPadding = (queryLength + numDeletions) - (targetAlignedLength + numInsertions) - leftPadding;
        }
        public void setQuery(String query) {
            this.query = query;
            this.queryLength = query.length();
        }

        public String toString() {
            StringBuffer result = new StringBuffer();
            result.append(String.format(
                "{%n  query=%s%n  queryLength=%d%n  minRefPosition=%d%n  maxRefPosition=%d%n  minReadIndex=%d%n" +
                "  maxReadIndex=%d%n  numDeletions=%d%n  numInsertions=%d%n  numMismatches=%d%n  refLength=%d%n" +
                "  readLength=%d%n  targetPosition=%d%n  queryPosition=%d%n  firstReadIndex=%d%n" +
                "  targetAlignedLength=%d%n  queryAlignedLength=%d%n  reverseStrand=%s%n" +
                "  leftPadding=%d%n  rightPadding=%d%n",
                    query, queryLength, minRefPosition, maxRefPosition, minReadIndex,
                    maxReadIndex, numDeletions, numInsertions, numMismatches, refLength,
                    readLength, targetPosition, queryPosition, firstReadIndex,
                    targetAlignedLength, queryAlignedLength,
                    reverseStrand ? "true" : "false", leftPadding, rightPadding));
            result.append("  seqvars=[");
            for (Map.Entry<String,String> entry : refPositionReadIndexToBaseMap.entrySet()) {
                result.append(String.format("%s %s, ", entry.getKey(), entry.getValue()));
            }
            result.append(String.format("]}"));
            return result.toString();
        }
    }

    private class CountsAtPosition {
    }

    private class MyIterateSortedAlignments extends IterateSortedAlignments<CountsAtPosition> {

        Int2ObjectMap<PerQueryAlignmentData> queryIndexToAlignmentDataMap;

        public MyIterateSortedAlignments() {
            queryIndexToAlignmentDataMap = new Int2ObjectOpenHashMap<PerQueryAlignmentData>();
        }

        @Override
        public void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry,
                                         Int2ObjectMap<CountsAtPosition> positionToBases,
                                         int currentReferenceIndex, int currentRefPosition, int currentReadIndex) {
            if (LOG.isDebugEnabled()) {
                LOG.debug(String.format("RB: queryIndex=%d\tref_position=%d\tread_index=%d",
                    alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex));
            }
            int queryIndex = alignmentEntry.getQueryIndex();
            if (currentReadIndex >= 1) {
                PerQueryAlignmentData alignmentData = queryIndexToAlignmentDataMap.get(queryIndex);
                if (alignmentData == null) {
                    alignmentData = new PerQueryAlignmentData();
                    queryIndexToAlignmentDataMap.put(queryIndex, alignmentData);
                }
                if (alignmentData.firstReadIndex == -1) {
                    alignmentData.firstReadIndex = currentReadIndex;
                    alignmentData.queryPosition = alignmentEntry.getQueryPosition();
                    alignmentData.targetPosition = alignmentEntry.getPosition();
                    alignmentData.queryLength = alignmentEntry.getQueryLength();
                    alignmentData.queryAlignedLength = alignmentEntry.getQueryAlignedLength();
                    alignmentData.targetAlignedLength = alignmentEntry.getTargetAlignedLength();
                    alignmentData.reverseStrand = alignmentEntry.getMatchingReverseStrand();
                }
                alignmentData.observe(currentRefPosition, currentReadIndex);
            } else {
                throw new RuntimeException(String.format("readIndex=%d should be >=1, queryIndex=",
                        currentReadIndex, queryIndex));
            }
        }

        @Override
        public void observeVariantBase(
                ConcatSortedAlignmentReader sortedReaders, Alignments.AlignmentEntry alignmentEntry,
                Int2ObjectMap<CountsAtPosition> positionToBases, Alignments.SequenceVariation var,
                char toChar, char fromChar, byte toQual, int currentReferenceIndex,
                int currentRefPosition, int currentReadIndex) {
            if (LOG.isDebugEnabled()) {
                LOG.debug(String.format("VB: queryIndex=%d\tref_position=%d\tread_index=%d\tfromChar=%c\ttoChar=%c",
                        alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex, fromChar, toChar));
            }

            int queryIndex = alignmentEntry.getQueryIndex();
            if (currentReadIndex >= 1) {
                PerQueryAlignmentData alignmentData = queryIndexToAlignmentDataMap.get(queryIndex);
                if (alignmentData == null) {
                    alignmentData = new PerQueryAlignmentData();
                    queryIndexToAlignmentDataMap.put(queryIndex, alignmentData);
                }
                if (alignmentData.firstReadIndex == -1) {
                    alignmentData.firstReadIndex = currentReadIndex;
                    alignmentData.queryPosition = alignmentEntry.getQueryPosition();
                    alignmentData.targetPosition = alignmentEntry.getPosition();
                    alignmentData.queryLength = alignmentEntry.getQueryLength();
                    alignmentData.queryAlignedLength = alignmentEntry.getQueryAlignedLength();
                    alignmentData.targetAlignedLength = alignmentEntry.getTargetAlignedLength();
                    alignmentData.reverseStrand = alignmentEntry.getMatchingReverseStrand();
                }
                alignmentData.observe(currentRefPosition, currentReadIndex, fromChar, toChar);
            } else {
                throw new RuntimeException(String.format("readIndex=%d should be >=1, queryIndex=",
                        currentReadIndex, queryIndex));
            }
        }

        @Override
        public void processPositions(int referenceIndex, int intermediatePosition, CountsAtPosition positionBaseInfos) {
        }

        /**
         * Remove any items from the map that don't have sequence variations.
         */
        public void removeWithoutSeqvars() {
            IntSet toRemoveQueryIndexes = new IntArraySet(queryIndexToAlignmentDataMap.size());
            for (Map.Entry<Integer, PerQueryAlignmentData> entry : queryIndexToAlignmentDataMap.entrySet()) {
                if (entry.getValue().refPositionReadIndexToBaseMap.size() == 0) {
                    toRemoveQueryIndexes.add(entry.getKey());
                }
            }
            for (int queryIndex : toRemoveQueryIndexes) {
                queryIndexToAlignmentDataMap.remove(queryIndex);
            }
        }
    }
}
