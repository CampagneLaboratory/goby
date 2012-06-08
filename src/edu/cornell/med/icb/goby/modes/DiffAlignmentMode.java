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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.chars.CharArrayList;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.util.Collections;
import java.util.Map;

/**
 * Report differences between two alignment files.  Note that this mode cannot handle pair alignments because it does
 * not take fragment indices into account when comparing entries. On such alignments, comparisons will be unreliable.
 *
 * @author Fabien Campagne
 *         Date: Apr 28, 2009
 *         Time: 6:03:56 PM
 */
public class DiffAlignmentMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "diff-alignments";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Report differences between two alignment files. Note that this mode cannot handle paired alignments " +
                    "because it does not take fragment indices into account when comparing entries. " +
                    "On such alignments, comparisons will be unreliable.";

    private String[] inputFilenames;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        inputFilenames = jsapResult.getStringArray("input");
        return this;
    }

    /**
     * Perform the concatenation.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        final String[] basenames = AlignmentReaderImpl.getBasenames(inputFilenames);
        if (basenames.length != 2) {
            System.out.println("You can only compare two alignments.");
            System.exit(1);
        }
        System.out.println("First alignment: " + basenames[0]);
        System.out.println("Second alignment: " + basenames[1]);

        final ProgressLogger progress = new ProgressLogger();
        progress.start("comparing alignment entries");

        final AlignmentReaderImpl readerA = new AlignmentReaderImpl(basenames[0]);
        final AlignmentReaderImpl readerB = new AlignmentReaderImpl(basenames[1]);
        readerA.readHeader();
        readerB.readHeader();
        final AlignmentReaderImpl[] readers = new AlignmentReaderImpl[2];
        readers[0] = readerA;
        readers[1] = readerB;

        final IntSet[] queryIndices = new IntSet[2];
        queryIndices[0] = new IntOpenHashSet();
        queryIndices[1] = new IntOpenHashSet();

        final IntArrayList[] queryMultiplicity = new IntArrayList[2];
        queryMultiplicity[0] = new IntArrayList();
        queryMultiplicity[1] = new IntArrayList();

        final IntArrayList[] targetIndices = new IntArrayList[2];
        targetIndices[0] = new IntArrayList();
        targetIndices[1] = new IntArrayList();

        final FloatArrayList[] scores = new FloatArrayList[2];
        scores[0] = new FloatArrayList();
        scores[1] = new FloatArrayList();

        final IntArrayList[] positions = new IntArrayList[2];
        positions[0] = new IntArrayList();
        positions[1] = new IntArrayList();

        final CharArrayList[] strands = new CharArrayList[2];
        strands[0] = new CharArrayList();
        strands[1] = new CharArrayList();

        int readerIndex = 0;

        // Map the target indexes of the SECOND alignment to the first
        Object2IntMap<MutableString> masterTargetNameToIndexMap = new Object2IntOpenHashMap<MutableString>();
        Int2ObjectOpenHashMap<MutableString> masterTargetIndexToNameMap = new Int2ObjectOpenHashMap<MutableString>();
        Int2IntMap[] targetsIndexesToMasterMap = new Int2IntMap[2];
        targetsIndexesToMasterMap[0] = addIdentiersToMaster(
                masterTargetNameToIndexMap, masterTargetIndexToNameMap, readerA.getTargetIdentifiers());
        targetsIndexesToMasterMap[1] = addIdentiersToMaster(
                masterTargetNameToIndexMap, masterTargetIndexToNameMap, readerB.getTargetIdentifiers());

        for (final AlignmentReaderImpl reader : readers) {
            queryMultiplicity[readerIndex].size(reader.getNumberOfQueries());
            targetIndices[readerIndex].size(reader.getNumberOfQueries());
            positions[readerIndex].size(reader.getNumberOfQueries());
            scores[readerIndex].size(reader.getNumberOfQueries());
            strands[readerIndex].size(reader.getNumberOfQueries());

            for (final Alignments.AlignmentEntry entry : reader) {
                final int queryIndex = entry.getQueryIndex();
                queryIndices[readerIndex].add(queryIndex);

                queryMultiplicity[readerIndex].set(queryIndex, entry.getMultiplicity());
                targetIndices[readerIndex].set(
                        queryIndex, targetsIndexesToMasterMap[readerIndex].get(entry.getTargetIndex()));
                scores[readerIndex].set(queryIndex, entry.getScore());
                positions[readerIndex].set(queryIndex, entry.getPosition());
                strands[readerIndex].set(queryIndex, entry.getMatchingReverseStrand() ? '-' : '+');
            }
            readerIndex++;
        }
        if (readerA.getNumberOfQueries() == 0 && readerB.getNumberOfQueries() == 0) {
            System.out.println("Both alignments have no entries. They are equal.");
            System.exit(0);
        }
        final IntSet commonQueryIndices = new IntOpenHashSet();
        commonQueryIndices.addAll(queryIndices[0]);
        commonQueryIndices.retainAll(queryIndices[1]);
        double percentCommon = commonQueryIndices.size();
        percentCommon /= (double) Math.max(queryIndices[0].size(), queryIndices[1].size());
        percentCommon *= 100;
        System.out.printf("%3.5g %% of query indices match between the two alignments.%n", percentCommon);
        // check multiplicity and target Index for common indices:

        for (final int commonQueryIndex : commonQueryIndices) {
            final int firstMultiplicity = queryMultiplicity[0].getInt(commonQueryIndex);
            final int secondMultiplicity = queryMultiplicity[1].getInt(commonQueryIndex);
            if (firstMultiplicity != secondMultiplicity) {
                System.out.printf("Entry queryIndex=%d has different multiplicity between "
                        + "the two alignments (first=%d, second=%d).%n", commonQueryIndex,
                        firstMultiplicity, secondMultiplicity);
            }
            final int firstTargetIndex = targetIndices[0].getInt(commonQueryIndex);
            final int secondTargetIndex = targetIndices[1].getInt(commonQueryIndex);
            if (firstTargetIndex != secondTargetIndex) {
                System.out.printf("Entry queryIndex=%d has different targetIndex between "
                        + "the two alignments (first=%s, second=%s).%n", commonQueryIndex,
                        masterTargetIndexToNameMap.get(firstTargetIndex),
                        masterTargetIndexToNameMap.get(secondTargetIndex));
            }

            final float firstScore = scores[0].getFloat(commonQueryIndex);
            final float secondScore = scores[1].getFloat(commonQueryIndex);
            if (firstScore != secondScore) {
                System.out.printf("Entry queryIndex=%d has different score between "
                        + "the two alignments (first=%g, second=%g).%n", commonQueryIndex,
                        firstScore, secondScore);
            }
        }

        if (commonQueryIndices.size() != Math.max(queryIndices[0].size(), queryIndices[1].size())) {
            System.out.println("--- Unique to A ----------");
            final IntSet setA = queryIndices[0];
            final IntList uniqueToFirst = uniqueTo(commonQueryIndices, setA, "first");
            printEntries(uniqueToFirst, targetIndices[0], queryMultiplicity[0], scores[0],
                    positions[0], strands[0], masterTargetIndexToNameMap);
            System.out.println("--- Unique to B ----------");
            final IntSet setB = queryIndices[1];
            final IntList uniqueToSecond = uniqueTo(commonQueryIndices, setB, "second");
            printEntries(uniqueToSecond, targetIndices[1], queryMultiplicity[1], scores[1],
                    positions[1], strands[1], masterTargetIndexToNameMap);
            System.out.println("--- Common to A and B ----");
            printCommonEntries(commonQueryIndices,  targetIndices, queryMultiplicity, scores, positions, strands,
                    masterTargetIndexToNameMap);

        }
        progress.stop();
    }

    private Int2IntMap addIdentiersToMaster(
            final Object2IntMap<MutableString> masterTargetNameToIndexMap,
            final Int2ObjectMap<MutableString> masterTargetIndexToNameMap,
            final IndexedIdentifier targetIdents) {
        Int2IntMap targetsToMasterMap = new Int2IntOpenHashMap();
        int nextIndex = masterTargetNameToIndexMap.size();
        for (Map.Entry<MutableString, Integer> entry : targetIdents.entrySet()) {
            Integer existingIndex = masterTargetNameToIndexMap.get(entry.getKey());
            if (existingIndex == null) {
                masterTargetNameToIndexMap.put(entry.getKey(), nextIndex);
                masterTargetIndexToNameMap.put(nextIndex, entry.getKey());
                targetsToMasterMap.put((int) entry.getValue(), nextIndex);
                nextIndex++;
            } else {
                targetsToMasterMap.put((int) entry.getValue(), (int) existingIndex);
            }
        }
        return targetsToMasterMap;
    }

    private void printEntries(
            final IntList queryIndices, final IntArrayList targetIndices, final IntArrayList multiplicity,
            final FloatArrayList scores, final IntArrayList positions, final CharArrayList strands,
            final Int2ObjectMap<MutableString> masterTargetIndexToNameMap) {
        for (final int queryIndex : queryIndices) {
            System.out.printf("[entry queryIndex=%d target=%s multiplicity=%d score=%g position=%d strand=%c]%n",
                    queryIndex, masterTargetIndexToNameMap.get(targetIndices.get(queryIndex)),
                    multiplicity.get(queryIndex), scores.get(queryIndex),
                    positions.get(queryIndex), strands.get(queryIndex));
        }
    }

    private void printCommonEntries(
            final IntSet queryIndicesSet,
            final IntArrayList[] targetIndices, final IntArrayList[] multiplicity,
            final FloatArrayList[] scores, final IntArrayList[] positions, final CharArrayList[] strands,
            final Int2ObjectMap<MutableString> masterTargetIndexToNameMap) {
        IntList queryIndices = new IntArrayList(queryIndicesSet.size());
        queryIndices.addAll(queryIndicesSet);
        Collections.sort(queryIndices);
        long exactMatches = 0;
        for (final int queryIndex : queryIndices) {
            // Check for exact match, disregarding score
            if ((targetIndices[0].get(queryIndex).equals(targetIndices[1].get(queryIndex))) &&
                    (multiplicity[0].get(queryIndex).equals(multiplicity[1].get(queryIndex))) &&
                    (positions[0].get(queryIndex).equals(positions[1].get(queryIndex))) &&
                    (strands[0].get(queryIndex) == strands[1].get(queryIndex))) {
                exactMatches++;
            } else {
                System.out.printf("[entry queryIndex=%d target=%s multiplicity=%d score=%g position=%d strand=%c%n",
                        queryIndex, masterTargetIndexToNameMap.get(targetIndices[0].get(queryIndex)),
                        multiplicity[0].get(queryIndex), scores[0].get(queryIndex),
                        positions[0].get(queryIndex), strands[0].get(queryIndex));
                System.out.printf(" entry queryIndex=%d target=%s multiplicity=%d score=%g position=%d strand=%c]%n",
                        queryIndex, masterTargetIndexToNameMap.get(targetIndices[1].get(queryIndex)),
                        multiplicity[1].get(queryIndex), scores[1].get(queryIndex),
                        positions[1].get(queryIndex), strands[1].get(queryIndex));
            }
        }
        System.out.printf("Disregarding score, there were %d exact matches and %d inexact matches.%n",
                exactMatches, queryIndices.size() - exactMatches);
    }

    private IntList uniqueTo(final IntSet commonQueryIndices, final IntSet a, final String label) {
        final IntSet queryIndicesUniqueToA = new IntOpenHashSet();
        queryIndicesUniqueToA.addAll(a);
        queryIndicesUniqueToA.removeAll(commonQueryIndices);
        final IntList sorted = new IntArrayList();
        sorted.addAll(queryIndicesUniqueToA);
        Collections.sort(sorted);
        System.out.printf("The %s alignment contained %d unique query indicies%n", label, sorted.size());
        return sorted;
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new DiffAlignmentMode().configure(args).execute();
    }
}
