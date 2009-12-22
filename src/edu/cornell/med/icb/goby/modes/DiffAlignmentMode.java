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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.alignments.AlignmentReader;
import edu.cornell.med.icb.alignments.Alignments;
import it.unimi.dsi.fastutil.chars.CharArrayList;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.io.TextIO;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.IOException;
import java.util.Collections;

/**
 * Report differences between two alignment files.
 * <p/>
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
    private static final String MODE_DESCRIPTION = "Report differences between two alignment files.";

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
        final String[] basenames = AlignmentReader.getBasenames(inputFilenames);
        if (basenames.length != 2) {
            System.out.println("You can only compare two alignments.");
            System.exit(1);
        }
        System.out.println("First alignment: " + basenames[0]);
        System.out.println("Second alignment: " + basenames[1]);

        final ProgressLogger progress = new ProgressLogger();
        progress.start("comparing alignment entries");

        final AlignmentReader readerA = new AlignmentReader(basenames[0]);
        final AlignmentReader readerB = new AlignmentReader(basenames[1]);
        readerA.readHeader();
        readerB.readHeader();
        final AlignmentReader[] readers = new AlignmentReader[2];
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

        for (final AlignmentReader reader : readers) {
            queryMultiplicity[readerIndex].size(reader.getNumberOfQueries());
            targetIndices[readerIndex].size(reader.getNumberOfQueries());
            positions[readerIndex].size(reader.getNumberOfQueries());
            scores[readerIndex].size(reader.getNumberOfQueries());
            strands[readerIndex].size(reader.getNumberOfQueries());

            for (final Alignments.AlignmentEntry entry : reader) {
                final int queryIndex = entry.getQueryIndex();
                queryIndices[readerIndex].add(queryIndex);

                queryMultiplicity[readerIndex].set(queryIndex, entry.getMultiplicity());
                targetIndices[readerIndex].set(queryIndex, entry.getTargetIndex());
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
        System.out.printf("%3.5g %% of query indices match between the two alignments.", percentCommon);
        // check multiplicity and target Index for common indices:

        for (final int commonQueryIndex : commonQueryIndices) {
            final int firstMultiplicity = queryMultiplicity[0].getInt(commonQueryIndex);
            final int secondMultiplicity = queryMultiplicity[1].getInt(commonQueryIndex);
            if (firstMultiplicity != secondMultiplicity) {
                System.out.printf("Entry queryIndex=%d has different multiplicity between " +
                        "the two alignments (first=%d, second=%d).  ", commonQueryIndex,
                        firstMultiplicity, secondMultiplicity);
            }
            final int firstTargetIndex = targetIndices[0].getInt(commonQueryIndex);
            final int secondTargetIndex = targetIndices[1].getInt(commonQueryIndex);
            if (firstTargetIndex != secondTargetIndex) {
                System.out.printf("Entry queryIndex=%d has different targetIndex between " +
                        "the two alignments (first=%d, second=%d).  ", commonQueryIndex,
                        firstTargetIndex, secondTargetIndex);
            }

            final float firstScore = scores[0].getFloat(commonQueryIndex);
            final float secondScore = scores[1].getFloat(commonQueryIndex);
            if (firstScore != secondScore) {
                System.out.printf("Entry queryIndex=%d has different score between " +
                        "the two alignments (first=%g, second=%g).  ", commonQueryIndex,
                        firstScore, secondScore);
            }
        }

        if (commonQueryIndices.size() != Math.max(queryIndices[0].size(), queryIndices[1].size())) {
            final IntSet A = queryIndices[0];
            final IntList uniqueToFirst = uniqueTo(commonQueryIndices, A, "first");
            printEntries(uniqueToFirst, targetIndices[0], queryMultiplicity[0], scores[0],
                    positions[0], strands[0]);
            final IntSet B = queryIndices[1];
            final IntList uniqueToSecond = uniqueTo(commonQueryIndices, B, "second");
            printEntries(uniqueToSecond, targetIndices[1], queryMultiplicity[1], scores[1],
                    positions[1], strands[1]);

        }
        progress.stop();
    }

    private void printEntries(final IntList queryIndices, final IntArrayList targetIndices, final IntArrayList multiplicity,
                              final FloatArrayList scores, final IntArrayList positions, final CharArrayList strands) {
        for (final int queryIndex : queryIndices) {
            System.out.printf("[entry queryIndex=%d targetIndex=%d multiplicity=%d score=%g position=%d strand=%c] %n", queryIndex,
                    targetIndices.get(queryIndex),
                    multiplicity.get(queryIndex), scores.get(queryIndex),
                    positions.get(queryIndex), strands.get(queryIndex));
        }
    }

    private IntList uniqueTo(final IntSet commonQueryIndices, final IntSet a, final String label) {
        final IntSet queryIndicesUniqueToA = new IntOpenHashSet();
        queryIndicesUniqueToA.addAll(a);
        queryIndicesUniqueToA.removeAll(commonQueryIndices);
        final IntList sorted = new IntArrayList();
        sorted.addAll(queryIndicesUniqueToA);
        Collections.sort(sorted);
        System.out.println("The following query indices are unique to the " + label + " alignment: ");
        TextIO.storeInts(sorted.toIntArray(), System.out);
        return sorted;
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new DiffAlignmentMode().configure(args).execute();
    }
}
