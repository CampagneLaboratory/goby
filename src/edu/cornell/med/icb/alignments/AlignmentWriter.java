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

package edu.cornell.med.icb.alignments;

import edu.cornell.med.icb.goby.reads.MessageChunksWriter;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import java.io.Closeable;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

/**
 * This class write alignments in a Protocol Buffer, compressed, binary and splitable format.
 * See Alignements.proto for the specification of this format.
 *
 * @author Fabien Campagne
 *         Date: Apr 30, 2009
 *         Time: 5:53:50 PM
 */
public class AlignmentWriter implements Closeable {
    private Alignments.AlignmentCollection.Builder collectionBuilder;
    private MessageChunksWriter entriesChunkWriter;
    private IndexedIdentifier queryIdentifiers;
    private IndexedIdentifier targetIdentifiers;
    private boolean headerWritten;
    private FileOutputStream headerOutput;
    private int[] queryLengths;

    private String[] queryIdentifiersArray;
    private String[] targetIdentifiersArray;
    private int maxTargetIndex = -1;
    private int maxQueryIndex = -1;
    private Properties stats;
    private boolean statsWritten;
    private FileWriter statsWriter;
    private String basename;
    private int numberOfAlignedReads = 0;

    public AlignmentWriter(final String outputBasename) throws IOException {
        final FileOutputStream alignmentEntries = new FileOutputStream(outputBasename + ".entries");
        headerOutput = new FileOutputStream(outputBasename + ".header");
        statsWriter = new FileWriter(outputBasename + ".stats");
        this.basename = outputBasename;
        collectionBuilder = Alignments.AlignmentCollection.newBuilder();
        entriesChunkWriter = new MessageChunksWriter(alignmentEntries);
        newEntry = Alignments.AlignmentEntry.newBuilder();
        queryIdentifiers = new IndexedIdentifier();
        targetIdentifiers = new IndexedIdentifier();
        stats = new Properties();
        // we assume stats were written until a client puts stats in this writer.
        statsWritten = true;

    }

    private Alignments.AlignmentEntry.Builder newEntry;

    public final void setQueryIndex(final int queryIndex) {
        newEntry.setQueryIndex(queryIndex);
    }

    public final void setTargetIndex(final int referenceIndex) {
        newEntry.setTargetIndex(referenceIndex);
    }

    public final void setTargetPosition(final int position) {
        newEntry.setPosition(position);
    }

    public final void setAlignmentScore(final float score) {
        newEntry.setScore(score);
    }

    public void setNumAlignmentEntriesPerChunk(final int numEntriesPerChunk) {
        entriesChunkWriter.setNumEntriesPerChunk(numEntriesPerChunk);
    }

    public final void setAlignmentEntry(final int queryIndex, final int referenceIndex,
                                        final int position,
                                        final float score, final boolean matchesReverseStrand) {
        newEntry.setQueryIndex(queryIndex);
        newEntry.setTargetIndex(referenceIndex);
        newEntry.setScore(score);
        newEntry.setPosition(position);
        newEntry.setMatchingReverseStrand(matchesReverseStrand);
    }

    /**
     * Obtain the alignment entry that is being prepared. Set values on the entry, then call appendAlignmentEntry()
     *
     * @return the current alignment entry.
     */
    public Alignments.AlignmentEntry.Builder getAlignmentEntry() {
        return newEntry;
    }

    /**
     * Append the current entry to the file being written.
     *
     * @throws IOException
     */
    public synchronized void appendEntry() throws IOException {

        final Alignments.AlignmentEntry builtEntry = newEntry.build();
        maxQueryIndex = Math.max(builtEntry.getQueryIndex(), maxQueryIndex);
        maxTargetIndex = Math.max(builtEntry.getTargetIndex(), maxTargetIndex);
        this.collectionBuilder.addAlignmentEntries(builtEntry);
        entriesChunkWriter.writeAsNeeded(collectionBuilder, builtEntry.getMultiplicity());
        newEntry = Alignments.AlignmentEntry.newBuilder();
    }

    /**
     * Append an entry to the file being written. The entry must have been prepared outside this writter.
     *
     * @param builtEntry
     * @throws IOException If an error occurs writting this entry.
     */
    public synchronized void appendEntry(final Alignments.AlignmentEntry builtEntry) throws IOException {
        this.collectionBuilder.addAlignmentEntries(builtEntry);
        entriesChunkWriter.writeAsNeeded(collectionBuilder, builtEntry.getMultiplicity());

        maxQueryIndex = Math.max(builtEntry.getQueryIndex(), maxQueryIndex);
        maxTargetIndex = Math.max(builtEntry.getTargetIndex(), maxTargetIndex);
        numberOfAlignedReads += 1;
    }

    /**
     * Replace the target entry, then append the entry.
     *
     * @param entry
     * @param replacementTargetIndex the replacement targetIndex.
     * @throws IOException
     */
    public void appendEntry(Alignments.AlignmentEntry entry, final int replacementTargetIndex) throws IOException {
        entry = entry.newBuilderForType().mergeFrom(entry).setTargetIndex(replacementTargetIndex).build();
        this.appendEntry(entry);
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws IOException {
        writeHeader();
        writeStats();
        IOUtils.closeQuietly(headerOutput);
        entriesChunkWriter.close(collectionBuilder);
    }

    private void writeHeader() throws IOException {
        if (!headerWritten) {
            final Alignments.AlignmentHeader.Builder headerBuilder = Alignments.AlignmentHeader.newBuilder();

            headerBuilder.setNumberOfTargets(maxTargetIndex + 1);
            headerBuilder.setNumberOfQueries(maxQueryIndex + 1);

            headerBuilder.setQueryNameMapping(getMapping(queryIdentifiers, queryIdentifiersArray));
            headerBuilder.setTargetNameMapping(getMapping(targetIdentifiers, targetIdentifiersArray));
            headerBuilder.setNumberOfAlignedReads(numberOfAlignedReads);
            if (queryLengths != null) {
                headerBuilder.addAllQueryLength(IntArrayList.wrap(queryLengths));
            }

            headerBuilder.build().writeTo(headerOutput);
            headerWritten = true;
        }
    }

    private void writeStats() throws IOException {
        if (!statsWritten) {
            stats.put("basename", FilenameUtils.getBaseName(basename));
            stats.put("basename.full", basename);
            stats.put("number.aligned.reads", Integer.toString(numberOfAlignedReads));
            stats.store(statsWriter, "Statistics for merged alignment. ");

            statsWritten = true;
        }
    }

    /**
     * Provide query identifiers as an array of string, where queryIndex is the index of the element in the array.
     *
     * @param queryIdentifiersArray
     */
    public void setQueryIdentifiersArray(final String[] queryIdentifiersArray) {
        this.queryIdentifiersArray = queryIdentifiersArray;
        maxQueryIndex = queryIdentifiersArray.length - 1;
    }

    /**
     * Provide target identifiers as an array of string, where targetIndex is the index of the element in the array.
     *
     * @param targetIdentifiersArray
     */
    public void setTargetIdentifiersArray(final String[] targetIdentifiersArray) {
        this.targetIdentifiersArray = targetIdentifiersArray;
    }

    private Alignments.IdentifierMapping getMapping(final IndexedIdentifier identifiers, final String[] targetIdentifiersArray) {
        final Alignments.IdentifierMapping.Builder result = Alignments.IdentifierMapping.newBuilder();

        final ObjectList<Alignments.IdentifierInfo> mappings =
                new ObjectArrayList<Alignments.IdentifierInfo>();

        if (targetIdentifiersArray == null) {
            for (final MutableString id : identifiers.keySet()) {

                mappings.add(Alignments.IdentifierInfo.newBuilder()
                        .setName(id.toString())
                        .setIndex(identifiers.get(id))
                        .build());

            }
        } else {
            // prefer to use the array when available, for speed:

            for (int index = 0; index < targetIdentifiersArray.length; index++) {
                mappings.add(Alignments.IdentifierInfo.newBuilder()
                        .setName(targetIdentifiersArray[index])
                        .setIndex(index)
                        .build());

            }

        }

        result.addAllMappings(mappings);
        return result.build();
    }

    public void printStats(final PrintStream out) {
        this.entriesChunkWriter.printStats(out);
        out.println("Number of queries: " + (maxQueryIndex + 1));
        out.println("Number of targets: " + (maxTargetIndex + 1));
    }


    public void setQueryIdentifiers(final IndexedIdentifier queryIdentifiers) {
        this.queryIdentifiers = queryIdentifiers;
        for (final int index : queryIdentifiers.values()) {
            maxQueryIndex = Math.max(maxQueryIndex, index);
        }

    }

    public void setTargetIdentifiers(final IndexedIdentifier targetIdentifiers) {
        this.targetIdentifiers = targetIdentifiers;


    }

    public void setQueryLengths(final int[] queryLengths) {
        assert queryLengths.length > maxQueryIndex :
                "The number of elements of queryLength is too small to accomodate queryIndex=" + maxQueryIndex;
        this.queryLengths = queryLengths;

    }

    /**
     * Set the total number of queries.
     *
     * @param numQueries The number of query sequences.
     */
    public void setNumQueries(final int numQueries) {
        maxQueryIndex = numQueries - 1;
    }

    /**
     * Set the total number of targets.
     *
     * @param numTargets The number of target sequences.
     */
    public void setNumTargets(final int numTargets) {
        maxTargetIndex = numTargets - 1;
    }

    public void putStatistic(String description, String value) {
        statsWritten = false;
        stats.put(description, value);
    }

    public void putStatistic(String description, double value) {
        putStatistic(description, String.format("%3.3g", value));
    }

    public void putStatistic(String description, int value) {
        putStatistic(description, String.format("%d", value));
    }

    public void setStatistics(Properties statistics) {
        for (Object key : statistics.keySet()) {
            putStatistic(key.toString(), statistics.get(key).toString());
        }
    }
}
