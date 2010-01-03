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

package edu.cornell.med.icb.goby.maq.filters;

import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.goby.maq.MaqMapEntry;
import edu.cornell.med.icb.goby.maq.MaqMapHeader;
import edu.cornell.med.icb.goby.maq.MaqMapReader;
import edu.cornell.med.icb.tissueinfo.similarity.GeneTranscriptRelationships;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * This class assists with filtering MaqMapEntry values when merging transcript-based
 * alignments.
 * @author Kevin Dorff
 */
public final class TranscriptsFilter extends AbstractMaqMapEntryFilter {
    /** The helps with relationships between genes and transcripts. */
    private final GeneTranscriptRelationships gtr;
    /** The transcript identifiers as represented in gtr. */
    private final IndexedIdentifier transcriptsIndexedIdentifiers;
    /** The map of read-name-index to gene's that read references. */
    private final Int2ObjectMap<IntSet> readIndexToGeneIdSetMap;
    /** The k value for the filter. */
    private final int k;
    /** An array of maq-seq-ref's to transcript-id. */
    private int[] maqMapSeqRefToTranscriptIndex;
    /** The number of times inspectEntry was called. */
    private int numInspected;
    /** How many entries SHOULD be written. */
    private int shouldWrite;

    /**
     * Constructor.
     * @param geneTranscriptFile the gene-transcripts-map file to read
     * @param kVal the k value for the filter
     * @throws FileNotFoundException if the gene-transcripts-map didn't exist
     */
    public TranscriptsFilter(final String geneTranscriptFile, final int kVal)
            throws FileNotFoundException {
        System.out.println("** TRANSCRIPT MODE **");
        this.gtr = new GeneTranscriptRelationships();
        transcriptsIndexedIdentifiers = gtr.load(geneTranscriptFile);
        readIndexToGeneIdSetMap = new Int2ObjectOpenHashMap<IntSet>();
        k = kVal;
        numInspected = 0;
        shouldWrite = -1;
    }

    /**
     * Post processing (after all inspectEntry's have been called).
     */
    @Override
    public void postProcessing() {
    }

    /**
     * Set the new / updated header that is being used when filtering these entries.
     * @param header the Maq Map header
     */
    @Override
    public void setHeader(final MaqMapHeader header) {
        // Use an int[] for this
        maqMapSeqRefToTranscriptIndex = new int[header.getNRef()];
        for (int i = 0; i < header.getNRef(); i++) {
            maqMapSeqRefToTranscriptIndex[i] =
                    transcriptsIndexedIdentifiers.registerIdentifier(
                            header.getRefNameMutable(i));
        }
    }

    /**
     * Inspect an entry (will be called during a first pass of reading the entries).
     * @param entry the entry
     */
    @Override
    public void inspectEntry(final MaqMapEntry entry) {
        numInspected++;
        final int readNameIndex = entry.getReadNameIndex();
        IntSet geneIdSet = readIndexToGeneIdSetMap.get(readNameIndex);
        if (geneIdSet != null && geneIdSet.size() > k) {
       //     System.out.println("read id found to match more than k genes..");
            return;
        }
        final int transcriptIndex = maqMapSeqRefToTranscriptIndex[(int) entry.getReferenceSequenceId()];
        final int geneIndex = gtr.transcript2Gene(transcriptIndex);
        if (geneIdSet == null) {
            // An array set is a fast implementation when the set is expected to be small
            // (a few genes on average per read). The ArraySet is also the smallest footprint implementation,
            // this matters most here where we keep one such set for each read.
            geneIdSet = new IntArraySet();
            readIndexToGeneIdSetMap.put(readNameIndex, geneIdSet);
        }
        if (geneIdSet.size() <= k) {
            // No need to add to the set if the size is already too big (save memory)
            geneIdSet.add(geneIndex);
        }
        if ((numInspected % 100000) == 0) {
            System.out.printf("inspectEntry calls = %,d, readIndex: %d memory used / total = %s%n",
                    numInspected, readNameIndex, getHeapSize());
        }
    }

    /**
     * Get the number of entries that should be written (sanity check).
     * @return the number of entries that should be written
     */
    @Override
    public int getShouldWrite() {
        return shouldWrite;
    }

    /**
     * Determine if this entry should be retained (will be called during a second
     * pass of reading the entries).
     * @param entry the entry
     * @return true if it should be retained
     */
    @Override
    public boolean shouldRetainEntry(final MaqMapEntry entry) {
        return (readIndexToGeneIdSetMap.get(entry.getReadNameIndex()).size() <= k);
    }

    /**
     * I included this to run with JProfiler so I could see
     * where time was spent when reading a large MAP file and filtering
     * the results.
     * @param args command line args
     * @throws IOException error reading
     */
    public static void main(final String[] args) throws IOException {
        final MaqMapReader reader = new MaqMapReader(
                "C:/save/IntelliJ/maq-java/Homo_sapiens.NCBI36.52.cdna.all.00-s_1_sequence.map",
                128, null);
        reader.setReadNameIdentifiers(new IndexedIdentifier());
        reader.setEntryFilter(new TranscriptsFilter(
                "C:/save/IntelliJ/maq-java/Homo_sapiens.NCBI36.52.cdna.all.config", 1));
        for (final MaqMapEntry entry : reader) {
        }
    }
}
