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

package edu.cornell.med.icb.alignments.filters;

import edu.cornell.med.icb.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.tissueinfo.similarity.GeneTranscriptRelationships;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.log4j.Logger;

import java.io.FileNotFoundException;

/**
 * This class assists with filtering MaqMapEntry values when merging transcript-based
 * alignments.
 *
 * @author Kevin Dorff
 * @author Fabien Campagne
 */
public final class TranscriptsAlignmentFilter extends AbstractAlignmentEntryFilter {
    private static final Logger LOG = Logger.getLogger(TranscriptsAlignmentFilter.class);

    /**
     * The helps with relationships between genes and transcripts.
     */
    private final GeneTranscriptRelationships gtr;
    /**
     * The transcript identifiers as represented in gtr.
     */
    private final IndexedIdentifier transcriptsIndexedIdentifiers;
    /**
     * The map of read-name-index to gene's that read references.
     */
    private final Int2ObjectMap<IntSet> readIndexToGeneIdSetMap;
    /**
     * The k value for the filter.
     */
    private final int k;
    /**
     * An array of reference sequence index to transcript index.
     */
    private int[] alignmentReferenceIndexToTranscriptIndex;
    /**
     * The number of times inspectEntry was called.
     */
    private int numInspected;

    /**
     * Constructor.
     *
     * @param geneTranscriptFile the gene-transcripts-map file to read
     * @param kVal the k value for the filter
     * @throws java.io.FileNotFoundException if the gene-transcripts-map didn't exist
     */
    public TranscriptsAlignmentFilter(final String geneTranscriptFile, final int kVal)
            throws FileNotFoundException {
        LOG.debug("** TRANSCRIPT MODE **");
        this.gtr = new GeneTranscriptRelationships();
        transcriptsIndexedIdentifiers = gtr.load(geneTranscriptFile);
        readIndexToGeneIdSetMap = new Int2ObjectOpenHashMap<IntSet>();
        k = kVal;
        numInspected = 0;

    }


    /**
     * Set the new / updated header that is being used when filtering these entries.
     *
     * @param targets targets of the merged alignment.
     */
    @Override
    public void setHeader(final IndexedIdentifier targets) {
        final ObjectSet<MutableString> targetIds = targets.keySet();
        // Use an int[] for this
        alignmentReferenceIndexToTranscriptIndex = new int[targetIds.size()];
        int i = 0;
        for (final MutableString id : targetIds) {
            assert transcriptsIndexedIdentifiers.containsKey(id) : "transcript id must be defined as a target sequence.";
            alignmentReferenceIndexToTranscriptIndex[i++] = transcriptsIndexedIdentifiers.get(id);
        }
    }

    /**
     * Post processing (after all inspectEntry's have been called).
     */
    @Override
    public void postProcessing() {
    }


    /**
     * Inspect an entry (will be called during a first pass of reading the entries).
     *
     * @param entry the entry
     */
    @Override
    public void inspectEntry(final Alignments.AlignmentEntry entry) {
        numInspected++;
        final int readNameIndex = entry.getQueryIndex();
        IntSet geneIdSet = readIndexToGeneIdSetMap.get(readNameIndex);
        if (geneIdSet != null && geneIdSet.size() > k) {
            //     System.out.println("read id found to match more than k genes..");
            return;
        }
        final int transcriptIndex = alignmentReferenceIndexToTranscriptIndex[entry.getTargetIndex()];
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
    }


    /**
     * Determine if this entry should be retained (will be called during a second
     * pass of reading the entries).
     *
     * @param entry the entry
     * @return true if it should be retained
     */
    @Override
    public boolean shouldRetainEntry(final Alignments.AlignmentEntry entry) {
        return (readIndexToGeneIdSetMap.get(entry.getQueryIndex()).size() <= k);
    }
}
