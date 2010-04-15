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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.tissueinfo.similarity.GeneTranscriptRelationships;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.FileNotFoundException;

/**
 * This class assists with filtering
 * {@link edu.cornell.med.icb.goby.alignments.Alignments.AlignmentEntry} values when merging
 * transcript-based alignments.
 *
 * @author Kevin Dorff
 * @author Fabien Campagne
 */
public final class TranscriptsAlignmentFilter extends AbstractAlignmentEntryFilter {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TranscriptsAlignmentFilter.class);

    /**
     * The helps with relationships between genes and transcripts.
     */
    private final GeneTranscriptRelationships gtr;
    /**
     * The transcript identifiers as represented in gtr.
     */
    private final IndexedIdentifier transcriptsIndexedIdentifiers;
    /**
     * The map of read index to gene identifiers.
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
   // private DoubleIndexedIdentifier backTranscriptIds;


    /**
     * Constructor.
     *
     * @param geneTranscriptFile the gene-transcripts-map file to read
     * @param kVal               the k value for the filter
     * @throws java.io.FileNotFoundException if the gene-transcripts-map didn't exist
     */
    public TranscriptsAlignmentFilter(final String geneTranscriptFile, final int kVal)
            throws FileNotFoundException {
        super();
        LOG.debug("** TRANSCRIPT MODE **");
        this.gtr = new GeneTranscriptRelationships();
        transcriptsIndexedIdentifiers = gtr.load(geneTranscriptFile);
        readIndexToGeneIdSetMap = new Int2ObjectOpenHashMap<IntSet>();
        k = kVal;
        numInspected = 0;
     //   backTranscriptIds = new DoubleIndexedIdentifier(transcriptsIndexedIdentifiers);

    }

    /**
     * Set the new / updated header that is being used when filtering these entries.
     *
     * @param targets targets of the merged alignment.
     */
    @Override
    public void setTargetIdentifiers(final IndexedIdentifier targets) {
        final ObjectSet<MutableString> targetIds = targets.keySet();
        // Use an int[] for this
        alignmentReferenceIndexToTranscriptIndex = new int[targetIds.size()];

        for (final MutableString id : targetIds) {
            assert transcriptsIndexedIdentifiers.containsKey(id) : "transcript id must be defined as a target sequence.";
            final int transcriptIndex = transcriptsIndexedIdentifiers.getInt(id);
            final int referenceIndex = targets.getInt(id);
            alignmentReferenceIndexToTranscriptIndex[referenceIndex] = transcriptIndex;

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
        final int queryIndex = entry.getQueryIndex();
        IntSet geneIdSet = readIndexToGeneIdSetMap.get(queryIndex);
        if (geneIdSet != null && geneIdSet.size() > k) {
            //     System.out.println("read id found to match more than k genes..");
            return;
        }
        final int transcriptIndex = alignmentReferenceIndexToTranscriptIndex[entry.getTargetIndex()];
        final int geneIndex = gtr.transcript2Gene(transcriptIndex);
        //   System.out.printf("transcript %s gene %s %n", backTranscriptIds.getId(transcriptIndex), gtr.getGeneId(geneIndex));
        if (geneIdSet == null) {
            // An array set is a fast implementation when the set is expected to be small
            // (a few genes on average per read). The ArraySet is also the smallest footprint implementation,
            // this matters most here where we keep one such set for each read.
            geneIdSet = new IntArraySet();
            readIndexToGeneIdSetMap.put(queryIndex, geneIdSet);
        }
        if (geneIdSet.size() <= k) {
            // No need to add to the set if the size is already too big (save memory)
            geneIdSet.add(geneIndex);
            final boolean verbose = false;
            if (verbose && geneIdSet.size() >= 5) {
                System.out.print("gene set: ");
                for (final int geneIndex2 : geneIdSet) {
                    if (geneIndex2 == geneIndex) {
                        System.out.printf("position: %d score: %g ", entry.getPosition(), entry.getScore());
                    }
                    System.out.printf("%s ", gtr.getGeneId(geneIndex2));
                }
                System.out.println("");
            }

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
