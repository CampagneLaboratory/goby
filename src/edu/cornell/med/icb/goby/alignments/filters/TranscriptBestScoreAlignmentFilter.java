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

import java.io.FileNotFoundException;

/**
 * Filter by best score and by number of genes matched.
 *
 * @author Fabien Campagne
 *         Date: May 19, 2009
 *         Time: 8:46:07 AM
 */
public class TranscriptBestScoreAlignmentFilter extends AbstractAlignmentEntryFilter {
    private final TranscriptsAlignmentFilter transcriptFilter;
    private final BestScoreOnlyAlignmentFilter bestScoreFilter;
    private int notBestScoreCount;
    private int geneAmbiguityCount;
    private int entryCount;

    @Override
    public void printStats() {
        System.out.printf("notBestScoreCount=%g %% geneAmbiguityCount=%g %% %n",
                div(notBestScoreCount, entryCount) * 100,
                div(geneAmbiguityCount, entryCount) * 100);
    }

    private double div(final int a, final int b) {
        return ((double) a / (double) b);
    }

    /**
     * Constructor.
     *
     * @param geneTranscriptFile the gene-transcripts-map file to read
     * @param kVal               the k value for the filter
     * @param maxNumberOfReads   Maximum number of reads.
     * @throws java.io.FileNotFoundException if the gene-transcripts-map didn't exist
     */
    public TranscriptBestScoreAlignmentFilter(final String geneTranscriptFile, final int kVal,
                                              final int maxNumberOfReads) throws FileNotFoundException {
        transcriptFilter = new TranscriptsAlignmentFilter(geneTranscriptFile, kVal);
        bestScoreFilter = new BestScoreOnlyAlignmentFilter(kVal, maxNumberOfReads);
    }

    @Override
    public void inspectEntry(final Alignments.AlignmentEntry entry) {
        bestScoreFilter.inspectEntry(entry);
        transcriptFilter.inspectEntry(entry);
    }

    @Override
    public void postProcessing() {
        bestScoreFilter.postProcessing();
        transcriptFilter.postProcessing();
    }

    @Override
    public boolean shouldRetainEntry(final Alignments.AlignmentEntry entry) {
        entryCount +=  entry.getMultiplicity();
        if (!bestScoreFilter.shouldRetainEntry(entry)) {
            // if the entry does not have the best score, skip
            notBestScoreCount += entry.getMultiplicity();
            return false;
        } else {
            // the entry has best score, let the transcript filter decide.
            final boolean notGeneAmbiguous = transcriptFilter.shouldRetainEntry(entry);
            if (!notGeneAmbiguous) {
                geneAmbiguityCount += entry.getMultiplicity();
            }
            return notGeneAmbiguous;
        }
    }

    /**
     * Give the filter access to targets of the merged alignment.
     *
     * @param targets
     */
    @Override
    public void setTargetIdentifiers(final IndexedIdentifier targets) {
        transcriptFilter.setTargetIdentifiers(targets);
        bestScoreFilter.setTargetIdentifiers(targets);
    }
}
