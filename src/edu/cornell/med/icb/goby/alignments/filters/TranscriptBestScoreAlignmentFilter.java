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

package edu.cornell.med.icb.goby.alignments.filters;

import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;

import java.io.FileNotFoundException;

/**
 * Filter by best score and by number of genes matched.
 * @author Fabien Campagne
 * Date: May 19, 2009
 * Time: 8:46:07 AM
 *
 */
public class TranscriptBestScoreAlignmentFilter extends AbstractAlignmentEntryFilter {
    TranscriptsAlignmentFilter transcriptFilter;
    BestScoreOnlyAlignmentFilter bestScoreFilter;

    /**
     * Constructor.
     *
     * @param geneTranscriptFile the gene-transcripts-map file to read
     * @param kVal               the k value for the filter
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
        if (!bestScoreFilter.shouldRetainEntry(entry)) {
            // if the entry does not have the best score, skip

            return false;
        } else {
            // let the transcript filter decide.
            return transcriptFilter.shouldRetainEntry(entry);
        }
    }
     /**
     * Give the filter access to targets of the merged alignment.
     *
     * @param targets
     */
    @Override
    public void setHeader(final IndexedIdentifier targets) {
         transcriptFilter.setHeader(targets);
         bestScoreFilter.setHeader(targets);
    }
}
