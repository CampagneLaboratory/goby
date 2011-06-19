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

import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntLinkedOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Arrays;

/**
 * A helper class to iterate through a set of alignments and process only a subset of
 * references in each alignment. Since Goby 1.7, this class uses the skipTo method to
 * scan alignments that are both sorted and indexed.
 *
 * @author Fabien Campagne
 *         Date: Mar 10, 2010
 *         Time: 11:12:23 AM
 */
public abstract class IterateAlignments {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(IterateAlignments.class);

    private boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();
    private IntSortedSet referencesToProcess;
    private DoubleIndexedIdentifier referenceIds;

    /**
     * Parse the string of reference sequences to process during the iteration. The JSAP
     * argument include-reference-names must be defined.
     *
     * @param jsapResult The jsapResult available to the mode.
     */
    public void parseIncludeReferenceArgument(final JSAPResult jsapResult) {

        final String includeReferenceNameCommas = jsapResult.getString("include-reference-names");
        parseIncludeReferenceArgument(includeReferenceNameCommas);

    }

    /**
     * Parse the string of reference sequences to include the iteration. The string must be a coma
     * separated list of reference identifiers.
     */
    public void parseIncludeReferenceArgument(final String includeReferenceNameCommas) {
        if (includeReferenceNameCommas != null) {
            includeReferenceNames = new ObjectOpenHashSet<String>();
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameCommas.split("[,]")));
            LOG.info("Will iterate through the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }
    }

    /**
     * Iterate through one alignment, restricting the iteration to the alignment records 'within'
     * the bytes startOffset and endOffset. In this context, within is defined by the chunk semantic
     * of the Goby files. See {@link edu.cornell.med.icb.goby.reads.FastBufferedMessageChunksReader}
     * for details about this semantic.
     *
     * @param startOffset Start of the allowed window, in bytes in the compressed entries file.
     * @param endOffset   End of the allowed window, in bytes in the compressed entries file.
     * @param basename    Basename of the alignment to iterate over.
     * @throws IOException If an error occured reading the input alignment.
     */
    public void iterate(final long startOffset, final long endOffset, final String basename) throws IOException {
        iterateOverOneAlignment(startOffset, endOffset, basename);
    }

    /**
     * Iterate through a set of alignments. Iterations are performed through these steps:
     * <UL>
     * <LI>Iterate will call processHeader on each input alignment, giving the opportunity
     * to the client to read the header of the alignment and extract information from it.
     * </LI>
     * </UL>
     *
     * @param basenames Basenames of the alignments to iterate over.
     * @throws IOException If an error occured reading the input alignment.
     */
    public void iterate(final String... basenames) throws IOException {
        for (final String basename : basenames) {
            iterateOverOneAlignment(0, Long.MAX_VALUE, basename);
        }
    }

    private void iterateOverOneAlignment(final long startOffset, final long endOffset, final String basename) throws IOException {
        final int numberOfReferences;
        {
            final AlignmentReader reader = alignmentReaderFactory.createReader(basename, startOffset, endOffset);
            reader.readHeader();
            numberOfReferences = reader.getNumberOfTargets();

            referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
            reader.close();
        }
        LOG.debug(String.format("Alignment contains %d reference sequences", numberOfReferences));
        processNumberOfReferences(basename, numberOfReferences);
        //  CountsWriter writers[] = new CountsWriter[numberOfReferences];
        referencesToProcess = new IntLinkedOpenHashSet();

        // setup referencesToProcess data structure according to the command line (filterByReferenceNames and includeReferenceNames)
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {

            final MutableString referenceId = referenceIds.getId(referenceIndex);
            assert referenceId != null : "reference id cannot be null for reference index=" + referenceIndex;
            final String referenceName = referenceId.toString();
            if (filterByReferenceNames) {
                if (includeReferenceNames.contains(referenceName)) {
                    // subset of reference names selected by the command line:
                    referencesToProcess.add(referenceIndex);
                }
            } else {
                // process each sequence:
                referencesToProcess.add(referenceIndex);
            }
        }

        final AlignmentReader alignmentReader = alignmentReaderFactory.createReader(basename, startOffset, endOffset);
        alignmentReader.readHeader();

        // Give the client the ability to prepare data structures for each reference that will be processed.
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
            if (referencesToProcess.contains(referenceIndex)) {
                prepareDataStructuresForReference(alignmentReader, referenceIndex);
            }
        }

        // read the alignment:

        LOG.debug("Loading the alignment " + basename);
        if (alignmentReader.isSorted()) {
            LOG.debug("The alignment is sorted, iteration will use the faster skipTo method.");
            // the alignment is not sorted, we leverage skipTo to get directly to the sequence of interest.:

            Alignments.AlignmentEntry alignmentEntry = null;
            // the first reference that we should skip to:
            int currentMinTargetIndex = referencesToProcess.firstInt();
            // skip to will go to the next entry in or after currentMinTargetIndex with at least position 0
            while ((alignmentEntry = alignmentReader.skipTo(currentMinTargetIndex, 0)) != null) {
                final int referenceIndex = alignmentEntry.getTargetIndex();
                if (referencesToProcess.contains(referenceIndex)) {
                    processAlignmentEntry(alignmentReader, alignmentEntry);
                }
                if (referenceIndex > currentMinTargetIndex) {
                    // we are past the reference sequence we use to skip to
                    // Check if we are done:
                    final boolean success = referencesToProcess.remove(currentMinTargetIndex);
                    assert success : "removing an element from referencesToProcess must succeed. ";
                    if (referencesToProcess.isEmpty()) {
                        // we are done.
                        break;
                    }
                    // Not done, we now look for the next requested reference sequence:
                    currentMinTargetIndex = referencesToProcess.firstInt();
                }
            }
        } else {
            // the alignment is not sorted, we cannot use skipTo:

            for (final Alignments.AlignmentEntry alignmentEntry : alignmentReader) {
                final int referenceIndex = alignmentEntry.getTargetIndex();
                if (referencesToProcess.contains(referenceIndex)) {
                    processAlignmentEntry(alignmentReader, alignmentEntry);
                }
            }
        }
        alignmentReader.close();
    }

    /**
     * Process one alignment entry.
     *
     * @param alignmentReader The reader that parsed this entry.
     * @param alignmentEntry  The parsed entry.
     */
    public abstract void processAlignmentEntry(AlignmentReader alignmentReader, Alignments.AlignmentEntry alignmentEntry);

    /**
     * Called to let the subclass prepare some datastructure for each reference sequence. The method is called
     * exactly once for each valid reference index.
     *
     * @param alignmentReader The reader currently being iterated over.
     * @param referenceIndex  The index of the reference sequence for which data structures should be initialized.
     */
    public void prepareDataStructuresForReference(final AlignmentReader alignmentReader, final int referenceIndex) {

    }

    /**
     * Will be called to let the client know how many references will be processed in a given alignment.
     *
     * @param basename           The basename of the alignment being processed.
     * @param numberOfReferences The number of references in this alignment.
     * @throws IOException If an error occcurs.
     */
    public void processNumberOfReferences(final String basename, final int numberOfReferences) throws IOException {
    }


    /**
     * Return the reference sequence id given its index.
     *
     * @param targetIndex The index of the desired reference sequence
     * @return The id of the reference sequence
     */
    protected CharSequence getReferenceId(final int targetIndex) {
        return referenceIds.getId(targetIndex);
    }

    private AlignmentReaderFactory alignmentReaderFactory;

    public void setAlignmentReaderFactory(AlignmentReaderFactory factory) {
        this.alignmentReaderFactory = factory;
    }

}
