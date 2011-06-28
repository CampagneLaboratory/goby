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

import com.google.protobuf.ByteString;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.processors.AlignmentProcessorFactory;
import edu.cornell.med.icb.goby.alignments.processors.AlignmentProcessorInterface;
import edu.cornell.med.icb.goby.alignments.processors.DefaultAlignmentProcessorFactory;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * A helper class to iterate through a set of sorted alignments in position order. The class supports processing
 * only a subset of positions from each alignment. Since Goby 1.8.
 *
 * @author Fabien Campagne
 *         Date: Sept 3, 2010
 */
public abstract class IterateSortedAlignments<T> {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(IterateSortedAlignments.class);

    private boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();
    private DoubleIndexedIdentifier referenceIds;
    protected int lastRemovedPosition = -1;
    private int numAlignmentEntries;
    private String startOffsetArgument;
    private String endOffsetArgument;
    private int startFlapLength;

    private AlignmentReaderFactory alignmentReaderFactory = new DefaultAlignmentReaderFactory();


    /**
     * Configure the alignment factory used by this class. If this setter is not called, the DefaultAlignmentProcessorFactory is used, resulting in
     * no change to the alignment.
     *
     * @param alignmentProcessorFactory what to set the factory to.
     */
    public void setAlignmentProcessorFactory(final AlignmentProcessorFactory alignmentProcessorFactory) {
        this.alignmentProcessorFactory = alignmentProcessorFactory;
    }

    private AlignmentProcessorFactory alignmentProcessorFactory = new DefaultAlignmentProcessorFactory();

    /**
     * Set the length of the start flap. If length is larger than zero, the iterator will start reading at position
     * start - length.
     *
     * @param length Length of the start flap.
     */
    public void setStartFlapLength(int length) {
        this.startFlapLength = length;
    }

    /**
     * Parse the string of reference sequences to process during the iteration. The JSAP
     * argument include-reference-names must be defined.
     *
     * @param jsapResult The jsapResult available to the mode.
     */
    public void parseIncludeReferenceArgument(final JSAPResult jsapResult) {

        final String includeReferenceNameCommas = jsapResult.getString("include-reference-names");
        parseIncludeReferenceArgument(includeReferenceNameCommas);
        startOffsetArgument = jsapResult.getString("start-position");
        endOffsetArgument = jsapResult.getString("end-position");
    }

    /**
     * Set the start position argument. The iterator will start iterating at the specified position.
     * Format is either abolute byte position or ref-id,position-in-ref
     *
     * @param arg start position argument
     */
    public void setStartPositionArgument(String arg) {
        startOffsetArgument = arg;
    }

    /**
     * Set the factory that should be used when creating alignment readers. Use this setter to install a factory that
     * filters ambiguous reads
     *
     * @param factory alignment reader factory.
     */
    public void setAlignmentReaderFactory(AlignmentReaderFactory factory) {
        alignmentReaderFactory = factory;
    }

    /**
     * Set the end position argument. The iterator will stop iterating after the specified position.
     * Format is either abolute byte position or ref-id,position-in-ref
     *
     * @param arg end position argument
     */
    public void setEndPositionArgument(String arg) {
        endOffsetArgument = arg;
    }

    /**
     * Parse the string of reference sequences to include the iteration. The string must be a coma
     * separated list of reference identifiers.
     *
     * @param includeReferenceNameCommas names of references, separated by commas, to include in the iteration.
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

    public int getNumAlignmentEntries() {
        return numAlignmentEntries;
    }

    protected int startPosition;
    protected int endPosition;
    protected int startReferenceIndex;
    protected int endReferenceIndex;

    /**
     * Determine if a position is within the start flap (defined by startFlapLength and the slice start position).
     *
     * @param referenceIndex Index of the reference sequence for the position.
     * @param position       Position within the sequence identified by referenceIndex
     * @return True if a position is in the flap, false otherwise.
     */
    public boolean isWithinStartFlap(int referenceIndex, int position) {
        if (referenceIndex == startReferenceIndex) {
            if (position < startPosition - startFlapLength) {
                System.out.println("returning flap true");
                System.out.println("returning flap true");
                return true;
            }
        }
        return false;
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
     * @throws java.io.IOException If an error occured reading the input alignment.
     */
    public void iterate(final String... basenames) throws IOException {
        ConcatSortedAlignmentReader sortedReaders = new ConcatSortedAlignmentReader(
                false, basenames);

        final int numberOfReferences = sortedReaders.getNumberOfTargets();

        referenceIds = new DoubleIndexedIdentifier(sortedReaders.getTargetIdentifiers());

        if (referenceIds == null) {
            // no reference Ids, process all the sequences.
            filterByReferenceNames = false;

        }
        sortedReaders.close();
        LOG.info(String.format("Alignment contains %d reference sequences", numberOfReferences));
        processNumberOfReferences(numberOfReferences);
        //  CountsWriter writers[] = new CountsWriter[numberOfReferences];
        IntSortedSet referencesToProcess = new IntLinkedOpenHashSet();


        // setup referencesToProcess data structure according to the command line (filterByReferenceNames and includeReferenceNames)
        for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {


            if (filterByReferenceNames) {
                final MutableString referenceId = referenceIds.getId(referenceIndex);
                assert referenceId != null : "reference id cannot be null for reference index=" + referenceIndex;
                final String referenceName = referenceId.toString();

                if (includeReferenceNames.contains(referenceName)) {
                    // subset of reference names selected by the command line:
                    referencesToProcess.add(referenceIndex);
                }
            } else {
                // process each sequence:
                referencesToProcess.add(referenceIndex);
            }
        }

        try {
            if (StringUtils.isEmpty(startOffsetArgument) && StringUtils.isEmpty(endOffsetArgument)) {
                sortedReaders = new ConcatSortedAlignmentReader(alignmentReaderFactory, false, basenames);
            } else {
                final String[] startTokens = startOffsetArgument.split("[,]");
                final String[] endTokens = endOffsetArgument.split("[,]");
                startPosition = Integer.parseInt(startTokens[1]);
                endPosition = Integer.parseInt(endTokens[1]);

                startReferenceIndex = referenceIds.getIndex(startTokens[0]);
                endReferenceIndex = referenceIds.getIndex(endTokens[0]);
                sortedReaders = new ConcatSortedAlignmentReader(alignmentReaderFactory,
                        false,
                        basenames,
                        startReferenceIndex,
                        Math.max(0, startPosition - startFlapLength),
                        endReferenceIndex,
                        endPosition);


                // adjust referenceIndex to contain only integers between start and end (inclusive):
                for (int referenceIndex = 0; referenceIndex < referencesToProcess.size(); referenceIndex++) {
                    if (referenceIndex < startReferenceIndex || referenceIndex > endReferenceIndex) {
                        referencesToProcess.rem(referenceIndex);
                    }
                }


            }
        } catch (NumberFormatException e) {
            System.err.println("An error occured parsing --start-position or --end-position. These arguments expect \n" +
                    "a string in the format ref-id,ref-position, where ref-id is a reference identifier \n" +
                    "string and ref-position in an integer that encodes a position within the reference sequence.");
            throw e;
        }
        // track the origin of each sample entry to the reader of origin:
        sortedReaders.setAdjustSampleIndices(true);


        Alignments.AlignmentEntry alignmentEntry;
        // the first reference that we should skip to:
        int currentMinTargetIndex = referencesToProcess.firstInt();
        // skip to will go to the next entry in or after currentMinTargetIndex with at least position 0
        int lastPosition = -1;
        int lastTarget = -1;
        Int2ObjectMap<T> positionToBases = new Int2ObjectOpenHashMap<T>();


        int currentPosition;
        boolean first = true;
        ProgressLogger pg = new ProgressLogger(LOG);
        pg.start();

        final AlignmentProcessorInterface realigner = alignmentProcessorFactory.create(sortedReaders);

        realigner.setGenome(getGenome());
        while ((alignmentEntry = realigner.nextRealignedEntry(currentMinTargetIndex, 0)) != null) {

            pg.lightUpdate();
            numAlignmentEntries = advanceReference(numAlignmentEntries);
            final int referenceIndex = alignmentEntry.getTargetIndex();
            if (lastTarget != -1 && referenceIndex != lastTarget) {
                // we switch to a new reference. Cleanup any previous
                processAllPreviousPositions(lastTarget, positionToBases);
            }
            int queryLength = alignmentEntry.getQueryLength();
            assert queryLength != 0 : "queryLength should never be zero";
            currentPosition = alignmentEntry.getPosition();
            boolean forwardStrand = !alignmentEntry.getMatchingReverseStrand();
            if (lastRemovedPosition == -1) {
                lastRemovedPosition = currentPosition - 1;
            }
            if (!first && currentPosition != lastPosition) {
                // we process the list of PositionBaseInfo that map to the previously visited position:
                //      processPositions(lastPosition, positionToBases.get(lastPosition));

                processAndCleanup(referenceIndex, lastPosition, positionToBases);

                lastPosition = currentPosition;

            }
            {
                first = false;
                assert queryLength != 0 : "queryLength cannot be zero to iterate sorted alignments.";
                int currentReadIndex = forwardStrand ? 0 : (queryLength + 1);
                int currentRefPosition = alignmentEntry.getPosition() - alignmentEntry.getQueryPosition();

                int numInsertions = 0;
                int numDeletions = 0;
                List<Alignments.SequenceVariation> seqVars = alignmentEntry.getSequenceVariationsList();
                for (Alignments.SequenceVariation var : seqVars) {
                    final String from = var.getFrom();
                    final int fromLength = from.length();
                    final String to = var.getTo();
                    final int toLength = from.length();
                    final int sequenceVariationLength = Math.max(fromLength, toLength);

                    for (int i = 0; i < sequenceVariationLength; i++) {
                        final char fromChar = i >= fromLength ? '-' : from.charAt(i);
                        final char toChar = i >= toLength ? '-' : to.charAt(i);
                        if (fromChar == '-') {
                            numInsertions++;
                        }
                        if (toChar == '-') {
                            numDeletions++;
                        }
                    }
                }

                final int leftPadding = alignmentEntry.getQueryPosition();
                final int rightPadding = (queryLength + numDeletions) -
                        (alignmentEntry.getTargetAlignedLength() + numInsertions) - leftPadding;

                if (leftPadding > 0) {
                    if (LOG.isDebugEnabled()) {
                        LOG.debug(String.format("queryIndex=%d, left padding, %d bases",
                                alignmentEntry.getQueryIndex(), leftPadding));
                    }
                    for (int i = 0; i < leftPadding; i++) {
                        currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                        currentRefPosition = advanceReference(currentRefPosition);
                        // Don't observe during padding chars
                    }
                }
                int numObservedBases = 0;
                for (final Alignments.SequenceVariation var : seqVars) {
                    final String to = var.getTo();
                    final String from = var.getFrom();
                    final ByteString qualityScores = var.getToQuality();

                    final int fromLength = from.length();
                    final int toLength = to.length();
                    final int qualLength = qualityScores.size();
                    final int sequenceVariationLength = Math.max(fromLength, toLength);


                    final int preSeqvarBases;
                    if (from.charAt(0) == '-') {
                        preSeqvarBases = var.getPosition() - numObservedBases;
                    } else {
                        preSeqvarBases = var.getPosition() - numObservedBases - 1;
                    }
                    for (int i = 0; i < preSeqvarBases; i++) {
                        // Bases before the next sequence variation
                        currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                        currentRefPosition = advanceReference(currentRefPosition);
                        observeReferenceBase(sortedReaders, alignmentEntry, positionToBases,
                                referenceIndex, currentRefPosition, currentReadIndex);
                        numObservedBases++;
                    }
                    for (int i = 0; i < sequenceVariationLength; i++) {
                        /*------------------------------------------------------------------
                         * For details on how to count refPosition and readIndex, especially
                         * with respect to DELETIONS and INSERTIONS see
                         *
                         *    http://tinyurl.com/goby-sequence-variations
                         *
                         *------------------------------------------------------------------*/

                        // Bases within the sequence variation
                        final char toChar = i >= toLength ? '-' : to.charAt(i);
                        final char fromChar = i >= fromLength ? '-' : from.charAt(i);
                        final byte toQual = i >= qualLength ? 0 : qualityScores.byteAt(i);
                        if (fromChar == '-') {
                            // During an insert, do not increment refPosition 
                        } else {
                            numObservedBases++;
                            currentRefPosition = advanceReference(currentRefPosition);
                        }

                        if (toChar == '-') {
                            if (forwardStrand) {
                                // Do no increment readIndex during a delete on forward strand
                            } else if (i == 0) {
                                // On reverse strand delete, decrement readIndex for the first base ONLY
                                currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                            } else {
                                // On reverse strand deletion, after the first base, do not decrement readIndex
                            }
                        } else {
                            currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                        }

                        observeVariantBase(sortedReaders, alignmentEntry, positionToBases,
                                var, toChar, fromChar, toQual,
                                referenceIndex, currentRefPosition, currentReadIndex);

                        if (toChar == '-' && !forwardStrand && i == sequenceVariationLength - 1) {
                            // The logic in this algorithm is (increment/decrement) before observe().
                            // After a deletion on reverse strand we want the next base after the deletion
                            // to have the same readIndex as the deletion bases. Here we will
                            // increment the readIndex by one. When the next base is processed the readIndex
                            // will be decremented and thus be the same readIndex as those for the deletion.
                            currentReadIndex = advanceReadIndex(!forwardStrand, currentReadIndex);
                        }

                    }
                    //
                    if (var.getFrom().indexOf('-') >= 0 || var.getTo().indexOf('-') >= 0) {
                        observeIndel(positionToBases, referenceIndex,
                                alignmentEntry.getPosition()+var.getPosition(),
                                var.getFrom(), var.getTo(),
                                alignmentEntry.getSampleIndex()
                        );
                    }
                }

                while (forwardStrand ? currentReadIndex < queryLength - rightPadding :
                        currentReadIndex > 1 + rightPadding) {

                    // match stretch before next variation / end of read
                    currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                    currentRefPosition = advanceReference(currentRefPosition);

                    assert currentReadIndex >= 1 && currentReadIndex < queryLength + 1 :
                            String.format("currentReadIndex %d is out of range.", currentReadIndex);

                    observeReferenceBase(sortedReaders, alignmentEntry, positionToBases,
                            referenceIndex, currentRefPosition, currentReadIndex);
                }
                if (rightPadding > 0) {
                    LOG.debug(String.format("queryIndex=%d, right padding, %d bases",
                            alignmentEntry.getQueryIndex(), rightPadding));
                }
            }

            if (referencesToProcess.contains(referenceIndex)) {


                lastPosition = alignmentEntry.getPosition();
                lastTarget = alignmentEntry.getTargetIndex();
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
        int minPos = Integer.MAX_VALUE;
        int maxPos = Integer.MIN_VALUE;

        for (int pos : positionToBases.keySet()) {
            minPos = Math.min(pos, minPos);
            maxPos = Math.max(pos, maxPos);
        }

        for (int position = minPos; position <= maxPos; position++) {
            processAndCleanup(lastTarget, position, positionToBases);
        }

        sortedReaders.close();
        pg.stop();
    }


    private int advanceReadIndex(boolean forwardStrand, int currentReadIndex) {
        currentReadIndex += forwardStrand ? 1 : -1;
        // System.out.printf(" read-index=%d ", currentReadIndex);

        return currentReadIndex;
    }

    private int advanceReference(int currentRefPosition) {
        currentRefPosition += 1;
        //   System.out.printf(" ref-pos=%d %n", currentRefPosition);

        return currentRefPosition;
    }

    private void processAndCleanup(final int lastReferenceIndex,
                                   final int lastPosition, final
    Int2ObjectMap<T> positionToBases) {

        for (int intermediatePosition = lastRemovedPosition + 1;
             intermediatePosition <= lastPosition; intermediatePosition++) {

            if (positionToBases.containsKey(intermediatePosition)) {

                processPositions(lastReferenceIndex, intermediatePosition, positionToBases.get(intermediatePosition));
                positionToBases.remove(intermediatePosition);
            }

        }
        lastRemovedPosition = lastPosition;
    }

    /**
     * Temporary position used for sorting in method processAllPreviousPositions.
     */
    IntArrayList tmpPositions = new IntArrayList();

    /**
     * Process positions on the previous target, which may still be in positionToBases. Note that this method is no
     * re-entrant.
     *
     * @param lastReferenceIndex the last referenceIndex?
     * @param positionToBases    positionToBases?
     */
    private void processAllPreviousPositions(final int lastReferenceIndex, final Int2ObjectMap<T> positionToBases) {

        tmpPositions.clear();
        tmpPositions.addAll(positionToBases.keySet());
        Collections.sort(tmpPositions);


        for (final int intermediatePosition : tmpPositions) {
            if (positionToBases.containsKey(intermediatePosition)) {

                processPositions(lastReferenceIndex, intermediatePosition, positionToBases.get(intermediatePosition));
                positionToBases.remove(intermediatePosition);
                lastRemovedPosition = intermediatePosition;
            }
        }
        positionToBases.clear();
    }

    /**
     * Implement this call-back method to observe reference bases.
     *
     * @param sortedReaders         The concat read that contains the variation.
     * @param alignmentEntry        The alignment entry that contains the variation
     * @param positionToBases       Map keyed by reference position, used to accumulate information for each position.
     * @param currentReferenceIndex Index of the reference sequence where the variant occurs.
     * @param currentRefPosition    Position where the variant occurs in the reference.
     * @param currentReadIndex      Index in the read where the variant occurs.
     */
    public abstract void observeReferenceBase(ConcatSortedAlignmentReader sortedReaders,
                                              Alignments.AlignmentEntry alignmentEntry,
                                              Int2ObjectMap<T> positionToBases,
                                              int currentReferenceIndex,
                                              int currentRefPosition,
                                              int currentReadIndex);

    /**
     * Implement this call-back method to observe variant bases.
     *
     * @param sortedReaders         The concat read that contains the variation.
     * @param alignmentEntry        The alignment entry that contains the variation
     * @param positionToBases       Map keyed by reference position, used to accumulate information for each position.
     * @param var                   The sequence variation from the alignment entry that triggered emiting this observation.
     * @param toChar                The base character in the read
     * @param fromChar              The base character in the reference.
     * @param toQual                The base quality value in the read (if it exists, otherwise 0)
     * @param currentReferenceIndex Index of the reference sequence where the variant occurs.
     * @param currentRefPosition    Position where the variant occurs in the reference.
     * @param currentReadIndex      Index in the read where the variant occurs.
     */
    public abstract void observeVariantBase(ConcatSortedAlignmentReader sortedReaders,
                                            Alignments.AlignmentEntry alignmentEntry,
                                            Int2ObjectMap<T> positionToBases,
                                            Alignments.SequenceVariation var,
                                            char toChar, char fromChar,
                                            byte toQual, int currentReferenceIndex,
                                            int currentRefPosition,
                                            int currentReadIndex);


    public abstract void processPositions(int referenceIndex, int intermediatePosition, T positionBaseInfos);

    /**
     * Implement this call-back method to observe a candidate indel that begins at startPosition.
     *
     * @param positionToBases map from eir start positions to info about indels whose eir start at the location.
     * @param referenceIndex  The reference sequence where the indel candidate is observed.
     * @param startPosition   The position where the indel start, according to pair-wise sequence alignment (before eir calculation)
     * @param from            from bases (before eir calculation)
     * @param to              to bases (before eir calculation)
     * @param sampleIndex     Index of the sample where the indel was observed.
     */
    public void observeIndel(final Int2ObjectMap<T> positionToBases,
                             final int referenceIndex,
                             final int startPosition,
                             final String from, final String to,
                             final int sampleIndex) {
    }


    /**
     * Will be called to let the client know how many references will be processed.
     *
     * @param numberOfReferences The number of references in this alignment.
     * @throws java.io.IOException If an error occcurs.
     */
    public void processNumberOfReferences(final int numberOfReferences) throws IOException {
    }


    /**
     * Return the reference sequence id given its index.
     *
     * @param targetIndex The index of the desired reference sequence
     * @return The id of the reference sequence
     */
    public CharSequence getReferenceId(final int targetIndex) {
        return referenceIds.getId(targetIndex);
    }

    public RandomAccessSequenceInterface getGenome() {
        return null;
    }
}