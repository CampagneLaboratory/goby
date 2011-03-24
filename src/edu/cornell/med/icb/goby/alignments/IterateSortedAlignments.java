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
import com.google.protobuf.ByteString;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

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
    private IntSortedSet referencesToProcess;
    private DoubleIndexedIdentifier referenceIds;
    protected int lastRemovedPosition = -1;
    private Int2IntMap readIndexVariationTally;
    private int numAlignmentEntries;
    private String startOffsetArgument;
    private String endOffsetArgument;
    private int startFlapLength;

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
        ConcatSortedAlignmentReader sortedReaders = new ConcatSortedAlignmentReader(basenames);

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
        referencesToProcess = new IntLinkedOpenHashSet();


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

        sortedReaders = null;
        try {
            if (StringUtils.isEmpty(startOffsetArgument) && StringUtils.isEmpty(endOffsetArgument)) {
                sortedReaders = new ConcatSortedAlignmentReader(basenames);
            } else {
                final String[] startTokens = startOffsetArgument.split("[,]");
                final String[] endTokens = endOffsetArgument.split("[,]");
                startPosition = Integer.parseInt(startTokens[1]);
                endPosition = Integer.parseInt(endTokens[1]);

                startReferenceIndex = referenceIds.getIndex(startTokens[0]);
                endReferenceIndex = referenceIds.getIndex(endTokens[0]);
                sortedReaders = new ConcatSortedAlignmentReader(basenames, startReferenceIndex,
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


        Alignments.AlignmentEntry alignmentEntry = null;
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
        while ((alignmentEntry = sortedReaders.skipTo(currentMinTargetIndex, 0)) != null) {
            pg.lightUpdate();
            numAlignmentEntries = advanceReference(numAlignmentEntries);
            final int referenceIndex = alignmentEntry.getTargetIndex();
            if (lastTarget != -1 && referenceIndex != lastTarget) {
                // we switch to a new reference. Cleanup any previous
                processAllPreviousPositions(lastTarget, positionToBases);
            }
            int queryLength = alignmentEntry.getQueryLength();

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
                assert queryLength!=0 : "queryLength cannot be zero to iterate sorted alignments.";
                int currentReadIndex = forwardStrand ? 0 : queryLength +1;
                int currentRefPosition = alignmentEntry.getPosition() - 1;
                int lastMatchIndex = 0;

                // add to the list at this position:
                for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                    for (int i = lastMatchIndex; i < var.getPosition() - 1; i++) {
                        // match stretch before variation:
                        currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                        currentRefPosition = advanceReference(currentRefPosition);
                        observeReferenceBase(sortedReaders, alignmentEntry, positionToBases,

                                referenceIndex, currentRefPosition, currentReadIndex);
                    }

                    final String to = var.getTo();
                    final String from = var.getFrom();
                    final ByteString qualityScores= var.getToQuality();

                    final int fromLength = from.length();
                    final int toLength = to.length();
                    final int qualLength=qualityScores.size();
                    int length = Math.max(fromLength, toLength);

                    lastMatchIndex = var.getPosition() + length - 1;

                    for (int i = 0; i < length; i++) {

                        final char toChar = i >= toLength ? '-' : to.charAt(i);
                        final char fromChar = i >= fromLength ? '-' : from.charAt(i);
                    //    System.out.printf("tL: %d qualScores: %s%n", toLength, ByteArrayList.wrap(qualityScores.toByteArray()));
                        final byte toQual= i >= qualLength ? 0 : qualityScores.byteAt(i);
                        if (toChar != '-') {
                            currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                        }
                        if (fromChar != '-') {
                            currentRefPosition = advanceReference(currentRefPosition);

                        }

                        observeVariantBase(sortedReaders, alignmentEntry, positionToBases,
                                var,
                                toChar, fromChar, toQual,
                                referenceIndex, currentRefPosition, currentReadIndex);

                        lastMatchIndex = i + var.getPosition();

                    }


                }
                while (forwardStrand ? currentReadIndex < queryLength :
                        currentReadIndex > 1) {

                    // match stretch before variation:
                    currentReadIndex = advanceReadIndex(forwardStrand, currentReadIndex);
                    currentRefPosition = advanceReference(currentRefPosition);


                    assert currentReadIndex >= 1 && currentReadIndex < queryLength + 1 :
                            String.format("currentReadIndex %d is out of range.", currentReadIndex);

                    observeReferenceBase(sortedReaders, alignmentEntry, positionToBases,
                            referenceIndex, currentRefPosition, currentReadIndex);
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

    private void processAndCleanup(int lastReferenceIndex, int lastPosition, Int2ObjectMap<T> positionToBases) {
        for (int intermediatePosition = lastRemovedPosition + 1; intermediatePosition <= lastPosition; intermediatePosition++) {
            if (positionToBases.containsKey(intermediatePosition)) {
                processPositions(lastReferenceIndex, intermediatePosition, positionToBases.get(intermediatePosition));
                positionToBases.remove(intermediatePosition);
            }

        }
        lastRemovedPosition = lastPosition;
    }

    /**
     * Process positions on the previous target, which may still be in positionToBases.
     *
     * @param lastReferenceIndex
     * @param positionToBases
     */
    private void processAllPreviousPositions(int lastReferenceIndex, Int2ObjectMap<T> positionToBases) {
        IntArrayList positions = new IntArrayList();
        positions.addAll(positionToBases.keySet());
        Collections.sort(positions);


        for (int intermediatePosition : positions) {
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
     * @param toQual
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


}