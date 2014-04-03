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

import edu.cornell.med.icb.goby.util.WarningCounter;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import org.apache.log4j.Logger;

/**
 * @author Fabien Campagne
 *         Date: Sep 7, 2010
 *         Time: 2:14:38 PM
 */
public abstract class IterateSortedAlignmentsListImpl
        extends IterateSortedAlignments<DiscoverVariantPositionData> {


    protected IterateSortedAlignmentsListImpl() {
        this.SUB_SAMPLE_SIZE = doc.getInteger("sub-sample-size");

    }

    public static final DynamicOptionClient doc() {
        return doc;
    }

    @RegisterThis
    public static final DynamicOptionClient doc = new DynamicOptionClient(IterateSortedAlignmentsListImpl.class,

            /*"max-coverage:integer, The maximum number of bases allowed before triggering sub-sampling. Some positions just " +
                    "have too much coverage, usually because of some artefactual effect. We use this parameter to determine " +
                    "when too much is too much. Any position with more covering bases will be sub-sampled to reduce computational load " +
                    "and make performance predictable. The default parameter is set at half a million bases.:500000" +*/
            "sub-sample-size:integer, The number of bases to keep when coverage exceeds the maximum and sub-sampling " +
                    "is needed to improve performance.:10000"
    );
    /**
     * Used to log debug and informational messages.
     */

    private static final Logger LOG = Logger.getLogger(IterateSortedAlignmentsListImpl.class);

    /**
     * Process a list of bases at a given reference position.
     *
     * @param referenceIndex       Index of the reference sequence where these bases align.
     * @param intermediatePosition Position is zero-based
     * @param positionBaseInfos    List of base information for bases that aligned.
     */
    @Override
    public abstract void processPositions(int referenceIndex, int intermediatePosition, DiscoverVariantPositionData positionBaseInfos);


    @Override
    public final void observeReferenceBase(final ConcatSortedAlignmentReader sortedReaders,
                                           final Alignments.AlignmentEntry alignmentEntry,
                                           final PositionToBasesMap<DiscoverVariantPositionData> positionToBases,
                                           final int currentReferenceIndex, final int currentRefPosition, final int currentReadIndex) {
        /*if (LOG.isTraceEnabled()) {
            LOG.trace(String.format("RB: queryIndex=%d\tref_position=%d\tread_index=%d",
                    alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex));
        } */

        final PositionBaseInfo info = new PositionBaseInfo();

        info.readerIndex = alignmentEntry.getSampleIndex();
        //     System.out.printf("observing ref readerIndex=%d%n",info.readerIndex);
        /* if (currentRefPosition==3661601-1 && !alignmentEntry.getMatchingReverseStrand()) {
            System.out.println("STOP1");
        } */
        info.readIndex = currentReadIndex;
        info.from = '\0';
        info.to = '\0';
        info.matchesReference = true;
        info.position = currentRefPosition - 1; // store 0-based position
        info.qualityScore = (byte) (alignmentEntry.hasMappingQuality() ? alignmentEntry.getMappingQuality() : 40);
        info.matchesForwardStrand = !alignmentEntry.getMatchingReverseStrand();
        //System.out.printf("position=%d %s%n", currentRefPosition, info);
        addToFuture(positionToBases, info);
    }


    @Override
    public void observeVariantBase(final ConcatSortedAlignmentReader sortedReaders,
                                   final Alignments.AlignmentEntry alignmentEntry,
                                   final PositionToBasesMap<DiscoverVariantPositionData> positionToBases,
                                   final Alignments.SequenceVariation variation,
                                   final char toChar, final char fromChar,
                                   final byte toQual, final int currentReferenceIndex,
                                   final int currentRefPosition,
                                   final int currentReadIndex) {

      /*  if (LOG.isTraceEnabled()) {
            LOG.trace(String.format("VB: queryIndex=%d\tref_position=%d\tread_index=%d\tfromChar=%c\ttoChar=%c",
                    alignmentEntry.getQueryIndex(), currentRefPosition, currentReadIndex, fromChar, toChar));
        }
        */
        final PositionBaseInfo info = new PositionBaseInfo();
        info.readerIndex = alignmentEntry.getSampleIndex();
        //    System.out.printf("observing var readerIndex=%d%n",info.readerIndex);

        info.readIndex = currentReadIndex;
        info.from = fromChar;
        info.to = toChar;
        info.matchesReference = false;
        info.position = currentRefPosition - 1; // store 0-based position
        final int readMappingQuality = (alignmentEntry.hasMappingQuality() ? alignmentEntry.getMappingQuality() : 40);
        info.qualityScore = (byte) Math.min(toQual, readMappingQuality);
        info.matchesForwardStrand = !alignmentEntry.getMatchingReverseStrand();
        addToFuture(positionToBases, info);
    }

    private final WarningCounter moreVariantsThanThreshold = new WarningCounter(10);
    private int SUB_SAMPLE_SIZE = 10000;

    private final void addToFuture(final PositionToBasesMap<DiscoverVariantPositionData> positionToBases,
                                   final PositionBaseInfo info) {
        final int position = info.position;
        DiscoverVariantPositionData list = positionToBases.get(position);
        if (list == null) {
            list = new DiscoverVariantPositionData(position);
            positionToBases.put(position, list);
        } else {
            assert list.getZeroBasedPosition() == position : "info position must match list position.";
        }
        boolean isIgnoredPosition = positionToBases.isIgnoredPosition(position);
        if (isIgnoredPosition || list.size() >= maxThreshold) {

            // stop accumulating if the position has more than half a million bases!
            // also sub-sample the already collection bases to reduce coverage to 10,000.

            if (!isIgnoredPosition) {
                IntOpenHashSet positions = new IntOpenHashSet();
                for (PositionBaseInfo element : list) {
                    positions.add(element.position);
                }
                moreVariantsThanThreshold.warn(LOG, "position=%d has more variants %d than max threshold=%d. Stopped recording. Distinct position#%d ",
                        info.position, maxThreshold, list.size(),
                        positions.size());
                list.subSample(SUB_SAMPLE_SIZE);
            }
            positionToBases.markIgnoredPosition(position);

            return;
        }

        list.add(info);

    }


}
