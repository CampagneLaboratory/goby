/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.alignments.Alignments;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.Collections;

/**
 * @author Fabien Campagne
 *         Date: Aug 30, 2010
 *         Time: 4:26:57 PM
 */

public class SequenceVariationPool {

    private Int2ObjectMap<ObjectArrayList<Variation>> positionToVariation = new Int2ObjectOpenHashMap<ObjectArrayList<Variation>>();
    private Int2IntMap referenceAlleleCounts[];
    private Int2ObjectMap<IntArraySet> distinctReadIndices[];
    private Int2ObjectMap<IntArrayList> variantQualityScoresAtPosition[];
    Int2ObjectMap<ObjectArrayList<ReadIndexInfo>> readIndicesAtPosition[];

    private int numberOfGroups;
    private int lastRemovedPosition;

    public SequenceVariationPool(int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.referenceAlleleCounts = new Int2IntMap[numberOfGroups];
        this.distinctReadIndices = new Int2ObjectMap[numberOfGroups];
        this.variantQualityScoresAtPosition = new Int2ObjectMap[numberOfGroups];
        this.readIndicesAtPosition = new Int2ObjectMap[numberOfGroups];

        for (int i = 0; i < numberOfGroups; i++) {
            referenceAlleleCounts[i] = new Int2IntOpenHashMap();
            referenceAlleleCounts[i].defaultReturnValue(0);


            distinctReadIndices[i] = new Int2ObjectOpenHashMap<IntArraySet>();
            variantQualityScoresAtPosition[i] = new Int2ObjectOpenHashMap<IntArrayList>();
            readIndicesAtPosition[i] = new Int2ObjectOpenHashMap<ObjectArrayList<ReadIndexInfo>>();
        }
    }

    /**
     * Reset the pool, removing all variations and associated information.
     */
    public void reset() {
        for (Int2IntMap counts : referenceAlleleCounts) {
            counts.clear();
        }

        for (Int2ObjectMap<IntArraySet> readIndices : distinctReadIndices) {
            readIndices.clear();
        }
        for (Int2ObjectMap<IntArrayList> variantScores : variantQualityScoresAtPosition) {
            variantScores.clear();
        }
        for (Int2ObjectMap<ObjectArrayList<ReadIndexInfo>> readIndices : readIndicesAtPosition) {
            readIndices.clear();
        }
        positionToVariation.clear();

    }

    /**
     * Get the count for the reference allele. This is the number of times the reference allele is observed in the alignments
     * of the specified group.
     *
     * @param position
     * @return
     */
    public int getReferenceAlleleCount(int position, int groupIndex) {
        return referenceAlleleCounts[groupIndex].get(position);
    }

    /**
     * Remove variation and associated information for the positions between the previoulsy removed position (exclusive)
     * and this new position (inclusive).
     *
     * @param position for which information will be removed from the pool.
     */
    public void removePosition(int position) {
        for (int intermediatePosition = lastRemovedPosition + 1; intermediatePosition <= position; intermediatePosition++) {
            positionToVariation.remove(intermediatePosition);
            for (Int2IntMap counts : referenceAlleleCounts) {
                counts.remove(intermediatePosition);
            }
            for (Int2ObjectMap<IntArraySet> readIndices : distinctReadIndices) {
                readIndices.remove(intermediatePosition);
            }
            for (Int2ObjectMap<IntArrayList> map : variantQualityScoresAtPosition) {
                map.remove(intermediatePosition);
            }
            for (Int2ObjectMap<ObjectArrayList<ReadIndexInfo>> readIndices : readIndicesAtPosition) {
                readIndices.remove(intermediatePosition);
            }
        }

        lastRemovedPosition = position;
    }

    /**
     * Return the number of distinct read indices supporting variations at this position.
     *
     * @param position   position of the variation
     * @param groupIndex Index of the group in which the alignment was found.
     * @return The number of distinct positions in the read supporting the variation at the specified position.
     */
    public int getNumDistinctReadIndices(int position, int groupIndex) {
        IntArraySet readIndices = distinctReadIndices[groupIndex].get(position);
        return readIndices == null ? 0 : readIndices.size();
    }

    /**
     * Return the distinct read indices supporting variations at this position.
     *
     * @param position   position of the variation
     * @param groupIndex Index of the group in which the alignment was found.
     * @return The distinct positions in the read supporting the variation at the specified position.
     */
    public IntArraySet getDistinctReadIndices(int position, int groupIndex) {
        IntArraySet readIndices = distinctReadIndices[groupIndex].get(position);
        return readIndices;
    }

    /**
     * Return the list of variant quality scores for this position in this group of samples.
     *
     * @param position   Position of the variation.
     * @param groupIndex Index of the group that contain the variation associated with these quality scores.
     * @return
     */
    public IntArrayList getVariantQualityScores(int position, int groupIndex) {
        return variantQualityScoresAtPosition[groupIndex].get(position);

    }

    /**
     * Return the average variant quality score for this position in this group of samples.
     *
     * @param position   Position of the variation.
     * @param groupIndex Index of the group that contain the variation associated with these quality scores.
     * @return Average quality score (on a Phred scale).
     */
    public double getAverageVariantQualityScore(int position, int groupIndex) {
        IntArrayList qualityScores = getVariantQualityScores(position, groupIndex);
        if (qualityScores == null) return Double.NaN;
        double average = 0;
        for (int qs : qualityScores) {
            average += qs;
        }
        return average / qualityScores.size();
    }

    /**
     * Return the list of read indices that aligned to this position of the alignment,
     * in the samples belonging to the specified group.
     *
     * @param position   Position of the alignment that is being queried.
     * @param groupIndex Index of the group to which the alignment must belong to be considered.
     * @return
     */
    public ObjectArrayList<ReadIndexInfo> getReadIndicesAt(int position, int groupIndex) {
        return readIndicesAtPosition[groupIndex].get(position);
    }

    public class Variation {
        public int targetIndex;
        public int positionOnReference;
        public String from;
        public String to;

        public int readIndex;
        public int groupIndex;
        public byte[] qualityScores;
    }

    public ObjectArrayList<Variation> getAtPosition(int position) {
        return positionToVariation.get(position);

    }

    public void store(Alignments.AlignmentEntry entry, int groupIndex, int readerIndex) {
        // increment the count of the minor alleles:
        final int alignStart = entry.getPosition();
        final int alignEnd = alignStart + entry.getQueryLength() + entry.getNumberOfIndels();
        final Int2IntMap referenceAlleleCounts = this.referenceAlleleCounts[groupIndex];
        //   System.out.printf("Increment ref allele [%d-%d] for group %d %n", alignStart, alignEnd, groupIndex);

        int matchIndex = 0;
        final int queryLength = alignEnd - alignStart;
        for (int position = alignStart; position <= alignEnd; position++) {
            incrementReferenceAlleleCount(referenceAlleleCounts, position);
            addReadIndex(entry.getMatchingReverseStrand(), position, groupIndex,
                    matchIndex, queryLength, readerIndex);
            ++matchIndex;

        }
        for (Alignments.SequenceVariation sequenceVariation : entry.getSequenceVariationsList()) {
            String from;
            String to;
            from = sequenceVariation.getFrom();
            to = sequenceVariation.getTo();
            for (int i = 0; i < Math.max(from.length(), to.length()); i++) {

                int indexPosition = entry.getPosition() + sequenceVariation.getPosition() + i;
                ObjectArrayList<Variation> list = positionToVariation.get(indexPosition);
                if (list == null) {
                    list = new ObjectArrayList<Variation>();
                    positionToVariation.put(indexPosition, list);
                }

                Variation var = new Variation();
                var.targetIndex = entry.getTargetIndex();
                var.positionOnReference = indexPosition;
                var.readIndex = sequenceVariation.getReadIndex();
                var.from = i < from.length() ? from.substring(i) : "";
                var.to = i < to.length() ? to.substring(i) : "";
                var.groupIndex = groupIndex;
                if (sequenceVariation.hasToQuality()) {
                    var.qualityScores = sequenceVariation.getToQuality().toByteArray();
                    addVariantQualityScoreAtPosition(indexPosition, groupIndex, var.qualityScores[i]);

                }
                addReadIndexAtPosition(indexPosition, groupIndex, var.readIndex);

                // since we observe a variation at this position, we decrement the reference allele count that already
                // counted all positions over the read:
                decrementReferenceAlleleCount(referenceAlleleCounts, var.positionOnReference);
                list.add(var);
            }


        }
    }


    private void addReadIndex(boolean matchingReverseStrand, int position, int groupIndex,
                              int matchIndex, int queryLength, int alignmentReaderIndex) {
        ObjectArrayList<ReadIndexInfo> readIndicesList = readIndicesAtPosition[groupIndex].get(position);
        if (readIndicesList == null) {
            readIndicesList = new ObjectArrayList<ReadIndexInfo>();
            readIndicesAtPosition[groupIndex].put(position, readIndicesList);
        }
        int readIndex = matchingReverseStrand ? queryLength - matchIndex - 1 : matchIndex;
        readIndicesList.add(new ReadIndexInfo(readIndex, alignmentReaderIndex));
    }

    public class ReadIndexInfo {
        ReadIndexInfo(int readIndex, int alignmnentReaderIndex) {
            this.readIndex = readIndex;
            this.alignmnentReaderIndex = alignmnentReaderIndex;
        }

        /**
         * The position in the read that mapped to the reference at a given genomic location.
         */
        public int readIndex;
        /**
         * The alignment reader index, that is the alignment file that contained the alignment.
         */
        public int alignmnentReaderIndex;
    }

    private void addVariantQualityScoreAtPosition(int indexPosition, int groupIndex, byte qualityScore) {
        IntArrayList qualityScores = variantQualityScoresAtPosition[groupIndex].get(indexPosition);
        if (qualityScores == null) {
            qualityScores = new IntArrayList();
            variantQualityScoresAtPosition[groupIndex].put(indexPosition, qualityScores);
        }
        qualityScores.add(qualityScore);
    }

    private void addReadIndexAtPosition(int indexPosition, int groupIndex, int readIndex) {
        IntArraySet readIndices = distinctReadIndices[groupIndex].get(indexPosition);
        if (readIndices == null) {
            readIndices = new IntArraySet();
            distinctReadIndices[groupIndex].put(indexPosition, readIndices);
        }
        readIndices.add(readIndex);
    }

    private void decrementReferenceAlleleCount(final Int2IntMap referenceAlleleCounts, final int position) {
        //   System.out.printf("Decrementing position %d%n", position);
        int newCount = referenceAlleleCounts.get(position) - 1;
        assert newCount >= 0 : "reference allele count must be positive.";
        /* if (newCount < 0) {
            System.err.println("STOP!");
        }*/
        referenceAlleleCounts.put(position, newCount);
    }

    private void incrementReferenceAlleleCount(Int2IntMap referenceAlleleCounts, int position) {
        int newCount = referenceAlleleCounts.get(position) + 1;
        referenceAlleleCounts.put(position, newCount);
    }


}
