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

import edu.cornell.med.icb.goby.Release1_9_7_2;
import edu.cornell.med.icb.goby.algorithmic.algorithm.EquivalentIndelRegionCalculator;
import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import edu.cornell.med.icb.goby.alignments.processors.ObservedIndel;
import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.modes.SequenceVariationOutputFormat;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.util.WarningCounter;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.PrintWriter;

/**
 * Helper class to implement the logic of discovering sequence variations in and across groups of samples.
 * Implements most of the work done by DiscoverSequenceVariantsMode
 * TODO: estimate reference frequency for indels, or in other words, how many reads overlap the location
 * of an indel and match exactly the reference for the length of the indel.
 *
 * @author Fabien Campagne
 *         Date: Sep 7, 2010
 *         Time: 2:14:38 PM
 * @see edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode
 */
public class DiscoverVariantIterateSortedAlignments extends IterateSortedAlignmentsListImpl {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(DiscoverVariantIterateSortedAlignments.class);

    private int thresholdDistinctReadIndices = 10;
    private int minimumVariationSupport = 3;

    private SequenceVariationOutputFormat format;
    private int numberOfGroups;
    private int[] readerIndexToGroupIndex;

    private GenotypeFilter[] genotypeFilters;
    private int genomeRefIndex;

    public void setMinimumVariationSupport(int minimumVariationSupport) {
        this.minimumVariationSupport = minimumVariationSupport;
    }

    public void setThresholdDistinctReadIndices(int thresholdDistinctReadIndices) {
        this.thresholdDistinctReadIndices = thresholdDistinctReadIndices;
    }

    public DiscoverVariantIterateSortedAlignments(SequenceVariationOutputFormat format) {
        this.format = format;
    }


    Object statWriter;
    String[] samples;

    public void initialize(DiscoverSequenceVariantsMode mode,
                           PrintWriter outWriter,
                           ObjectArrayList<GenotypeFilter> filters) {
        readerIndexToGroupIndex = mode.getReaderIndexToGroupIndex();

        format.defineColumns(outWriter, mode);
        format.setGenome(getGenome());
        genotypeFilters = filters.toArray(new GenotypeFilter[filters.size()]);

    }

    public void finish() {
        format.close();

    }

    private RandomAccessSequenceInterface genome;
    private EquivalentIndelRegionCalculator equivalentIndelRegionCalculator;

    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
        equivalentIndelRegionCalculator = new EquivalentIndelRegionCalculator(genome);
    }

    @Override
    protected void checkGenomeMatchAlignment(final ConcatSortedAlignmentReader sortedReaders,
                                             final RandomAccessSequenceInterface genome) {
        super.checkGenomeMatchAlignment(sortedReaders, genome);
        // permutation was built by previous method, set on eir calculator:
        equivalentIndelRegionCalculator.setReferenceIndexPermutation(alignmentToGenomeTargetIndices);
    }

    @Override
    public RandomAccessSequenceInterface getGenome() {
        return genome;
    }

    public void setOverrideReferenceWithGenome(boolean overrideReferenceWithGenome) {
        this.overrideReferenceWithGenome = overrideReferenceWithGenome;
    }

    public class PositionBaseInfo {
        public int readIndex;
        public int readerIndex;
        public byte qualityScore;
        public boolean matchesReference;
        public char from;
        public char to;
        public int position;
    }

    public void allocateStorage(int numberOfSamples, int numberOfGroups) {
        this.numberOfSamples = numberOfSamples;
        this.numberOfGroups = numberOfGroups;
        sampleCounts = new SampleCountInfo[numberOfSamples];
        for (int i = 0; i < numberOfSamples; i++) {
            sampleCounts[i] = new SampleCountInfo();
        }
        format.allocateStorage(numberOfSamples, numberOfGroups);
        filteredList = new ObjectOpenHashSet<edu.cornell.med.icb.goby.alignments.PositionBaseInfo>();

    }

    private ObjectSet<edu.cornell.med.icb.goby.alignments.PositionBaseInfo> filteredList;
    private SampleCountInfo[] sampleCounts;

    private int numberOfSamples;
    private int previousReference = -1;
    WarningCounter refBaseWarning = new WarningCounter();

    /**
     * Enable or disable calling indels.
     * @param callIndels True will enable calling indels.
     */
    public void setCallIndels(final boolean callIndels) {
        this.callIndels = callIndels;
    }

    private boolean callIndels=false;
    @Override
    public void observeIndel(final Int2ObjectMap<DiscoverVariantPositionData> positionToBases,
                             final int referenceIndex,
                             final int startPosition, final String from, final String to,
                             final int sampleIndex, final int readIndex) {

        if (callIndels) {
            final ObservedIndel indel = new ObservedIndel(startPosition, from, to, readIndex);
            int flankLeftSize = 1;
            equivalentIndelRegionCalculator.setFlankLeftSize(flankLeftSize); // VCF output requires one base before the indel
            equivalentIndelRegionCalculator.setFlankRightSize(0);
            final EquivalentIndelRegion indelCandidateRegion = equivalentIndelRegionCalculator.determine(referenceIndex, indel);
            if (indelCandidateRegion == null) {
                return;
            }

            // subtracts -1 to yield keyPos: the position of the first base in the indel (includes the 1 base left flank, as
            // per VCF spec.) keyPos is zero-based
            final int keyPos = indelCandidateRegion.startPosition - flankLeftSize + 1;
            assert genome.get(referenceIndex, keyPos) == indelCandidateRegion.fromInContext().charAt(0) :
                    "first base of context must match genome at key position";
            indelCandidateRegion.sampleIndex = sampleIndex;

            DiscoverVariantPositionData positionBaseInfos = positionToBases.get(keyPos);
            //   System.out.printf("Observing indel at position %d %n", keyPos);
            if (positionBaseInfos == null) {
                positionBaseInfos = new DiscoverVariantPositionData(keyPos);
                positionToBases.put(keyPos, positionBaseInfos);
            }

            positionBaseInfos.observeCandidateIndel(indelCandidateRegion);
            //printBasesAround(keyPos, positionToBases);
        }

    }

    private void printBasesAround(int keyPos, Int2ObjectMap<DiscoverVariantPositionData> positionBaseInfos) {
        System.out.println("keyPos=" + keyPos);
        for (int i = keyPos - 10; i < keyPos + 10; i++) {
            DiscoverVariantPositionData list = positionBaseInfos.get(i);
            System.out.printf("%c", list != null ? getReferenceAllele(genome, i, list) : ' ');

        }
        System.out.println();
    }

    /**
     * This method is not re-entrant. The behavior is unpredictable if multiple threads call this method. This
     * constraint is added to avoid allocating data structures inside this method, because it is usually called in a loop
     * over all the positions of a genome, and the data structures can be reused, saving garbage collection time.
     *
     * @param referenceIndex Index of the reference sequence where these bases align.
     * @param position       position where the bases are observed
     * @param list           list of bases observed at position.
     */
    @Override
    public void processPositions(final int referenceIndex,
                                 final int position,
                                 final DiscoverVariantPositionData list) {
        int sumVariantCounts = 0;

        if (referenceIndex != previousReference && genome != null) {
            genomeRefIndex = genome.getReferenceIndex(getReferenceId(referenceIndex).toString());
            previousReference = referenceIndex;
        }

        final char referenceBase = getReferenceAllele(genome, position, list);

        for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++) {
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_A_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_T_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_C_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_G_INDEX] = 0;
            sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_OTHER_INDEX] = 0;
            sampleCounts[sampleIndex].referenceBase = referenceBase;
            sampleCounts[sampleIndex].distinctReadIndices.clear();
            sampleCounts[sampleIndex].sampleIndex = sampleIndex;
            sampleCounts[sampleIndex].varCount = 0;
            sampleCounts[sampleIndex].refCount = 0;
            sampleCounts[sampleIndex].failedCount = 0;
            sampleCounts[sampleIndex].clearIndels();
        }

        if (list != null) {
            final IntSet distinctReadIndices = new IntArraySet();

            boolean hasIndel = false;
            if (list.getIndels() != null) {
                for (final EquivalentIndelRegion indel : list.getIndels()) {

                    if (indel.matchesReference()) {
                        sampleCounts[indel.sampleIndex].refCount += indel.getFrequency();

                    } else {
                        sumVariantCounts++;
                        sampleCounts[indel.sampleIndex].varCount += indel.getFrequency();
                        sumVariantCounts += indel.getFrequency();
                    }
                    sampleCounts[indel.sampleIndex].distinctReadIndices.add(indel.readIndex);
                    sampleCounts[indel.sampleIndex].addIndel(indel);
                    hasIndel = true;
                }
            }


            for (final edu.cornell.med.icb.goby.alignments.PositionBaseInfo info : list) {

                if (info.matchesReference && referenceBase != '\0') {
                    // from and to have to be set if the position matches the reference.
                    info.from = referenceBase;
                    info.to = referenceBase;
                }
                final int sampleIndex = info.readerIndex;
                distinctReadIndices.add(info.readIndex);
                if (info.matchesReference) {

                    sampleCounts[sampleIndex].referenceBase = referenceBase;
                    sampleCounts[sampleIndex].refCount++;
                    incrementBaseCounter(info.from, sampleIndex);

                } else {
                    sampleCounts[sampleIndex].varCount++;
                    sumVariantCounts++;
                    if (info.from != referenceBase && info.from != '.' && info.from != '-') {

                        refBaseWarning.warn(LOG, "reference base differ between variation (%c) and genome (%c) at chr %s position %d",
                                info.from, referenceBase, getReferenceId(referenceIndex),
                                position);

                    }
                    sampleCounts[sampleIndex].referenceBase = referenceBase;
                    sampleCounts[sampleIndex].distinctReadIndices.add(info.readIndex);
                    incrementBaseCounter(info.to, sampleIndex);
                }
            }

            if (distinctReadIndices.size() >= thresholdDistinctReadIndices && sumVariantCounts > minimumVariationSupport) {
                final int groupIndexA = 0;
                final int groupIndexB = 1;
                // Do not write statistics for positions in the start flap. The flap start is used to accumulate
                // base counts for reads that can overlap with the window under consideration.

                if (inRegionToWrite(referenceIndex, position)) {

                    if (genotypeFilters.length != 0) {
                        filteredList.clear();
                        for (final GenotypeFilter filter : genotypeFilters) {
                            filter.filterGenotypes(list, sampleCounts, filteredList);
                            //       System.out.printf("filter %s removed %3g %% %n", filter.getName(), filter.getPercentFilteredOut());
                        }
                        final CountFixer fixer = new CountFixer();
                        fixer.fix(list, sampleCounts, filteredList);

                    }

                    // make genotypes comparable across all samples:
                    SampleCountInfo.alignIndels(sampleCounts);

                    format.writeRecord(this, sampleCounts, referenceIndex, position, list, groupIndexA, groupIndexB);

                }

            }
        }
    }

    private boolean inRegionToWrite(final int referenceIndex, final int position) {
        if (useWindow) {
            return !isWithinStartFlap(referenceIndex, position) && !isPastEnd(referenceIndex, position);
        } else {
            return true;
        }
    }

    /**
     * Return true when position is beyond the window length.
     *
     * @param referenceIndex
     * @param position
     * @return
     */
    private boolean isPastEnd(int referenceIndex, int position) {
        return (referenceIndex > endReferenceIndex || referenceIndex == endReferenceIndex && position > endPosition);
    }

    boolean overrideReferenceWithGenome;

    /**
     * Instead of this method, use the random access genome to find the reference base for all bases.
     *
     * @param genome
     * @param position
     * @param list
     * @return
     */
    private char getReferenceAllele(RandomAccessSequenceInterface genome,
                                    int position,
                                    DiscoverVariantPositionData list) {


        if (overrideReferenceWithGenome) {
            char refBase = genome.get(genomeRefIndex, position);
            return refBase;
        } else {
            // We will find some referenceBase among the variations that do not match the reference:
            // this procedure will not be able to determine the refBase if all samples are homzygotes matching the reference
            final ObjectIterator<edu.cornell.med.icb.goby.alignments.PositionBaseInfo> iterator = list.iterator();
            char refBase = '\0';
            // find the reference base from any variant:
            while (iterator.hasNext()) {
                edu.cornell.med.icb.goby.alignments.PositionBaseInfo positionBaseInfo = iterator.next();
                if (!positionBaseInfo.matchesReference) {
                    if (positionBaseInfo.from != '-' && positionBaseInfo.from != '.') {
                        // skip the variant if this was an insertion in the read and we don't know the reference.
                        refBase = positionBaseInfo.from;
                        break;
                    }
                }
            }
            if (refBase == '\0' && genome != null) {
                // look up the reference base since we have a genome:
                refBase = genome.get(genomeRefIndex, position);

            }
            return refBase;
        }

    }

    private void incrementBaseCounter(char base, int sampleIndex) {
        switch (base) {
            case 'A':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_A_INDEX] += 1;
                break;
            case 'T':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_T_INDEX] += 1;
                break;
            case 'C':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_C_INDEX] += 1;
                break;
            case 'G':
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_G_INDEX] += 1;
                break;
            default:
                sampleCounts[sampleIndex].counts[SampleCountInfo.BASE_OTHER_INDEX] += 1;
                break;
        }
    }


}
