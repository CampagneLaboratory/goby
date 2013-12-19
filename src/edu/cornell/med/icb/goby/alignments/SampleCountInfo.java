/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.EquivalentIndelRegion;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.lang.MutableString;

import java.util.Collections;
import java.util.Comparator;

/**
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 11:37:42 AM
 */
public class SampleCountInfo {
    public static final int BASE_A_INDEX = 0;
    public static final int BASE_T_INDEX = 1;
    public static final int BASE_C_INDEX = 2;
    public static final int BASE_G_INDEX = 3;
    public static final int BASE_OTHER_INDEX = 4;
    public static final int BASE_MAX_INDEX = BASE_OTHER_INDEX + 1;

    public char referenceBase;
    public IntSet distinctReadIndices = new IntArraySet();
    public int sampleIndex;
    public int varCount;
    public int refCount;


    /**
     * Number of bases that failed base filters, for any reason.
     */

    public int failedCount;
    private int[][] counts = new int[2][5];

    public SampleCountInfo() {
        this.counts[0] = new int[5];
        this.counts[1] = new int[5];
    }

    /**
     * While base genotypes were flagged for removal by some filter in this sample.
     */
    public boolean filtered[] = new boolean[5];

    /**
     * List of indel eirs that start at the position in this sample.
     */
    private ObjectArrayList<EquivalentIndelRegion> indels;
    private static final Comparator<? super EquivalentIndelRegion> INDEL_COMPARATOR = new Comparator<EquivalentIndelRegion>() {
        /**
         * Simply order eir by to string, alphabetically.
         * @param eir1 first indel
         * @param eir2 second indel
         * @return sort order
         */
        @Override
        public int compare(final EquivalentIndelRegion eir1, final EquivalentIndelRegion eir2) {
            return eir1.to.compareTo(eir2.to);
        }
    };
    private int sumOfCounts;


    /**
     * Add an indel to this sample info. This method lazily creates the indels collection when the first
     * indel is added. If the same indel was already added (equality by range of eir start-end), the frequency
     * is incremented.
     *
     * @param indel eir observed starting at the position this sampleCountInfo is associated with.
     */
    public void addIndel(final EquivalentIndelRegion indel) {
        if (!indel.isFiltered()) {
            // only add indels that were not filtered.
            if (indels == null) {
                indels = new ObjectArrayList<EquivalentIndelRegion>();
            }
            int previousStartPosition = -1;
            for (final EquivalentIndelRegion prevIndel : indels) {
                if (prevIndel.equals(indel)) {
                    prevIndel.incrementFrequency();
                }
                previousStartPosition = prevIndel.startPosition;

            }
            if (previousStartPosition != -1) {
                assert indel.startPosition == previousStartPosition : "You can only add indels with the same start position in any given SampleCountInfo ";
            }

            indels.add(indel);
        }
    }


    /**
     * Return the maximum index for genotypes observed in this sample.
     *
     * @return the maximum index that will retrieve a genotype count.
     */
    public final int getGenotypeMaxIndex() {
        return BASE_MAX_INDEX + (hasIndels() ? indels.size() : 0);
    }


    /**
     * Align indel genotypes so that genotype indices refer to the same genotype across a set of samples.
     *
     * @param samples a set of samples under comparison.
     */
    public static void alignIndels(final SampleCountInfo[] samples) {
        final ObjectArraySet<EquivalentIndelRegion> dummyIndels = new ObjectArraySet<EquivalentIndelRegion>();
        boolean hasSomeIndels = false;
        for (final SampleCountInfo sample : samples) {
            if (sample.hasIndels()) {
                hasSomeIndels = true;
                for (final EquivalentIndelRegion indel : sample.indels) {
                    if (!dummyIndels.contains(indel)) {

                        dummyIndels.add(indel.copy());
                    }
                }
            }
        }
        if (!hasSomeIndels) {
            // no indels, done aligning.
            return;
        }

        // dummyIndels contains the canonical set of indels seen across all the samples.
        for (final SampleCountInfo sample : samples) {
            for (final EquivalentIndelRegion indel : dummyIndels) {
                if (!sample.hasMatching(indel)) {
                    // make a copy again for this specific sample, set the sampleIndex:
                    final EquivalentIndelRegion copy = indel.copy();
                    copy.sampleIndex = sample.sampleIndex;
                    //  we just create a dummy observation to align genotypes across samples and set its filtered state so that the frequency is always zero

                    copy.markFiltered();
                    if (sample.indels == null) {
                        sample.indels = new ObjectArrayList<EquivalentIndelRegion>();
                    }
                    sample.indels.add(copy);
                }
            }
            sample.sortGenotypes();
        }
    }


    /**
     * Sort genotypes in a predefined order. Indels are sorted by increasing to field (alphabetically).
     */
    private void sortGenotypes() {
        Collections.sort(indels, INDEL_COMPARATOR);
    }

    /**
     * Determine if the sample has an indel matching the argument.
     *
     * @param indel argument to look for.
     * @return True if the indels collection for this sample has an indel at the same position, with the same to sequence. False otherwise.
     */
    private boolean hasMatching(final EquivalentIndelRegion indel) {
        if (!hasIndels()) {
            return false;
        }
        for (final EquivalentIndelRegion eir : indels) {
            if (eir.startPosition == indel.startPosition && eir.endPosition == indel.endPosition &&
                    eir.referenceIndex == indel.referenceIndex &&
                    eir.to.equals(indel.to)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Return true if this sample has indel observations (at least one indel with frequency>0).
     *
     * @return True or false.
     */
    public boolean hasIndels() {
        return !(indels == null || indels.isEmpty()/*|| allMarkedRemoved() */);
    }
    /**
     * Return true if this sample has any indel (even with frequency=0).
     *
     * @return True or false.

    private boolean hasSomeIndels() {
    return !(indels == null || indels.isEmpty() );
    }
     */
    /**
     * Determine if the indels are all marked removed.
     *
     * @return True when all indels were marked removed. False if any indel exists such that markRemoved=false.
     */
    private boolean allMarkedRemoved() {
        boolean allRemoved = true;
        for (EquivalentIndelRegion eir : indels) {
            allRemoved &= eir.isFiltered();
        }
        return allRemoved;
    }

    public final char base(final int baseIndex) {
        switch (baseIndex) {
            case BASE_A_INDEX:
                return 'A';
            case BASE_C_INDEX:
                return 'C';
            case BASE_T_INDEX:
                return 'T';
            case BASE_G_INDEX:
                return 'G';
            default:
                return 'N';
        }
    }

    public final int baseIndex(final char to) {
        switch (to) {
            case 'A':
                return BASE_A_INDEX;
            case 'C':
                return BASE_C_INDEX;
            case 'T':
                return BASE_T_INDEX;
            case 'G':
                return BASE_G_INDEX;
            default:
                return BASE_OTHER_INDEX;
        }
    }

    public ObjectArrayList<EquivalentIndelRegion> getEquivalentIndelRegions() {
        return indels;
    }

    /**
     * Return the indel at a given genotype index.
     *
     * @param genotypeIndex index of the indel genotype
     * @return indel EIR.
     */
    public EquivalentIndelRegion getIndelGenotype(int genotypeIndex) {

        assert genotypeIndex >= BASE_MAX_INDEX : "genotype index must not refer to a base when fetching an indel.";
        assert hasIndels() : "SampleCountInfo must have indels to fetch an indel";

        final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
        return indels.get(indelIndex);
    }

    public final void getGenotype(final int genotypeIndex, final MutableString destination) {
        destination.setLength(0);
        if (genotypeIndex < BASE_MAX_INDEX) {
            destination.append(base(genotypeIndex));
            return;
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                destination.append(indels.get(indelIndex).to);
                return;
            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);
    }


    public final String getGenotypeString(final int genotypeIndex) {
        if (genotypeIndex < BASE_MAX_INDEX) {

            return STRING[genotypeIndex];
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                return indels.get(indelIndex).toInContext();

            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);
    }


    public boolean isReferenceGenotype(final int genotypeIndex) {
        if (genotypeIndex < BASE_MAX_INDEX) {

            return base(genotypeIndex) == referenceBase;
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                final EquivalentIndelRegion equivalentIndelRegion = indels.get(indelIndex);
                return equivalentIndelRegion.to.equals(equivalentIndelRegion.from);
            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);
    }

    public boolean isIndel(int genotypeIndex) {
        if (genotypeIndex < BASE_MAX_INDEX) {

            return false;
        } else {
            return hasIndels() ? genotypeIndex - BASE_MAX_INDEX < indels.size() : false;
        }

    }

    public String getReferenceGenotype() {
        if (hasIndels()) {
            return indels.get(0).fromInContext();
        } else {
            return Character.toString(referenceBase);

        }
    }

    /**
     * Remove all previously recorded indels.
     */
    public void clearIndels() {
        if (indels != null) {
            indels.clear();
        }

    }

    static final String A_BASE = "A";
    static final String T_BASE = "T";
    static final String C_BASE = "C";
    static final String G_BASE = "G";
    static final String N_BASE = "N";
    static final String[] STRING = {A_BASE, T_BASE, C_BASE, G_BASE, N_BASE};

    public final String baseString(final int genotypeIndex) {
        if (genotypeIndex < BASE_MAX_INDEX) {
            return STRING[genotypeIndex];
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                return indels.get(indelIndex).toString();
            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);

    }

    /**
     * Remove an indel from the list of candidate indels.
     *
     * @param indel to remove.
     */
    public void removeIndel(final EquivalentIndelRegion indel) {

        for (final EquivalentIndelRegion eir : indels) {
            if (!eir.isFiltered() && indel.equals(eir)) {
                eir.markFiltered();
            }
        }
    }

    @Override
    public String toString() {
        return String.format("sample: %d counts %s %s %s %s %s FB=%d indels={ %s }%n",
                sampleIndex,
                toString(BASE_A_INDEX),
                toString(BASE_T_INDEX),
                toString(BASE_C_INDEX),
                toString(BASE_G_INDEX),
                toString(BASE_OTHER_INDEX),
                failedCount,
                indels);
    }

    private String toString(int baseIndex) {
        StringBuffer buffer = new StringBuffer();

        buffer.append(base(baseIndex));
        buffer.append("=");
        buffer.append(getGenotypeCount(baseIndex));

        return buffer.toString();

    }


    public void suggestRemovingGenotype(int baseIndex, boolean matchesForwardStrand) {

        if (getGenotypeCount(baseIndex, matchesForwardStrand) - 1 >= 0) {

            decrementGenotypeCount(baseIndex, matchesForwardStrand);
        }

    }


    /**
     * Return the genotype frequency in the sample. The number of times the genotype was observed in the sample.
     *
     * @param genotypeIndex Index of the genotype for which count/frequency is sought.
     * @return the frequency.
     */
    public final int getGenotypeCount(final int genotypeIndex) {

        if (genotypeIndex < BASE_MAX_INDEX) {
            return Math.max(0, counts[POSITIVE_STRAND][genotypeIndex]) + Math.max(0, counts[NEGATIVE_STRAND][genotypeIndex]);
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                return indels.get(indelIndex).getFrequency();
            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);
    }

    /**
     * Return the genotype frequency in the sample. The number of times the genotype was observed in the sample.
     *
     * @param genotypeIndex Index of the genotype for which count/frequency is sought.
     * @return the frequency.
     */
    public final int getGenotypeCount(final int genotypeIndex, boolean matchesForwardStrand) {

        if (genotypeIndex < BASE_MAX_INDEX) {
            return Math.max(0, counts[recodeStrand(matchesForwardStrand)][genotypeIndex]);
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                return indels.get(indelIndex).getFrequency();
            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);
    }

    private int recodeStrand(boolean matchesForwardStrand) {
        return matchesForwardStrand ? POSITIVE_STRAND : NEGATIVE_STRAND;
    }

    /**
     * Reset a genotype count to a given value.
     *
     * @param genotypeIndex
     * @param count
     */
    public void setGenotypeCount(int genotypeIndex, int count, boolean matchesForwardStrand) {
        if (genotypeIndex < BASE_MAX_INDEX) {
            counts[recodeStrand(matchesForwardStrand)][genotypeIndex] = count;
            return;
        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                EquivalentIndelRegion equivalentIndelRegion = indels.get(indelIndex);
                equivalentIndelRegion.setFrequency(count);
                equivalentIndelRegion.removeFiltered();
                return;
            }
        }
        throw new IllegalArgumentException("The genotype index argument was out of range: " + genotypeIndex);

    }

    /**
     * Reset a genotype count to a given value.
     *
     * @param genotypeIndex
     * @param count
     */
    public void setGenotypeCount(int genotypeIndex, int count) {
        setGenotypeCount(genotypeIndex, count, true);
    }

    private int decrementGenotypeCount(int baseIndex, boolean matchesForwardStrand) {
        return --(counts[recodeStrand(matchesForwardStrand)][baseIndex]);
    }

    public int incrementGenotypeCount(int baseIndex, boolean matchesForwardStrand) {
        return (++counts[recodeStrand(matchesForwardStrand)][baseIndex]);
    }

    private static final int POSITIVE_STRAND = 1;
    private static final int NEGATIVE_STRAND = 0;


    public int getSumCounts() {
        int sum = 0;
        for (int i = 0; i < getGenotypeMaxIndex(); i++) {
            sum += getGenotypeCount(i);
        }
        return sum;
    }

    public boolean isFiltered(int genotypeIndex) {
        if (genotypeIndex < BASE_MAX_INDEX) {
            return filtered[genotypeIndex];

        } else {
            if (hasIndels()) {
                final int indelIndex = genotypeIndex - BASE_MAX_INDEX;
                EquivalentIndelRegion equivalentIndelRegion = indels.get(indelIndex);

                return equivalentIndelRegion.isFiltered();
            } else return false;
        }
    }

    /**
     * Return the base frequency of this genotype in this one sample.
     *
     * @param genotypeIndex Index of a genotype
     * @return the frequency of this genotype (count of the genotype divided by the sum of all other counts).
     */
    public float frequency(int genotypeIndex) {
        int sum = 0;

        for (int index = 0; index < getGenotypeMaxIndex(); index++) {
            sum += getGenotypeCount(index);
        }
        if (sum == 0) return 0;
        float count = (float) getGenotypeCount(genotypeIndex);
        return count / ((float) sum);
    }


    public void clearGenotypeCount(int baseIndex) {
        counts[0][baseIndex] = 0;
        counts[1][baseIndex] = 0;
    }

    /**
     * return the number of base genotypes.
     *
     * @return
     */
    public int getCountsSize() {
        return counts.length;
    }

    public boolean anyCountNegative() {
        for (int i = 0; i < getCountsSize(); i++) {
            if (getGenotypeCount(i) < 0) return true;
        }
        return false;
    }

    public int getSumOfCounts() {
        int sum = 0;
        for (int count : counts[0]) {
            sum += count;
        }
        for (int count : counts[1]) {
            sum += count;
        }
        return sum;
    }

    /**
     * Return the number of bases that cover this site, irrespective of genotype.
     *
     * @return the coverage at the site.
     */
    public int coverage() {
        int coverage = 0;
        for (int genotypeIndex = 0; genotypeIndex < getGenotypeMaxIndex(); genotypeIndex++) {
            coverage += getGenotypeCount(genotypeIndex);

        }
        return coverage;
    }

    /**
     * Return the number of bases that cover this site, irrespective of genotype.
     *
     * @return the coverage at the site.
     */
    public int coverage() {
       int coverage=0;
        for (int genotypeIndex=0; genotypeIndex<getGenotypeMaxIndex(); genotypeIndex++){
            coverage+=getGenotypeCount(genotypeIndex);

        }
        return coverage;
    }
}
