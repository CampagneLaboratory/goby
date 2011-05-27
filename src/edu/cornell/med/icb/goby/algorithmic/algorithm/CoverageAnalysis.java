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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import edu.cornell.med.icb.goby.counts.AnyTransitionCountsIterator;
import edu.cornell.med.icb.goby.counts.CountsAggregatorI;
import edu.cornell.med.icb.goby.counts.CountsReaderI;
import edu.cornell.med.icb.goby.counts.UnionDumpIterator;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.io.IOException;

/**
 * A class to calculate coverage statistics from Counts (histogram base level) data.
 *
 * @author Fabien Campagne
 *         Date: May 21, 2011
 *         Time: 10:44:00 AM
 */
public class CoverageAnalysis {
    CountsAggregatorI orIterator;
    /**
     * Arrays where each element is the number of bases observed for which exactly i reads match span the base. The index of the array is i.
     */
    LongArrayList depthTallyBasesInAnnotation = new LongArrayList(10000);
    LongArrayList depthTallyBasesOutsideAnnotation = new LongArrayList(10000);
    /**
     * Arrays where each element is the number of genomic position/site observed for which exactly i reads match span the base. The index of the array is i.
     */
    LongArrayList depthTallySitesInAnnotation = new LongArrayList(10000);
    LongArrayList depthTallySitesOutsideAnnotation = new LongArrayList(10000);

    private long sumDepth = 0;
    private long countDepth = 0;
    private long sumDepthAnnot = 0;
    /**
     * The number of annotated/captured sites:
     */
    private long countDepthAnnot = 0;
    private boolean statsEstimated = false;
    private long[] cumulativeSitesNotCaptured;
    private long[] cumulativeSitesCaptured;


    public LongArrayList getDepthTallyBasesInAnnotation() {
        return depthTallyBasesInAnnotation;
    }

    public LongArrayList getDepthTallyBasesOutsideAnnotation() {
        return depthTallyBasesOutsideAnnotation;
    }

    private long countAllBases = 0;

    public long getSumDepthAnnot() {
        return sumDepthAnnot;
    }

    public long getCountDepthAnnot() {
        return countDepthAnnot;
    }

    public long getCountAllBases() {
        return countAllBases;
    }

    public CoverageAnalysis() throws IOException {

    }

    private void grow(LongArrayList update, int depth) {
        for (int k = update.size(); k <= depth; k++) {
            update.add(0);
        }
    }

    public void process(CountsReaderI annotationReader, CountsReaderI reader) throws IOException {
        orIterator = new UnionDumpIterator(reader, annotationReader);

        while (orIterator.hasNextTransition()) {
            orIterator.nextTransition();
            int readerCount = orIterator.getCount(0);
            int annotationCount = orIterator.getCount(1);
            int position = orIterator.getPosition();
            int length = orIterator.getLength();
            int end = position + length;

            boolean inAnnotation = annotationCount == 1;
            LongArrayList updateBases = inAnnotation ? depthTallyBasesInAnnotation : depthTallyBasesOutsideAnnotation;
            LongArrayList updateSites = inAnnotation ? depthTallySitesInAnnotation : depthTallySitesOutsideAnnotation;
            int depth = readerCount;
            final int numBases = depth * length;
            if (depth != 0) {
               // System.out.printf("Position %d Adding %d to sumDepth.%n",position, numBases);
              //  System.out.printf("Position %d Adding %d to countDepth.%n",position, length);
                sumDepth += numBases;
                countDepth += length;
                if (inAnnotation) {
                    sumDepthAnnot += numBases;
                    countDepthAnnot += length;

                }
            }
            grow(depthTallyBasesInAnnotation, depth);
            grow(depthTallyBasesOutsideAnnotation, depth);
            grow(depthTallySitesInAnnotation, depth);
            grow(depthTallySitesOutsideAnnotation, depth);
            // count bases over constant count segment: depth time length

            updateSites.set(depth, updateSites.get(depth) + length);
            updateBases.set(depth, updateBases.get(depth) + numBases);
            /*     if (false) System.out.printf("position=%d depth=%d #sites=%d length=%d %n",
                       position, depth,
                       depthTallySitesInAnnotation.get(depth) + depthTallySitesOutsideAnnotation.get(depth),
                       length);
           if (false)    System.out.printf("annotations only position=%d depth=%d #sites=%d length=%d %n",
                       position, depth,
                       depthTallySitesInAnnotation.get(depth),
                       length);
            */
            countAllBases += depth * length;

        }
        orIterator.close();
        orIterator = null;
        Runtime.getRuntime().gc();
    }

    private double sum(long[] array) {
        double sum = 0;
        int o = 0;
        for (long value : array) {

            sum += value;

        }
        return sum;
    }

    private double divide(double v1, double v2) {
        return v1 / v2;
    }

    private double sum(LongArrayList array, int offset) {
        double sum = 0;
        int o = 0;
        for (long value : array) {
            if (o >= offset) {
                sum += value;
            }
            ++o;
        }
        return sum;
    }

    private long[] cumulativeBasesCaptured;
    private long[] cumulativeBasesNotCaptured;

    public void estimateStatistics() {

        final int length = getNumberOfDepths();
        cumulativeBasesCaptured = new long[length];
        cumulativeBasesNotCaptured = new long[length];
        cumulativeSitesCaptured = new long[length];
        cumulativeSitesNotCaptured = new long[length];
        long sumBasesCaptured = (long) sum(depthTallyBasesInAnnotation, 0);
        long sumBasesNotCaptured = (long) sum(depthTallyBasesOutsideAnnotation, 0);
        long sumSitesCaptured = (long) sum(depthTallySitesInAnnotation, 0);
        long sumSitesNotCaptured = (long) sum(depthTallySitesOutsideAnnotation, 0);

        for (int depth = 0; depth < length; ++depth) {

            cumulativeBasesCaptured[depth] = sumBasesCaptured;
            sumBasesCaptured -= depthTallyBasesInAnnotation.get(depth);
            cumulativeBasesNotCaptured[depth] = sumBasesNotCaptured;
            sumBasesNotCaptured -= depthTallyBasesOutsideAnnotation.get(depth);


            cumulativeSitesCaptured[depth] = sumSitesCaptured;
            sumSitesCaptured -= depthTallySitesInAnnotation.get(depth);
            cumulativeSitesNotCaptured[depth] = sumSitesNotCaptured;
            sumSitesNotCaptured -= depthTallySitesOutsideAnnotation.get(depth);
        }
        statsEstimated = true;
    }

    /**
     * Return enrichment efficiency (ee). We define enrichment efficiency as the number of bases captured divided
     * by the total number of bases observed and mapped in the sequencing experiment.
     * Formally, ee=A/(A+O), where A is the number of bases mapped to the reference whithin an annotation
     * segment, and O is the number of bases mapped to the reference that do not overlap an annotation.
     * Enrichment efficiency is larger or equal to zero and less or equal to one.
     *
     * @return enrichment efficiency.
     */
    public double getEnrichmentEfficiency() {
        checkStatsEstimated();

        return divide(cumulativeBasesCaptured[1], cumulativeBasesCaptured[1] + cumulativeBasesNotCaptured[1]);

    }

    private void checkStatsEstimated() {
        if (!statsEstimated)
            throw new IllegalStateException("Statistics must be estimated before they can be accesed.");
    }

    public long getSumDepth() {
        return sumDepth;
    }

    public long getCountDepth() {
        return countDepth;
    }

    /**
     * The number of distinct base depths observed.
     *
     * @return number of depths.
     */
    public int getNumberOfDepths() {
        return depthTallyBasesInAnnotation.size();
    }

    /**
     * Returns the average depth, over annotations or non annotation sites. Only sites that have at least one
     * base mapped are considered.
     *
     * @return
     */
    public double getAverageDepth() {
        return divide(sumDepth, countDepth);
    }

    /**
     * Returns the average depth estimated exclusively over annotations. Only sites that have at least one
     * base mapped over are considered.
     *
     * @return Average depth over annotations.
     */
    public double getAnnotationAverageDepth() {
        return divide(sumDepthAnnot, countDepthAnnot);
    }

    /**
     * Returns the average depth estimated exclusively over nases outside of annotations.
     *
     * @return Average depth outside of annotations.
     */
    public double getNotAnnotationAverageDepth() {
        return divide(sumDepth - sumDepthAnnot, countDepth - countDepthAnnot);
    }

    /**
     * Return the total number of bases observed at sites with depth equal or larger to d
     *
     * @param d depth threshold
     * @return total number of bases
     */
    public long getNumBasesWithDepthAtLeast(int d) {
        checkStatsEstimated();

        if (d >= cumulativeBasesCaptured.length) return 0;
        return cumulativeBasesCaptured[d] + cumulativeBasesNotCaptured[d];
    }

    /**
     * Return the total number of sites with depth equal or larger to d
     *
     * @param d depth threshold
     * @return total number of sites with depth>=d
     */
    public double getNumSitesWithDepthAtLeast(int d) {
        checkStatsEstimated();

        if (d >= cumulativeSitesCaptured.length) return 0;
        return cumulativeSitesCaptured[d] + cumulativeSitesNotCaptured[d];
    }

    /**
     * Return the total number of captured sites with depth equal or larger to d
     *
     * @param d depth threshold
     * @return total number of captured sites with depth>=d
     */
    public double getNumSitesCapturedWithDepthAtLeast(int d) {
        checkStatsEstimated();

        if (d >= cumulativeSitesCaptured.length) return 0;
        return cumulativeSitesCaptured[d];
    }

    /**
     * Return the total number of non captured sites with depth equal or larger to d
     *
     * @param d depth threshold
     * @return total number of non-captured sites with depth>=d
     */
    public double getNumSitesNotCapturedWithDepthAtLeast(int d) {
        checkStatsEstimated();

        if (d >= cumulativeSitesNotCaptured.length) return 0;
        return cumulativeSitesNotCaptured[d];
    }

    /**
     * Return the depth observed for p percentile sites when sites are ordered by depth. If p=.9, returns
     * the d such that 90% of the sites have depth greater or equal to d.
     *
     * @param percentile
     * @return
     */
    public int depthAtPercentile(double percentile) {
        checkStatsEstimated();
        double numSites = countDepth;
        int size = cumulativeSitesCaptured.length;
        final double threshold = numSites * percentile;
        for (int depth = 0; depth < size; ++depth) {

            long sum = cumulativeSitesCaptured[depth] + cumulativeSitesNotCaptured[depth];

            if (sum <= threshold) return depth;
        }
        return size - 1;
    }

    /**
     * Return the depth observed for p percentile captured sites when sites are ordered by depth. If p=.9, returns
     * the d such that 90% of the captured sites have depth greater or equal to d.
     *
     * @param percentile a number between 0 and 1.
     * @return
     */
    public int depthCapturedAtPercentile(double percentile) {
        checkStatsEstimated();
        double numSites = countDepthAnnot;
        int size = cumulativeSitesCaptured.length;
        final double threshold = numSites * percentile;
        for (int depth = 0; depth < size; ++depth) {

            long sum = cumulativeSitesCaptured[depth];
            //   System.out.printf("depth=%d sum %d threshold %g%n",depth, sum,threshold);
            if (sum <= threshold) return depth;
        }
        return size - 1;
    }

    public LongArrayList getDepthTallySitesInAnnotation() {
        return depthTallySitesInAnnotation;
    }


    public long[] getCumulativeBasesCaptured() {
        return cumulativeBasesCaptured;
    }

    public long[] getCumulativeSitesCaptured() {
        return cumulativeSitesCaptured;
    }

    public long[] getCumulativeSitesTotal() {
        long[] result = new long[cumulativeBasesCaptured.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = cumulativeSitesCaptured[i] + cumulativeSitesNotCaptured[i];
        }
        return result;
    }

    public long[] getCumulativeSitesNotCaptured() {
        return cumulativeSitesNotCaptured;
    }

    /**
     * Returns the percentage of captured sites whose depth d is larger or equal to d
     *
     * @param d
     * @return
     */
    public double percentSitesCaptured(int d) {
        return divide(cumulativeSitesCaptured[d], countDepthAnnot);
    }
}
