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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import edu.cornell.med.icb.goby.R.GobyRengine;
import org.rosuda.JRI.Rengine;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.collections.map.LRUMap;

/**
 * Removes likely sequencing errors from the list of variantions. Either filter the errors completely or tries to determine
 * which base should have been sequenced.
 *
 * @author Fabien Campagne
 *         Date: Mar 21, 2011
 *         Time: 11:39:40 AM
 */
public class FisherBaseFilter extends BaseFilter {

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(FisherBaseFilter.class);

    private boolean fisherRInstalled;
    private boolean firstReport = true;
    private ObjectArrayList<ReadIndexStats> readIndexStats;
    private LRUMap fisherPCache;
    private Double pValueThreshold = 0.05;

    @Override
    public String describe() {
        return "fisher p>" + pValueThreshold;
    }

    public void setPValueThreshold(Double pValueThreshold) {
        this.pValueThreshold = pValueThreshold;
    }

    public FisherBaseFilter(ObjectArrayList<ReadIndexStats> readIndexStats) {
        fisherPCache = new LRUMap(10000);

        this.readIndexStats = readIndexStats;
        //activate R only if we need it:
        final Rengine rEngine = GobyRengine.getInstance().getRengine();
        fisherRInstalled = rEngine != null && rEngine.isAlive();

    }

    public void filterBases(ObjectArrayList<PositionBaseInfo> list,
                            SampleCountInfo[] sampleCounts,
                            ObjectArrayList<PositionBaseInfo> filteredList) {


        if (!fisherRInstalled) {
            if (firstReport) {
                System.err.println("R connection is needed to detect errors in sequence variations and adjust counts.");
                firstReport = false;
            }
            return;
        }

        for (SampleCountInfo sci : sampleCounts) {
            final ReadIndexStats stats = readIndexStats.get(sci.sampleIndex);
            ObjectArrayList<PositionBaseInfo> considered = new ObjectArrayList<PositionBaseInfo>();
            for (int baseIndex = SampleCountInfo.BASE_A_INDEX;
                 baseIndex < SampleCountInfo.BASE_MAX_INDEX;
                 baseIndex++) {
                considered.clear();
                int observedReferenceCount = 0;
                double expectedVariationRate = 0;
                int expectedVariationCount = 0;
                int expectedReferenceCount = 0;
                int observedVariationCount = 0;
                int totalCount = 0;
                for (int count : sci.counts) {
                    totalCount += count;
                }

                observedVariationCount = sci.counts[baseIndex];
                if (observedVariationCount == 0) {
                    // no filtering needed for this base.
                    continue;
                }
                // determine if the count for this allele is the likely result of sequencing errors.

                long sum = 0;

                for (PositionBaseInfo info : list) {
                    if (info.readerIndex != sci.sampleIndex) continue;

                    final int readIndex = info.readIndex;
                    if (readIndex < 1 ||
                            readIndex > stats.countVariationBases.length ||
                            readIndex > stats.countReferenceBases.length) {

                        // TODO this test is a kludge: readIndex should always be within the bounds of countVariationBases or countReferenceBases
                        // the fact that we need this test indicates that there is a bug in the calculation of readIndex, probably when read insertions or deletions are present.
                        //              assert false : "should not fail" + String.format("readIndex =%d || readIndex >= stats.countVariationBases.length || readIndex >= stats.countReferenceBases.length",
                        //                     readIndex);
                        continue;
                    }
                    long variationBases = stats.countVariationBases[readIndex - 1];
                    long referenceBases = stats.countReferenceBases[readIndex - 1];

                    //  System.out.printf("readIndex=%d variationBases=%d referenceBases=%d %n",readIndexInfo.readIndex, variationBases, referenceBases);
                    expectedVariationRate += variationBases;
                    sum += variationBases + referenceBases;
                    //   refBaseCountAnyPosition += referenceBases;
                    if (info.to == sci.base(baseIndex)) {
                        considered.add(info);
                    }
                }
                expectedVariationRate /= sum;
                //   System.out.printf("Expected variation rate: %f%n", expectedVariationRate);
                expectedVariationCount = (int) Math.round(expectedVariationRate * (double) (totalCount));
                expectedReferenceCount = (int) Math.round((1 - expectedVariationRate) * (double) (totalCount));

                // if the allele is observed with more counts than expected from sequencing errors, we keep this allele count.
                int count00 = expectedVariationCount;
                int count10 = observedVariationCount;
                int count01 = totalCount - count00;
                int count11 = totalCount - count01;
                Double pValue;
                {
                    contingencyValue value = new contingencyValue(count00, count10, count01, count11);
                    pValue = (Double) fisherPCache.get(value);
                    if (pValue == null) {
                        // if (LOG.isTraceEnabled()) {
                        pValue = estimatePValue(count00, count10, count01, count11);

                        fisherPCache.put(value, pValue);
                    }
                }

                if (pValue > pValueThreshold) {
                    int numErroneouslyCalledBases = sci.counts[baseIndex];

                    //   System.out.printf("filtering out %d counts %n", numErroneouslyCalledBases);
                    /*  System.out.printf("p=%f filtering out <%d> %d %d %d %d %n", pValue,  numErroneouslyCalledBases,
                            expectedVariationCount, observedVariationCount,
                            expectedReferenceCount, observedReferenceCount);
                    */
                    sci.counts[baseIndex] = 0;
                    filteredList.addAll(considered);

                }
            }
        }
        numScreened += list.size();
        numFiltered += filteredList.size();
    }


    public final Double estimatePValue(int count00, int count10, int count01, int count11) {
        Double pValue;

        pValue = fisherRInstalled ? FisherExactRCalculator.getFisherOneTailedLesserPValue(
                count00, count10,
                count01, count11
        ) : Double.NaN;
        /* System.out.println(String.format("contingency : %n" +
                "    [  exp    obs ]%n" +
                "ref [ %d       %d ]%n" +
                "var [ %d       %d ]%n" +
                 "P=%f",
                count00, count10, count01, count11,pValue));
                */
        return pValue;
    }


    private class contingencyValue {
        int count00;
        int count10;
        int count01;

        contingencyValue(int count00, int count10, int count01, int count11) {
            this.count00 = count00;
            this.count10 = count10;
            this.count01 = count01;
            this.count11 = count11;
        }

        int count11;
    }


}
