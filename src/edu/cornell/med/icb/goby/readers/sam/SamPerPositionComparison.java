/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.readers.sam;

import edu.cornell.med.icb.goby.alignments.Alignments;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Class for comparing sam records. Similar to SamComparison but this class collects all of the records from
 * a single position (assuming the inputs are sorted) and finds the records (based on the clipped read)
 * that match up so they can be compared. If two records within a start position have the exact same clipped
 * read, that position will not have any comparisons.
 */
public class SamPerPositionComparison extends SamComparison {

    private int currentPosition;

    private final List<SAMRecord> sources;
    private final List<SAMRecord> dests;
    private final List<Alignments.AlignmentEntry> gobyDests;
    private final Object2ObjectMap<String, LinkedList<Integer>> destclippedReadToIndexesMap;

    public SamPerPositionComparison() {
        sources = new ArrayList<SAMRecord>();
        dests = new ArrayList<SAMRecord>();
        gobyDests = new ArrayList<Alignments.AlignmentEntry>();
        destclippedReadToIndexesMap = new Object2ObjectOpenHashMap<String, LinkedList<Integer>>();
        currentPosition = -1;
    }

    /**
     * Compare expectedSamRecord vs actualSamRecord. Output details if differences are found.
     * Returns the number of differences found, 0 if the same. Unlike SamComparison, this class
     * compares all the records at a single alignment start position, so until alignment start
     * position changes, this will just emit zero.
     * @param source the expected sam record
     * @param dest the actual sam record
     * @param gobyDest the actual goby alignment record
     * @return the number of comparison failures. 0 is the preferred.
     */
    @Override
    public int compare(final SAMRecord source, final SAMRecord dest, final Alignments.AlignmentEntry gobyDest) {
        int result = 0;
        if (source.getAlignmentStart() != currentPosition) {
            result = makeComparisons();
        }
        final String destClippedRead = usableReadOf(dest);
        LinkedList<Integer> indexesList = destclippedReadToIndexesMap.get(destClippedRead);
        if (indexesList == null) {
            indexesList = new LinkedList<Integer>();
            destclippedReadToIndexesMap.put(destClippedRead, indexesList);
        }
        indexesList.add(dests.size());
        sources.add(source);
        dests.add(dest);
        if (gobyDest != null) {
            gobyDests.add(gobyDest);
        }
        currentPosition = source.getAlignmentStart();
        return result;
    }

    @Override
    public int finished() {
        return makeComparisons();
    }

    private void resetForPosition() {
        sources.clear();
        dests.clear();
        gobyDests.clear();
        destclippedReadToIndexesMap.clear();
    }

    /**
     * Make comparisons of every record at a specific alignment start position.
     * @return The total of the comparison failures for all of the records compared
     */
    private int makeComparisons() {
        int resultSum = 0;
        if (!sources.isEmpty()) {
            for (final SAMRecord source : sources) {
                SAMRecord dest = null;
                Alignments.AlignmentEntry gobyDest = null;

                final String sourceClippedRead = usableReadOf(source);
                final LinkedList<Integer> indexesList = destclippedReadToIndexesMap.get(sourceClippedRead);
                int leastDiffIndex = -1;
                int leastDiffValue = Integer.MAX_VALUE;
                if (indexesList != null) {
                    // Get the index of the first unused matching destination that contains the same
                    // exact usable read as the source. Run a comparison with all diffs that match with
                    // this same read. Keep the one with the fewest differences. If we find a perfect
                    // match, stop and use that one. Unfortunately, without the read name stored in goby,
                    // there isn't a faster way to do this that I can think of.
                    final boolean tempOutputFailedComparisons = outputFailedComparisons;
                    // Save comparisonFailureCount because we may do multiple comparisons
                    // but only one failure or success will count.
                    final int tempComparisonFailureCount = comparisonFailureCount;
                    // Never output failed comparisons during this phase
                    outputFailedComparisons = false;
                    for (final int currentDestIndex : indexesList) {
                        dest = dests.get(currentDestIndex);
                        if (!gobyDests.isEmpty()) {
                            gobyDest = gobyDests.get(currentDestIndex);
                        }
                        final int currentDiffValue = super.compare(source, dest, gobyDest);
                        if (currentDiffValue < leastDiffValue) {
                            leastDiffValue = currentDiffValue;
                            leastDiffIndex = currentDestIndex;
                            if (currentDiffValue == 0) {
                                // Stop if we get a diff of 0
                                break;
                            }
                        }
                    }
                    outputFailedComparisons = tempOutputFailedComparisons;
                    comparisonFailureCount = tempComparisonFailureCount;
                    if (leastDiffIndex != -1) {
                        // Found the best destination record
                        indexesList.removeFirstOccurrence(leastDiffIndex);
                        if (leastDiffValue > 0) {
                            // To have the correct comparison failure details, we need to re-compare or the
                            // dump will have an incorrect comparison failure list.
                            super.compare(source, dest, gobyDest);
                            if (outputFailedComparisons) {
                                // Output the failed but still best failure here.
                                dest = dests.get(leastDiffIndex);
                                if (!gobyDests.isEmpty()) {
                                    gobyDest = gobyDests.get(leastDiffIndex);
                                }
                                dumpComparison(source, dest, gobyDest);
                            }
                        }
                    }
                }
                if (leastDiffIndex == -1) {
                    // Couldn't find anything to diff with AT ALL, nothing with the same clipped read bases
                    readNum++;
                    comparisonFailureCount++;
                    System.out.println("ERROR: Couldn't find a dest record to compare against for readName=" + source.getReadName());
                    // Consider: When this is done, pick the sources that didn't get compared with the dest's that
                    // didn't get compared. Often there may just be 1 and 1??
                } else {
                    resultSum += leastDiffValue;
                }
            }
        }
        resetForPosition();
        return resultSum;
    }
}
