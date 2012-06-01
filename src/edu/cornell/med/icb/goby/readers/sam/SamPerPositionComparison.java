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
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Class for comparing sam records. Similar to SamComparison but this class collects all of the records from
 * a single position (assuming the inputs are sorted) and finds the records (based on the clipped read)
 * that match up so they can be compared. If two records within a start position have the exact same clipped
 * read, that position will not have any comparisons.
 */
public class SamPerPositionComparison extends SamComparison {

    private static final List<SAMRecord> EMPTY_SAM_LIST = new ArrayList<SAMRecord>(0);
    private static final List<Alignments.AlignmentEntry> EMPTY_GOBY_LIST = new ArrayList<Alignments.AlignmentEntry>(0);

    private int currentPosition;

    private final List<SAMRecord> sources;
    private final Map<String, List<SAMRecord>> destsMap;
    private final Map<String, List<Alignments.AlignmentEntry>> gobyDestsMap;

    /**
     * The default (true) means when a new position is found in the source record, comparisons will be made
     * with all the records that were found before.
     * If false, comparisons will be made after ALL records have been read. This will take a lot more memory
     * (enough to load two or three whole alignments in memory) as ALL records will need to be read into memory)
     * BUT, this will be very tolerant of the sorting of the input files.
     */
    private boolean compareAtEachNewPosition = true;

    public SamPerPositionComparison() {
        sources = new ArrayList<SAMRecord>();
        destsMap = new HashMap<String, List<SAMRecord>>();
        gobyDestsMap = new HashMap<String, List<Alignments.AlignmentEntry>>();
        currentPosition = -1;
    }

    /**
     * Get if a new comparison will be made at each new position. Setting this to false will require all alignments
     * to be loaded into memory before comparisons will be made which could consume a LOT of memory.
     * @return if ...
     */
    public boolean isCompareAtEachNewPosition() {
        return compareAtEachNewPosition;
    }

    /**
     * Set if a new comparison will be made at each new position. Setting this to false will require all alignments
     * to be loaded into memory before comparisons will be made which could consume a LOT of memory.
     * @param compareAtEachNewPosition if ...
     */
    public void setCompareAtEachNewPosition(final boolean compareAtEachNewPosition) {
        this.compareAtEachNewPosition = compareAtEachNewPosition;
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
        if (compareAtEachNewPosition && source.getAlignmentStart() != currentPosition) {
            result = makeComparisons();
        }
        sources.add(source);
        addToSamMap(dest);
        if (gobyDest != null) {
            addToGobyMap(gobyDest);
        }
        currentPosition = source.getAlignmentStart();
        return result;
    }

    private void addToSamMap(final SAMRecord dest) {
        final String key = dest.getReadName();
        List<SAMRecord> list = destsMap.get(key);
        if (list == null) {
            list = new LinkedList<SAMRecord>();
            destsMap.put(key, list);
        }
        list.add(dest);
    }

    private void addToGobyMap(final Alignments.AlignmentEntry dest) {
        final String key = dest.getReadName();
        List<Alignments.AlignmentEntry> list = gobyDestsMap.get(key);
        if (list == null) {
            list = new LinkedList<Alignments.AlignmentEntry>();
            gobyDestsMap.put(key, list);
        }
        list.add(dest);
    }

    private void removeFromSamMap(final SAMRecord dest) {
        final String key = dest.getReadName();
        final List<SAMRecord> list = destsMap.get(key);
        list.remove(dest);
    }

    private void removeFromGobyMap(final Alignments.AlignmentEntry dest) {
        if (dest != null) {
            final String key = dest.getReadName();
            final List<Alignments.AlignmentEntry> list = gobyDestsMap.get(key);
            list.remove(dest);
        }
    }

    @Override
    public int finished() {
        return makeComparisons();
    }

    private void resetForPosition() {
        sources.clear();
        destsMap.clear();
        gobyDestsMap.clear();
    }

    /**
     * Make comparisons of every record at a specific alignment start position.
     * This method REQUIRES readName() is preserved when making compact-reads file.
     * @return The total of the comparison failures for all of the records compared
     */
    private int makeComparisons() {
        int resultSum = 0;
        if (!sources.isEmpty()) {
            for (final SAMRecord source : sources) {
                final List<SamAndGobyDestPair> toCompares = findDestPairsForSource(source);
                if (toCompares.isEmpty()) {
                    // Found no comparisons
                    readNum++;
                    comparisonFailureCount++;
                    System.out.println("WARNING: Couldn't find any records to compare against source.readName=" +
                            source.getReadName());
                } else {
                    countComparisonFailures = false;
                    final boolean tempOutputFailedComparisons = outputFailedComparisons;
                    outputFailedComparisons = false;
                    int leastDiffIndex = -1;
                    int leastDiffScore = Integer.MAX_VALUE;
                    int currentDestIndex = 0;
                    for (final SamAndGobyDestPair toCompare : toCompares) {
                        final int currentDiffScore = super.compare(source, toCompare.dest, toCompare.gobyDest);
                        if (currentDiffScore < leastDiffScore) {
                            leastDiffScore = currentDiffScore;
                            leastDiffIndex = currentDestIndex;
                            if (currentDiffScore == 0) {
                                // Stop if we get a diff of 0, it won't get better
                                break;
                            }
                        }
                        currentDestIndex++;
                    }
                    outputFailedComparisons = tempOutputFailedComparisons;
                    countComparisonFailures = true;
                    if (leastDiffScore > 0) {
                        if (outputFailedComparisons) {
                            // Output the failed but still best failure here.
                            super.compare(source,
                                    toCompares.get(leastDiffIndex).dest,
                                    toCompares.get(leastDiffIndex).gobyDest);
                        }
                    }
                    removeFromSamMap(toCompares.get(leastDiffIndex).dest);
                    removeFromGobyMap(toCompares.get(leastDiffIndex).gobyDest);
                }
            }
        }
        for (final Map.Entry<String, List<SAMRecord>> entry : destsMap.entrySet()) {
            for (final SAMRecord dest : entry.getValue()) {
                System.out.println("WARNING: Remaining dest.readName=" + dest.getReadName());
            }
        }
        /*
        // For spliced alignments, there WILL BE remaining AlignmentEntries since a splice has 2 or more
        // alignment entries. Don't report this.
        for (final Map.Entry<String, List<Alignments.AlignmentEntry>> entry : gobyDestsMap.entrySet()) {
            for (final Alignments.AlignmentEntry dest : entry.getValue()) {
                System.out.println("WARNING: Remaining gobyDest.readName=" + dest.getReadName());
            }
        }
        */
        resetForPosition();
        return resultSum;
    }

    /**
     * A sam destination and goby destination pair for comparison (gobyDest may be null).
     */
    private static class SamAndGobyDestPair {
        SAMRecord dest;
        Alignments.AlignmentEntry gobyDest;
        SamAndGobyDestPair(final SAMRecord dest, final Alignments.AlignmentEntry gobyDest) {
            this.dest = dest;
            this.gobyDest = gobyDest;
        }
    }

    /**
     * Make a list of all the possible sam/goby destination comparisons that are possible to assist
     * with finding the best combination of expected/actual/gobyActual for the
     * given position.
     * @param source the source to find comparison pairs for 
     * @return the list of all actual/gobyActual combinations
     */
    private List<SamAndGobyDestPair> findDestPairsForSource(final SAMRecord source) {
        final List<SAMRecord> localDests = findDestsForSource(source);
        final List<Alignments.AlignmentEntry> localGobyDests = findGobyDestsForSource(source);
        final List<SamAndGobyDestPair> result = new LinkedList<SamAndGobyDestPair>();
        for (final SAMRecord dest : localDests) {
            if (localGobyDests.isEmpty()) {
                result.add(new SamAndGobyDestPair(dest, null));
            } else {
                for (final Alignments.AlignmentEntry gobyDest : localGobyDests) {
                    result.add(new SamAndGobyDestPair(dest, gobyDest));
                }
            }
        }
        return result;
    }

    /**
     * Of of the destinations (actual) that have the same read name as source.
     * @param source the source
     * @return the list of actuals that have the same read name.
     */
    private List<SAMRecord> findDestsForSource(final SAMRecord source) {
        final List<SAMRecord> result = destsMap.get(source.getReadName());
        if (result == null) {
            return EMPTY_SAM_LIST;
        } else {
            return result;
        }
    }

    /**
     * Of of the goby destinations (gobyActual) that have the same read name as source.
     * @param source the source
     * @return the list of gobyActuals that have the same read name. if there are no gobyDests this will return
     * and empty list.
     */
    private List<Alignments.AlignmentEntry> findGobyDestsForSource(final SAMRecord source) {
        final List<Alignments.AlignmentEntry> result = gobyDestsMap.get(source.getReadName());
        if (result == null) {
            return EMPTY_GOBY_LIST;
        } else {
            return result;
        }
    }
}
