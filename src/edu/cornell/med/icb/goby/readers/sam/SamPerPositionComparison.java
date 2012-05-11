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

    //@Override
    public boolean comparex(final SAMRecord source, final SAMRecord dest, final Alignments.AlignmentEntry gobyDest) {
        if (source.getAlignmentStart() != currentPosition) {
            System.out.println("-----------------");
        }
        currentPosition = source.getAlignmentStart();
        System.out.println("s:" + usableReadOf(source));
        System.out.println("d:" + usableReadOf(dest));
        return true;
    }

    @Override
    public boolean compare(final SAMRecord source, final SAMRecord dest, final Alignments.AlignmentEntry gobyDest) {
        boolean result = true;
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
    public void finished() {
        makeComparisons();
    }

    private void resetForPosition() {
        sources.clear();
        dests.clear();
        gobyDests.clear();
        destclippedReadToIndexesMap.clear();
    }

    private boolean makeComparisons() {
        boolean result = true;
        if (!sources.isEmpty()) {
            for (final SAMRecord source : sources) {
                SAMRecord dest = null;
                Alignments.AlignmentEntry gobyDest = null;
                if (dests.size() == 1) {
                    dest = dests.get(0);
                    if (!gobyDests.isEmpty()) {
                        gobyDest = gobyDests.get(0);
                    }
                } else {
                    //
                    // TODO: The logic here is a little tricky. When multiple items match at the same position,
                    // TODO: which is common, AND they have the exact same read since we don't have the read name.
                    // TODO: If we are keeping qualities, we can use the qualities + clippedReads to do a really
                    // TODO: good lineup of read-to-read but without qualities, we use cigar which is NOT good
                    // TODO: at lining them up.
                    final String sourceClippedRead = usableReadOf(source);
                    final LinkedList<Integer> indexesList = destclippedReadToIndexesMap.get(sourceClippedRead);
                    if (indexesList != null) {
                        // Get the index of the first unused matching destination that contains the same
                        // exact usable read as the source.
                        int foundDestIndex = -1;
                        for (final int localDestIndex : indexesList) {
                            dest = dests.get(localDestIndex);
                            if (isMappedQualitiesPreserved()) {
                                if (source.getBaseQualityString().equals(dest.getBaseQualityString())) {
                                    foundDestIndex = localDestIndex;
                                    break;
                                }
                            } else if (dest.getCigarString().equals(source.getCigarString())) {
                                foundDestIndex = localDestIndex;
                                break;
                            }
                        }
                        if (foundDestIndex != -1) {
                            // Found the best destination record
                            indexesList.removeFirstOccurrence(foundDestIndex);
                            if (!gobyDests.isEmpty()) {
                                gobyDest = gobyDests.get(foundDestIndex);
                            }
                        }
                    }
                }
                if (dest != null) {
                    if (!super.compare(source, dest, gobyDest)) {
                        result = false;
                    }
                } else {
                    readNum++;
                    comparisonFailureCount++;
                    System.out.println("ERROR: Couldn't find a dest record to compare against for readName=" + source.getReadName());
                    // Consider: When this is done, pick the sources that didn't get compared with the dest's that
                    // didn't get compared. Often there may just be 1 and 1??
                    result = false;
                }
            }
        }
        resetForPosition();
        return result;
    }
}
