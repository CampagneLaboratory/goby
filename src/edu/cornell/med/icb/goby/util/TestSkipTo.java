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

package edu.cornell.med.icb.goby.util;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.util.Random;

/**
 * @author Fabien Campagne
 *         Date: May 8, 2011
 *         Time: 9:43:09 AM
 */
public class TestSkipTo {
    public static void main(final String[] args) throws IOException {
        if (args.length == 0) {
            System.err.println("This utility requires exactly one arguments: the filename of an alignment file to use for testing.");
        }
        String filename = args[0];
        TestSkipTo processor = new TestSkipTo();
        processor.check(filename);
    }

    private void check(String filename) throws IOException {

        AlignmentReader reader = new AlignmentReaderImpl(filename);
        reader.readHeader();
        ObjectList<location> locations = new ObjectArrayList<location>(reader.getNumberOfAlignedReads());
        int chunkIndex = 0;
        int entryIndex = 0;
        int maxPerChunk = 10000;
        System.out.println("Loading locations from alignment");

        for (Alignments.AlignmentEntry entry : reader) {

            locations.add(new location(entry.getTargetIndex(), entry.getPosition(), chunkIndex, entryIndex,
                    entry.getQueryIndex()));
            entryIndex++;
            chunkIndex = entryIndex / maxPerChunk;
        }
        reader.close();
        checkSlices(locations, filename);
        checkSkipTo(locations, filename);
    }

    Random random = new Random();

    private void checkSkipTo(ObjectList<location> locations, String filename) throws IOException {
        int numTrials = 1000;
        System.out.printf("Checking random skipTo %n");

        AlignmentReader reader = null;
        reader = new AlignmentReaderImpl(filename);
        for (int i = 0; i < numTrials; i++) {


            int randomLocationIndex = chooseRandom(random, 0, locations.size() - 1);
            location loc = locations.get(randomLocationIndex);
            loc = adjustStart(locations, randomLocationIndex, loc);
            reader.reposition(loc.targetIndex, loc.position);
            Alignments.AlignmentEntry entry = reader.skipTo(loc.targetIndex, loc.position);
            assertLocation(entry, locations, loc, loc.locationIndex);

        }
        if (reader != null) reader.close();
    }

    /**
     * @param lo lower limit of range
     * @param hi upper limit of range
     * @return a random integer in the range <STRONG>lo</STRONG>,
     *         <STRONG>lo</STRONG>+1, ... ,<STRONG>hi</STRONG>
     */
    private int chooseRandom(Random random, final int lo, final int hi) {
        final double r = random.nextDouble();
        int result = (int) ((long) lo + (long) ((1L + (long) hi - (long) lo) * r));
        assert result >= lo && result <= hi;
        return result;
    }

    private void checkSlices(ObjectList<location> locations, String filename) throws IOException {
        for (int i = 0; i < 10; i++) {
            int startLocationIndex = i * locations.size() / 10;
            int endLocationIndex = startLocationIndex + locations.size() / 20;
            missingLocations = new WarningCounter(20);
            mismatchReference = new WarningCounter(20);
            mismatchPosition = new WarningCounter(20);
            mismatchQueryIndex = new WarningCounter(20);
            location start = locations.get(startLocationIndex);
            location end = locations.get(endLocationIndex);

            start = adjustStart(locations, startLocationIndex, start);
            end = adjustEnd(locations, endLocationIndex, end);
            startLocationIndex = start.locationIndex;
            endLocationIndex = end.locationIndex;

            System.out.printf("Checking slice %s to %s%n", start, end);
            int locationIndex = startLocationIndex;

            AlignmentReader reader = new AlignmentReaderImpl(filename,
                    start.targetIndex, start.position,
                    end.targetIndex, end.position);

            for (Alignments.AlignmentEntry entry : reader) {

                assertLocation(entry, locations, locations.get(locationIndex), locationIndex);
                locationIndex++;
            }
            if (locationIndex != endLocationIndex + 1) {
                missingLocations.warn(LOG, "Missed %d entries at end of slice, last location seen was %s index=%d expected %s index=%d %n",
                        endLocationIndex - locationIndex,
                        locations.get(locationIndex),
                        locationIndex,
                        locations.get(endLocationIndex), endLocationIndex);
            }
        }
    }

    WarningCounter missingLocations = new WarningCounter(20);
    WarningCounter mismatchQueryIndex = new WarningCounter(20);
    WarningCounter mismatchReference = new WarningCounter(20);
    WarningCounter mismatchPosition = new WarningCounter(20);

    private location adjustEnd(ObjectList<location> locations, int endLocationIndex, location end) {
        // advance 'end' to the last entry that has the same location as 'end'.
        location newPos;
        for (; ;) {
            newPos = locations.get(endLocationIndex);
            if (newPos.targetIndex != end.targetIndex) break;
            if (newPos.position != end.position) break;
            if (endLocationIndex == locations.size()) break;
            end = newPos;
            ++endLocationIndex;


            //  System.out.printf("adjusting locationIndex=%d end=%s%n", endLocationIndex, end);
        }

        return end;
    }

    private location adjustStart(ObjectList<location> locations, int startLocationIndex, location start) {
        // rewind 'start' to the first entry that has the same location as 'start'.

        location newPos;
        for (; ;) {
            newPos = locations.get(startLocationIndex);
            if (newPos.targetIndex != start.targetIndex) break;
            if (newPos.position != start.position) break;
            if (startLocationIndex == 0) break;
            start = newPos;
            --startLocationIndex;


            //    System.out.printf("adjusting locationIndex=%d start=%s%n", startLocationIndex, start);
        }


        return start;

    }

    private void assertLocation(Alignments.AlignmentEntry entry, ObjectList<location> locations, location location, int locationIndex) {
        if (location.queryIndex != entry.getQueryIndex()) {
            mismatchQueryIndex.warn(LOG, "queryIndex must match for entry (%d,%d) at location %d (%d,%d)", entry.getQueryIndex(),
                    entry.getPosition(),
                    locationIndex, location.targetIndex, location.position);
        }

        if (location.targetIndex != entry.getTargetIndex()) {
            mismatchReference.warn(LOG, "targetIndex must match for entry (%d,%d) at location %d (%d,%d)", entry.getTargetIndex(),
                    entry.getPosition(),
                    locationIndex, location.targetIndex, location.position);
        }
        if (location.position != entry.getPosition()) {
            mismatchPosition.warn(LOG, "position must match for entry (%d,%d) at location %d (%d,%d)", entry.getTargetIndex(),
                    entry.getPosition(),
                    locationIndex, location.targetIndex, location.position);
        }
    }

    class location {
        int targetIndex;
        int position;

        int chunkIndex;
        int locationIndex;
        private int queryIndex;

        location(int referenceIndex, int position, int chunkIndex, int locationIndex, int queryIndex) {
            this.targetIndex = referenceIndex;
            this.position = position;
            this.chunkIndex = chunkIndex;
            this.locationIndex = locationIndex;
            this.queryIndex = queryIndex;
        }

        @Override
        public String toString() {
            return String.format("(%d,%d) in chunk %d query-index=%d", targetIndex, position, chunkIndex, queryIndex);
        }

    }

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(TestSkipTo.class);

}


