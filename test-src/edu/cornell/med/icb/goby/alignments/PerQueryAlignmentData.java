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

import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * Class used to assist with observing sequence variations for an
 * individual query index. These can then be compared to make sure
 * two alignments have the same sequence variations.
 */
public class PerQueryAlignmentData {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(PerQueryAlignmentData.class);

    public String query;           // From compact-reads file
    public int queryLength;        // From compact-reads file, verified against alignment
    public int minRefPosition;     // From seq-var
    public int maxRefPosition;     // From seq-var
    public int minReadIndex;       // From seq-var
    public int maxReadIndex;       // From seq-var
    public int numDeletions;       // From seq-var, from display-seq-var TSV file
    public int numInsertions;      // From seq-var, from display-seq-var TSV file
    public int numMismatches;      // From seq-var, from display-seq-var TSV file
    public int refLength;          // Calculated, should equal queryLength
    public int readLength;         // Calculated, should equal queryLength
    public int targetPosition;     // From alignment
    public int queryPosition;      // From alignment
    public int targetAlignedLength; // From alignment
    public int queryAlignedLength;   // From alignment
    public boolean reverseStrand;  // From alignment
    public int firstReadIndex;     // From alignment
    public int leftPadding;        // Calculated
    public int rightPadding;       // Calculated
    public Object2ObjectMap<String, String> refPositionReadIndexToBaseMap;
    public PerQueryAlignmentData() {
        query = "";
        queryLength = 0;
        minRefPosition = Integer.MAX_VALUE;
        maxRefPosition = Integer.MIN_VALUE;
        minReadIndex = Integer.MAX_VALUE;
        maxReadIndex = Integer.MIN_VALUE;
        numDeletions = 0;
        numInsertions = 0;
        numMismatches = 0;
        refLength = 0;
        readLength = 0;
        targetPosition = 0;
        queryPosition = 0;
        targetAlignedLength = 0;
        queryAlignedLength = 0;
        reverseStrand = false;
        firstReadIndex = -1;
        queryPosition = -1;
        rightPadding = -1;
        refPositionReadIndexToBaseMap = new Object2ObjectOpenHashMap<String, String>();
    }
    public void observe(int refPosition, int readIndex) {
        minRefPosition = Math.min(refPosition, minRefPosition);
        maxRefPosition = Math.max(refPosition, maxRefPosition);
        minReadIndex = Math.min(readIndex, minReadIndex);
        maxReadIndex = Math.max(readIndex, maxReadIndex);
        refLength = maxRefPosition - minRefPosition + numInsertions + 1;
        readLength = maxReadIndex - minReadIndex + numDeletions + 1;
    }
    public void observe(int refPosition, int readIndex, char fromBase, char toBase) {
        if (toBase == '-') {
            numDeletions++;
        } else if (fromBase == '-') {
            numInsertions++;
        } else {
            numMismatches++;
        }
        observe(refPosition, readIndex);
        refPositionReadIndexToBaseMap.put(
                String.format("%d:%d", refPosition, readIndex),
                String.format("%c->%c", fromBase, toBase));
        leftPadding = queryPosition;
        rightPadding = (queryLength + numDeletions) - (targetAlignedLength + numInsertions) - leftPadding;
    }
    public void setQuery(String query) {
        this.query = query;
        this.queryLength = query.length();
    }

    public String toString() {
        StringBuffer result = new StringBuffer();
        result.append(String.format(
            "{%n  query=%s%n  queryLength=%d%n  minRefPosition=%d%n  maxRefPosition=%d%n  minReadIndex=%d%n" +
            "  maxReadIndex=%d%n  numDeletions=%d%n  numInsertions=%d%n  numMismatches=%d%n  refLength=%d%n" +
            "  readLength=%d%n  targetPosition=%d%n  queryPosition=%d%n  firstReadIndex=%d%n" +
            "  targetAlignedLength=%d%n  queryAlignedLength=%d%n  reverseStrand=%s%n" +
            "  leftPadding=%d%n  rightPadding=%d%n",
                query, queryLength, minRefPosition, maxRefPosition, minReadIndex,
                maxReadIndex, numDeletions, numInsertions, numMismatches, refLength,
                readLength, targetPosition, queryPosition, firstReadIndex,
                targetAlignedLength, queryAlignedLength,
                reverseStrand ? "true" : "false", leftPadding, rightPadding));
        result.append("  seqvars=[");
        final Set<String> keys = refPositionReadIndexToBaseMap.keySet();
        final List<String> keysList = new ArrayList<String>(keys.size());
        keysList.addAll(keys);
        Collections.sort(keysList);
        for (final String key : keysList) {
            final String value = refPositionReadIndexToBaseMap.get(key);
            result.append(String.format("%s %s, ", key, value));
        }
        result.append("]}");
        return result.toString();
    }
}

