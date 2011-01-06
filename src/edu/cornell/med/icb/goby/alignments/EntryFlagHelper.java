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

/**
 * @author Fabien Campagne
 *         Date: Jan 6, 2011
 *         Time: 5:02:53 PM
 */
public class EntryFlagHelper {
    /*
      000000001    paired
      000000010    properly paired
      000000100    read unmapped
      000001000    mate unmapped
      000010000    read reverse strand
      000100000    mate reverse strand
      001000000    first in pair
      010000000    second in pair
      100000000    not primary alignment
     */
    /**
     * If this alignment entry is part of a pair.
     *
     * @param entry The entry.
     * @return True if the entry is part of a pair.
     */
    public static boolean isPaired(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & 0x1L) != 0;
    }

    /**
     * If this alignment entry is properly paired. (e.g., mapped in the correct orientation and
     * within reasonable distance on the same reference sequence).
     *
     * @param entry The entry.
     * @return True if the entry is properly paired.
     */
    public static boolean isProperlyPaired(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 1)) != 0;
    }

    public static boolean isReadUnmapped(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 2)) != 0;
    }

    public static boolean isMateUnmapped(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 3)) != 0;
    }

    public static boolean isReadReverseStrand(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 4)) != 0;
    }

    public static boolean isMateReverseStrand(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 5)) != 0;
    }

    public static boolean isFirstInPair(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 6)) != 0;
    }

    public static boolean isSecondInPair(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 7)) != 0;
    }

    public static boolean isNotPrimaryAlignment(final Alignments.AlignmentEntry entry) {
        return (entry.getPairFlags() & (0x1L << 8)) != 0;
    }
}

