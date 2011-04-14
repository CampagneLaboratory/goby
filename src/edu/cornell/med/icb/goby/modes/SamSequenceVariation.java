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

package edu.cornell.med.icb.goby.modes;

import it.unimi.dsi.lang.MutableString;

import java.util.List;
import java.util.LinkedList;

/**
 * Store details of a SequenceVariation. Used by SamHelper.
 */
class SamSequenceVariation {
    private MutableString fromString;
    private MutableString toString;
    private MutableString qual;
    private boolean hasQual;
    private int readIndex;
    private int refPosition;
    private int lastRefPosition;

    public SamSequenceVariation(int refPosition, char refChar, int readIndex, char readChar, boolean hasQual, char qualChar) {
        this.refPosition = refPosition;
        this.lastRefPosition = refPosition;
        this.fromString = new MutableString();
        this.fromString.append(refChar);
        this.readIndex = readIndex;
        this.toString = new MutableString();
        this.toString.append(readChar);
        this.hasQual = hasQual;
        if (hasQual) {
            this.qual = new MutableString();
            this.qual.append(qualChar);
        }
    }

    public MutableString getFromString() {
        return fromString;
    }

    public MutableString getToString() {
        return toString;
    }

    public MutableString getQual() {
        return qual;
    }

    public boolean isHasQual() {
        return hasQual;
    }

    public int getReadIndex() {
        return readIndex;
    }

    public int getRefPosition() {
        return refPosition;
    }

    public String toString() {
        MutableString output = new MutableString();
        output.append(String.format("%s[%d]->%s[%d] / qual=",
                fromString.toString(), refPosition,
                toString.toString(), readIndex));
        if (hasQual) {
            for (int i = 0; i < qual.length(); i++) {
                if (i > 0) {
                    output.append(":");
                }
                output.append(String.format("%d", (int) qual.charAt(i)));
            }
        } else {
            output.append("none");
        }
        return output.toString();
    }

    /**
     * This takes a List[SequenceVariation] which are all fromString and toString are all single
     * base and merge them so contiguous sequence variations are in one SamSequenceVariation.
     * Insertions and Deletions always start a new SamSequenceVariation.
     * @param vars the SamSequenceVariation list
     */
    static void merge(final List<SamSequenceVariation> vars) {
        List<SamSequenceVariation> toRemoves = new LinkedList<SamSequenceVariation>();
        SamSequenceVariation current = null;
        for (final SamSequenceVariation var : vars) {
            // This should only be called if the sequence variations are
            // all a single base. Pre-check this assumption.
            if (var.fromString.length() > 1) {
                return;
            }
        }
        for (final SamSequenceVariation var : vars) {
            if (current == null) {
                current = var;
            } else {
                if ((var.refPosition == current.lastRefPosition) || ((var.refPosition) == current.lastRefPosition + 1)) {
                    final char currentFromChar = current.fromString.charAt(current.fromString.length() - 1);
                    final char currentToChar = current.toString.charAt(current.toString.length() - 1);
                    final char varFromChar = var.fromString.charAt(0);
                    final char varToChar = var.toString.charAt(0);
                    if ((varFromChar == '-' && currentFromChar != '-') ||
                        (varToChar == '-' && currentToChar != '-') ||
                        (var.hasQual != current.hasQual)) {
                        // Insertions / deletions start new sequence variations and shouldn't be merged with previous
                        current = var;
                    } else {
                        // Merge with previous
                        current.fromString.append(var.fromString);
                        current.toString.append(var.toString);
                        if (current.hasQual) {
                            current.qual.append(var.qual);
                        }
                        toRemoves.add(var);
                        current.lastRefPosition = var.refPosition;
                    }
                } else {
                    // Not continguous
                    current = var;
                }
            }
        }
        for (final SamSequenceVariation toRemove : toRemoves) {
            vars.remove(toRemove);
        }
    }

    public boolean equals(final int refPosition, final CharSequence refChars, final int readIndex,
                          final CharSequence readChars, final int[] qualChars) {

        if (refPosition != this.refPosition || readIndex != this.readIndex) {
            return false;
        }
        if (!this.fromString.equals(refChars) || !this.toString.equals(readChars)) {
            return false;
        }
        if (hasQual) {
            if (qualChars == null || qual.length() != qualChars.length) {
                return false;
            }
            for (int i = 0; i < qual.length(); i++) {
                if (qual.charAt(i) != qualChars[i]) {
                    return false;
                }
            }
        } else {
            if (qualChars != null && qualChars.length > 0){
                return false;
            }
        }
        return true;
    }

    public static boolean contains(final List<SamSequenceVariation> vars, final int refPosition,
                                   final CharSequence refChars, final int readIndex, final CharSequence readChars,
                                   final int[] qualChars) {

        for (final SamSequenceVariation var : vars) {
            if (var.equals(refPosition, refChars, readIndex, readChars, qualChars)) {
                return true;
            }
        }
        return false;
    }
}
