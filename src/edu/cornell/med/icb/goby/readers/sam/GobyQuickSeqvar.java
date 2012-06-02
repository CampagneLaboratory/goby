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

import edu.cornell.med.icb.goby.util.pool.Resettable;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.bytes.ByteList;
import it.unimi.dsi.lang.MutableString;

/**
 * Class to assist with parsing SAM to Goby, making sequence variations.
 * SequenceVariation objects will be created from these objects.
 */
public class GobyQuickSeqvar implements Resettable {
    MutableString from;
    MutableString to;
    private int readIndex;
    private int position;
    ByteList toQuals;
    int lastIndexPosition;  // bookkeeping

    public GobyQuickSeqvar() {
        from = new MutableString();
        to = new MutableString();
        toQuals = new ByteArrayList();
        reset();
    }

    @Override
    public void reset() {
        from.length(0);
        to.length(0);
        toQuals.clear();
        setReadIndex(1);
        setPosition(0);
        lastIndexPosition = Integer.MAX_VALUE;
    }

    public String getFrom() {
        return from.toString();
    }

    public String getTo() {
        return to.toString();
    }

    public ByteList getToQualities() {
        return toQuals;
    }

    public byte[] getToQualitiesAsBytes() {
        return toQuals.toByteArray();
    }

    public int getReadIndex() {
        return readIndex;
    }

    public int getPosition() {
        return position;
    }

    public String toString() {
        String s = "sequence_variations {\n";
        s += "  to: \"" + to.toString() + "\"\n";
        s += "  from: \"" + from.toString() + "\"\n";
        s += "  position: " + getPosition() + "\n";
        s += "  to_quality: \"" + toQuals.toString() + "\"\n";
        s += "  read_index: " + getReadIndex() + "\n";
        s += "}\n";
        return s;
    }

    public final void setPosition(final int position) {
        assert position>=0 : "position must be positive. seqVar.position must be >=0.";
        this.position = position;
    }

    public final void setReadIndex(final int readIndex) {
        assert readIndex>=1 : "readIndex cannot be zero. seqVar.readIndex is one-based.";

        this.readIndex = readIndex;
    }
}
