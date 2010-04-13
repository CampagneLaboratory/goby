/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.algorithmic.data;

public class Segment implements Comparable<Segment> {
    public final int start;
    public final int end;
    public final String id;
    public final String strand;

    public Segment(final int start, final int end, final String id, final String strand) {
        this.start = start;
        this.end = end;
        this.id = id;
        this.strand = strand;
    }

    public int compareTo(final Segment seg) {
        return start - seg.start;
    }

    public boolean equals(final Segment seg) {
        return start == seg.start;
    }

    /**
     * Return the length of this annotation segment (number of bases in the segment). End-start+1;
     *
     * @return
     */
    public int getLength() {
        return end - start + 1;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append(start);
        sb.append("-");
        sb.append(end);

        return sb.toString();
    }

    public boolean overlap(final Segment segment2) {
        boolean overlap = false;

        //4-cases to consider

        //1st: segment2's start is between this segment start and end AND segment2's end is after this end
        if ((this.start <= segment2.start && segment2.start <= this.start) &&
                (segment2.end >= this.end)) {
            overlap = true;
        }
        //2nd: this start is between segment2 start and end AND this end is after segment2 end
        else if ((segment2.start <= this.start && this.start <= segment2.end) &&
                (this.end >= segment2.end)) {
            overlap = true;

        }
        //3rd: this is contained within segment2
        else if ((segment2.start <= this.start && this.start <= segment2.end) &&
                (segment2.start <= this.end && this.end <= segment2.end)) {
            overlap = true;

        }
        //4th: segment2 is contained within this
        else if ((this.start <= segment2.start && segment2.start <= this.end) &&
                (this.start <= segment2.end && segment2.end <= this.end)) {
            overlap = true;
        }

        return overlap;
    }

    public Segment merge(final Segment mergee) {
        //this segment is the merger
        final int newStart = Math.min(this.start, mergee.start);
        final int newEnd = Math.max(this.end, mergee.end);
        return new Segment(newStart, newEnd, this.id + "_" + mergee.id, this.strand);
    }

    public int distance(final Segment segment2) {
        int dist = 0;
        if (this.overlap(segment2)) {
            dist = 0;
        } else if (segment2.start > this.start) {
            dist = segment2.start - this.end;
        } else {
            dist = this.start - segment2.end;
        }
        return dist;
    }
}
