/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
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
        StringBuffer sb = new StringBuffer();
        sb.append(start);
        sb.append("-");
        sb.append(end);

        return sb.toString();
    }
}
