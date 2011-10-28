/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.data;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.io.PrintWriter;
import java.util.Collections;

/**
 * Represent annotations along reference sequences. Each annotation can contain several
 * intervals along the reference (called annotation segments).
 */
public class Annotation implements Comparable<Annotation> {
    protected final String id;
    protected final String chromosome;
    protected final ObjectList<Segment> segments;
    protected final String strand;

    public Annotation(final String id, final String chromosome) {
        this(id, chromosome, "N/A");
    }

    public Annotation(final String id, final String chromosome, final String strand) {
        super();
        segments = new ObjectArrayList<Segment>();
        this.id = id;
        this.chromosome = chromosome;
        this.strand = strand;
    }

    public void sortSegments() {
        Collections.sort(segments);
    }

    public void addSegment(final Segment segment) {
        segments.add(segment);
    }

    public int compareTo(final Annotation annot) {
        return this.getStart() - annot.getStart();
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null) {
            return false;
        } else if (obj == this) {
            return true;
        } else if (!(obj instanceof Annotation)) {
            return false;
        } else {
            final Annotation annot = (Annotation) obj;
            return ((this.getStart() == annot.getStart()) && (this.getEnd() == annot.getEnd()));
        }
    }

    public int getStart() {
        return segments.size()==0? -1: segments.get(0).getStart();
    }

    public int getEnd() {
        return segments.size()==0? -1:segments.get(segments.size() - 1).getEnd();
    }
    int length = -1;
    public int getLength() {
        if (length != -1) {
            return length;
        }
        int result = 0;
        for (final Segment segment : segments) {
            result += segment.getEnd() - segment.getStart() + 1;
        }
        length = result;
        return result;
    }

    public String getId() {
        return id;
    }

    public String getChromosome() {
        return chromosome;
    }

    public ObjectList<Segment> getSegments() {
        return segments;
    }

    public String getStrand() {
        return strand;
    }

    public boolean overlap(final Annotation annotation2) {
        boolean overlap = false;
        // 4-cases to consider

        // 1st: annotation2's start is between this annotation start and end AND annotation2's end is after this end
        if ((this.getStart() <= annotation2.getStart() && annotation2.getStart() <= this.getStart()) &&
                (annotation2.getEnd() >= this.getEnd())) {
            overlap = true;
        }
        //2nd: this start is between annotation2 start and end AND this end is after annotation2 end
        else if ((annotation2.getStart() <= this.getStart() && this.getStart() <= annotation2.getEnd()) &&
                (this.getEnd() >= annotation2.getEnd())) {
            overlap = true;

        }
        //3rd: this is contained within annotation2
        else if ((annotation2.getStart() <= this.getStart() && this.getStart() <= annotation2.getEnd()) &&
                (annotation2.getStart() <= this.getEnd() && this.getEnd() <= annotation2.getEnd())) {
            overlap = true;

        }
        //4th: annotation2 is contained within this
        else if ((this.getStart() <= annotation2.getStart() && annotation2.getStart() <= this.getEnd()) &&
                (this.getStart() <= annotation2.getEnd() && annotation2.getEnd() <= this.getEnd())) {
            overlap = true;
        }

        return overlap;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append("[ ");
        sb.append(chromosome);
        sb.append(' ');
        sb.append(id);
        sb.append(' ');
        for (final Segment segment : segments) {
            sb.append(segment);
            sb.append(' ');
        }
        sb.append(" ]");
        return sb.toString();
    }

    public void write(final PrintWriter annotationWriter) {
        final char delimiter = '\t';
        write(annotationWriter, delimiter);
    }

    public void write(final PrintWriter annotationWriter, final char delimiter) {
        // Chromosome Name Strand  Ensembl Gene ID Ensembl Exon ID Exon Chr Start (bp) Exon Chr End (bp)
        for (final Segment segment : segments) {
            annotationWriter.write(chromosome);
            annotationWriter.write(delimiter);
            annotationWriter.write(strand);
            annotationWriter.write(delimiter);
            annotationWriter.write(id);
            annotationWriter.write(delimiter);
            annotationWriter.write(segment.getId());
            annotationWriter.write(delimiter);
            annotationWriter.write(Integer.toString(segment.getStart()));
            annotationWriter.write(delimiter);
            annotationWriter.write(Integer.toString(segment.getEnd()));
            annotationWriter.write('\n');
        }
    }

    public void remove(Segment element) {
        segments.remove(element);
    }
}
