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

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectListIterator;

import java.io.PrintWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.util.Collections;

import org.apache.commons.io.IOUtils;

public class Annotation implements Comparable<Annotation> {
    public final String id;
    public final String chromosome;
    public final ObjectList<Segment> segments;
    public String strand;

    public Annotation(final String id, final String chromosome) {
        segments = new ObjectArrayList<Segment>();
        this.id = id;
        this.chromosome = chromosome;
        strand = "N/A";
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

    public boolean equals(final Annotation annot) {
        return ((this.getStart() == annot.getStart()) && (this.getEnd() == annot.getEnd()));
    }

    public int getStart() {
        return segments.get(0).start;
    }

    public int getEnd() {
        return segments.get(segments.size() - 1).end;
    }

    public int getLength() {
        int result = 0;
        for (final Segment segment : segments) {
            result += segment.end - segment.start + 1;
        }
        return result;
    }

    public boolean overlap(Annotation annotation2) {
        boolean overlap = false;
        //4-cases to consider

        //1st: annotation2's start is between this annotation start and end AND annotation2's end is after this end
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
        sb.append(" ");
        for (final Segment segment : segments) {
            sb.append(segment);
            sb.append(" ");
        }
        sb.append(" ]");
        return sb.toString();
    }

    public void write(final PrintWriter annotationWriter) {
        char delimiter = '\t';
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
            annotationWriter.write(segment.id);
            annotationWriter.write(delimiter);
            annotationWriter.write(Integer.toString(segment.start));
            annotationWriter.write(delimiter);
            annotationWriter.write(Integer.toString(segment.end));
            annotationWriter.write('\n');
        }
    }
}
