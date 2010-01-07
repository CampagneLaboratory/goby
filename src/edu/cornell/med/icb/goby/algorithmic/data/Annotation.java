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

import java.io.PrintWriter;
import java.util.Collections;

public class Annotation {
    public final String id;
    public final String chromosome;
    public final ObjectList<Segment> segments;
    public String strand;

    public Annotation(final String id, final String chromosome) {
        segments = new ObjectArrayList<Segment>();
        this.id = id;
        this.chromosome = chromosome;
        strand="N/A";
    }

    public void sortSegments() {
        Collections.sort(segments);
    }

    public void addSegment(final Segment segment) {
        segments.add(segment);
    }

    public int getStart() {
        return segments.get(0).start;
    }

    public int getEnd() {
        return segments.get(segments.size() - 1).end;
    }

    public int getLength() {
        int result = 0;
        for (Segment s : segments) {
            result += s.end - s.start + 1;
        }
        return result;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("[ ");
        sb.append(chromosome);
        sb.append(' ');
        sb.append(id);
        sb.append(" ");
        for (Segment segment : segments) {
            sb.append(segment);
            sb.append(" ");
        }
        sb.append(" ]");
        return sb.toString();
    }

    public void write(PrintWriter annotationWriter) {
        //Chromosome Name Strand  Ensembl Gene ID Ensembl Exon ID Exon Chr Start (bp)     Exon Chr End (bp)
        for (Segment segment : segments) {
            annotationWriter.write(chromosome);
            annotationWriter.write('\t');
            annotationWriter.write(strand);
            annotationWriter.write('\t');
            annotationWriter.write(id);
            annotationWriter.write('\t');
            annotationWriter.write(segment.id);
            annotationWriter.write('\t');
            annotationWriter.write(Integer.toString(segment.start));
            annotationWriter.write('\t');
            annotationWriter.write(Integer.toString(segment.end));
            annotationWriter.write('\n');
        }
    }
}
