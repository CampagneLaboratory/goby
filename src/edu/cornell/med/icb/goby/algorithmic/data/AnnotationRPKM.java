/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import java.io.PrintWriter;

/**
 * A data structure for storing RPKM value together with an annotation
 *
 * @author Jaaved Mohammed
 */

public class AnnotationRPKM extends Annotation{
    public double RPKM;

    public AnnotationRPKM(final String id, final String chromosome, final double RPKM)
    {
        super(id, chromosome);
        this.RPKM = RPKM;
    }

    @Override
        public void write(final PrintWriter annotationWriter) {
        final char delimiter = '\t';
        write(annotationWriter, delimiter);
    }

    @Override
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
            annotationWriter.write(delimiter);
            annotationWriter.write(Double.toString(RPKM));
            annotationWriter.write('\n');
        }
    }
}
