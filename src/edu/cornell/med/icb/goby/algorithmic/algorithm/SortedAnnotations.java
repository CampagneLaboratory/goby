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

package edu.cornell.med.icb.goby.algorithmic.algorithm;

import com.sun.xml.internal.bind.v2.TODO;
import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.modes.CompactAlignmentToAnnotationCountsMode;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;

/**
 * A set of annotations sorted in genome order. Provide the capability of returning the set of annotations that
 * overlap with a given position.
 *
 * @author Fabien Campagne
 *         Date: 12/12/11
 *         Time: 4:48 PM
 */
public class SortedAnnotations {


    /**
     * Set of annotations ordered in increasing chromosome position. The order must be consistent with the traversal
     * order of the genome.
     */
    private Annotation[] annotations;

    private int annotationIndex;

    private RandomAccessSequenceInterface genome;


    public SortedAnnotations() {
        annotationIndices = new IntAVLTreeSet();
        annotationIndex=0;
    }

    public void setAnnotations(Annotation[] annotations) {
            this.annotations = annotations;
        }

    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome=genome;
    }

    /**
     * Index of the active annotations. Active annotations are those with start before position and end after position.
     * All active annotations must be considered for overlap with every position. It is expected that only a subset of
     * active annotations will have segments that overlap with the position.
     */
    private IntAVLTreeSet annotationIndices;

    /**
     * Load and sort annotations by their end position.
     *
     * @param filename
     * @throws java.io.IOException
     */
    public void loadAnnotations(final String filename) throws IOException {

        final ObjectArrayList<Annotation> result = new ObjectArrayList<Annotation>();
        final Object2ObjectMap<String, ObjectList<Annotation>> map = CompactAlignmentToAnnotationCountsMode.readAnnotations(filename);
        for (int referenceIndex = 0; referenceIndex < genome.size(); referenceIndex++) {
            final ObjectList<Annotation> list = map.get(genome.getReferenceName(referenceIndex));
            Collections.sort(list, compareAnnotationEnd);
            result.addAll(list);
        }
        annotations = new Annotation[result.size()];
        result.toArray(annotations);
    }

    /**
     * Determine if there are any more annotations to consider.
     *
     * @return
     */
    public boolean hasMoreAnnotations() {
        return annotationIndex < annotations.length;
    }

    public Annotation nextAnnotation() {
        return annotations[annotationIndex];
    }

    public void advanceToPosition(final String chromosome, final int pos) {
        while (hasMoreAnnotations()) {
            final Annotation ann = nextAnnotation();
            if (chromosome1StrictlyBefore2(ann.getChromosome(), chromosome)) {
                continue;
            }
            if (chromosome1StrictlyBefore2(chromosome, ann.getChromosome())) {
                // we are past chromosome. Stop immediately and backup.
                annotationIndex -= 1;
                break;
            }

            // ann is located on the query chromosome.
            if (ann.overlap(chromosome, pos)) {
                annotationIndices.add(annotationIndex);
            }

            advanceToNextAnnotation();
        }

    }


    public boolean hasOverlappingAnnotations(String chrom, int pos) {
        annotationIndices.clear();
        advanceToPosition(chrom, pos);
        return !(annotationIndices.isEmpty());
    }

    private ObjectArrayList<Annotation> set = new ObjectArrayList<Annotation>();

    /**
     * Return the set of annotations that overlap with the current position.
     *
     * @return
     */
    public ObjectArrayList<Annotation> currentAnnotations() {
        set.clear();
        for (final int index : annotationIndices) {
            set.add(annotations[index]);
        }
        return set;
    }

    /**
     * Determine if (currentChromosome,pos) is a position past the end of the current annotation.
     *
     * @param currentChromosome chromosome of VCF.
     * @param pos               position of VCF.
     * @return true when position of VCF is past the current annotation.
     */
    public boolean pastCurrentAnnotation(final String currentChromosome, final int pos) {
        final Annotation annotation = annotations[annotationIndex];

        if (!annotation.getChromosome().equals(currentChromosome)) {
            // chromosome differ:
            if (chromosome1StrictlyBefore2(annotation.getChromosome(), currentChromosome)) {
                // currentChromosome is past the chromosome of the current annotation.
                return true;
            } else {
                // position is past the end of the current annotation.
                return pos > annotation.getEnd();
            }
        } else {
            // chromosomes match
            // is position past the end of the current annotation?
            return pos > annotation.getEnd();
        }
    }

    /**
     * Returns true if chromosome c1 occurs strictly before chromosome c2 in the genome.
     *
     * @param c1 chromosome 1
     * @param c2 chromosome 2
     * @return true when c1 occurs strictly before c2
     */
    private boolean chromosome1StrictlyBefore2(final String c1, final String c2) {
        final int c1Index = genome.getReferenceIndex(c1);
        final int c2Index = genome.getReferenceIndex(c2);
        return c1Index < c2Index;
    }


    public void advanceToNextAnnotation() {
        annotationIndex += 1;
    }

    private final Comparator<? super Annotation> compareAnnotationEnd = new Comparator<Annotation>() {
        @Override
        public int compare(final Annotation annotation, final Annotation annotation1) {
            return annotation.getEnd() - annotation1.getEnd();
        }
    };

}
