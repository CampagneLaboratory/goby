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

import edu.cornell.med.icb.goby.algorithmic.data.Annotation;
import edu.cornell.med.icb.goby.modes.CompactAlignmentToAnnotationCountsMode;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Collections;
import java.util.Comparator;

/**
 * A set of annotations sorted in genome order. Provide the capability of returning the set of annotations that
 * overlap with a given position.
 *
 * @author Fabien Campagne
 * @author Nyasha Chambwe
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

    /**
     * Index of the active annotations. Active annotations are those with start before position and end after position.
     * All active annotations must be considered for overlap with every position. It is expected that only a subset of
     * active annotations will have segments that overlap with the position.
     */
    private IntAVLTreeSet annotationIndicesInRange;

    public IntAVLTreeSet getValidAnnotationIndices() {
        return validAnnotationIndices;
    }

    /**
     * Index of the active annotations that overlap with a given position.
     */
    private IntAVLTreeSet validAnnotationIndices;

    public SortedAnnotations() {
        annotationIndicesInRange = new IntAVLTreeSet();
        validAnnotationIndices = new IntAVLTreeSet();
    }

    public void setAnnotations(Annotation[] annotations) {
        this.annotations = annotations;
    }

    public void setGenome(RandomAccessSequenceInterface genome) {
        this.genome = genome;
    }

    public int getAnnotationIndex() {
        return annotationIndex;
    }


    /**
     * Load and sort annotations by their end position.
     *
     * @param filename
     * @throws java.io.IOException
     */
    public void loadAnnotations(final String filename) throws IOException {
        loadAnnotations(new FileReader(filename));
    }

    /**
     * Load and sort annotations by their start position.
     *
     * @param annotReader
     * @throws java.io.IOException
     */
    public void loadAnnotations(final Reader annotReader) throws IOException {

        final ObjectArrayList<Annotation> result = new ObjectArrayList<Annotation>();
        final Object2ObjectMap<String, ObjectList<Annotation>> map = CompactAlignmentToAnnotationCountsMode.readAnnotations(annotReader);
        for (int referenceIndex = 0; referenceIndex < genome.size(); referenceIndex++) {
            final ObjectList<Annotation> list = map.get(genome.getReferenceName(referenceIndex));
            if (list != null) {
                Collections.sort(list, COMPARE_ANNOTATION_START);
                result.addAll(list);
            }
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
        return annotationIndex < annotations.length - 1;
    }

    public Annotation nextAnnotation() {
        return annotations[annotationIndex];
    }

    /**
     * * Method that populates the validAnnotationIndices
     * with annotations that overlap a given position
     *
     * @param refIndex
     * @param pos
     */
    public void advanceToPosition(final int refIndex, final int pos) {
        Annotation annNew = nextAnnotation();
        int annNewRefIndex = genome.getReferenceIndex(annNew.getChromosome());

        final String referenceName = genome.getReferenceName(refIndex);
        if (annNewRefIndex != refIndex) {
            advanceToChromosome(refIndex);
            annNew = nextAnnotation();
            annNewRefIndex = genome.getReferenceIndex(annNew.getChromosome());
        }

        // both the annotation and query position are on the same chromosome
        while (annNew.getStart() <= pos && annNewRefIndex == refIndex) {
            if (annNew.overlap(referenceName, pos)) {
                annotationIndicesInRange.add(annotationIndex);
            }
            if (hasMoreAnnotations()) {
                advanceToNextAnnotation();
                annNew = nextAnnotation();
                annNewRefIndex = genome.getReferenceIndex(annNew.getChromosome());
            } else {
                // there are no more annotations
                break;
            }
        }
        validAnnotationIndices.clear();
        validAnnotationIndices.addAll(annotationIndicesInRange);

        if (!annotationIndicesInRange.isEmpty()) {
            for (int index : annotationIndicesInRange) {
                // annotation does not overlap the current position
                if (!(annotations[index].overlap(referenceName, pos))) {
                    validAnnotationIndices.remove(index);
                    // check that the segments are not out of range
                    if (!annotations[index].withinRange(referenceName, pos)) {
                        annotationIndicesInRange.remove(index);
                    }
                }
            }
        }

    }


    /**
     * This method sets the annotationIndex to the index of the first occurrence
     * of annotations on this chromosome
     *
     * @param refIndex
     */
    public void advanceToChromosome(int refIndex) {

        while (hasMoreAnnotations()) {
            final Annotation ann = nextAnnotation();

            if (chromosome1StrictlyBefore2(genome.getReferenceIndex(ann.getChromosome()), refIndex)) {
                advanceToNextAnnotation();
            } else {
                if (chromosome1StrictlyBefore2(refIndex, genome.getReferenceIndex(ann.getChromosome()))) {
                    // we are past chromosome. Stop immediately and backup.
                    annotationIndex -= 1;
                    annotationIndex=Math.max(0,annotationIndex);
                    break;
                } else {
                    break;
                }
            }
        }
    }


    /**
     * Returns whether or not the given chromosome and position has overlapping annotations
     *
     * @param refIndex
     * @param pos
     * @return
     */
    public boolean hasOverlappingAnnotations(int refIndex, int pos) {
        advanceToPosition(refIndex, pos);
        return !(validAnnotationIndices.isEmpty());
    }

    private ObjectArrayList<Annotation> set = new ObjectArrayList<Annotation>();

    /**
     * Return the set of annotations that overlap with the current position.
     *
     * @return
     */
    public ObjectArrayList<Annotation> currentAnnotations() {
        set.clear();
        for (final int index : validAnnotationIndices) {
            set.add(annotations[index]);
        }
        return set;
    }

    public Annotation getAnnotation(int annoIndex) {
        return annotations[annoIndex];
    }

    /**
     * Determine if (currentChromosome,pos) is a position past the end of the given annotation
     */
    public boolean pastChosenAnnotation(int annoIndex, final String currentChromosome, final int pos) {
        final Annotation annotation = annotations[annoIndex];

        if (!annotation.getChromosome().equals(currentChromosome)) {
            // chromosome differ:
            if (chromosome1StrictlyBefore2(genome.getReferenceIndex(annotation.getChromosome()), genome.getReferenceIndex(currentChromosome))) {
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
            if (chromosome1StrictlyBefore2(genome.getReferenceIndex(annotation.getChromosome()), genome.getReferenceIndex(currentChromosome))) {
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
     * @param c1Index index of chromosome 1
     * @param c2Index index of chromosome 2
     * @return true when c1 occurs strictly before c2
     */
    private boolean chromosome1StrictlyBefore2(final int c1Index, final int c2Index) {
        return c1Index < c2Index;
    }


    public void advanceToNextAnnotation() {
        annotationIndex += 1;
    }

    public static final Comparator<? super Annotation> COMPARE_ANNOTATION_START = new Comparator<Annotation>() {
        @Override
        public int compare(final Annotation annotation, final Annotation annotation1) {
            return annotation.getStart() - annotation1.getStart();
        }
    };

    public String getAnnotationsLastChromosome() {
        int size = annotations.length;
        return annotations[size - 1].getChromosome();
    }

}
