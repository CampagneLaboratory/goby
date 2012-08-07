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
import edu.cornell.med.icb.goby.algorithmic.data.Interval;
import edu.cornell.med.icb.goby.modes.CompactAlignmentToAnnotationCountsMode;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.Object2ObjectArrayMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.lang.MutableString;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

/**
 * @author Fabien Campagne
 *         Date: 7/27/12
 *         Time: 11:24 AM
 */
public class RandomAccessAnnotations {

    private Comparator<Interval> compareIntervals = new Comparator<Interval>() {
        @Override
        public int compare(final Interval a, final Interval b) {
            return a.start - b.start;
        }
    };

    /**
     * Load and sort annotations by their end position.
     *
     * @param filename
     * @throws java.io.IOException
     */
    public void loadAnnotations(final String filename) throws IOException {
        loadAnnotations(new FileReader(filename));
    }

    IndexedIdentifier references = new IndexedIdentifier();

    /**
     * Load and sort annotations by their start position.
     *
     * @param annotReader
     * @throws java.io.IOException
     */
    public void loadAnnotations(final Reader annotReader) throws IOException {

        final ObjectArrayList<Interval> result = new ObjectArrayList<Interval>();
        final Object2ObjectMap<String, ObjectList<Annotation>> map = CompactAlignmentToAnnotationCountsMode.readAnnotations(annotReader);
        for (final String key : map.keySet()) {
            final ObjectList<Annotation> list = map.get(key);

            if (list != null) {
                Collections.sort(list, compareAnnotationStart);
                for (final Annotation element : list) {
                    final Interval interval = new Interval();
                    interval.referenceIndex = references.registerIdentifier(new MutableString(element.getChromosome()));
                    interval.start = element.getStart();
                    interval.end = element.getEnd();
                    interval.id = element.getId();
                    result.add(interval);
                }
                final Interval[] intervals = new Interval[list.size()];

                chromosomeToMap.put(key, result.toArray(intervals));
            }
        }

    }

    Object2ObjectMap<String, Interval[]> chromosomeToMap = new Object2ObjectArrayMap<String, Interval[]>();
    private final Comparator<? super Annotation> compareAnnotationStart = new Comparator<Annotation>() {
        @Override
        public int compare(final Annotation annotation, final Annotation annotation1) {
            return annotation.getStart() - annotation1.getStart();
        }
    };

    protected void addAnnotation(String elementId, String chr, int start, int end) {
        Interval newInterval = new Interval();
        newInterval.id = elementId;
        newInterval.referenceIndex = references.registerIdentifier(new MutableString(chr));
        newInterval.start = start;
        newInterval.end = end;
        Interval[] array = chromosomeToMap.containsKey(chr) ? chromosomeToMap.get(chr) : new Interval[0];
        Interval[] dest = new Interval[array.length + 1];
        System.arraycopy(array, 0, dest, 0, array.length);
        dest[array.length] = newInterval;
        chromosomeToMap.put(chr, dest);


    }

    /**
     * Find the annotation that overlaps an interval.
     *
     * @param chromosome
     * @param start
     * @param end
     * @return The annotation or null is none overlaps the query interval.
     */
    public Interval find(final String chromosome, final int start, final int end) {
        final Interval[] intervals = chromosomeToMap.get(chromosome);
        if (intervals==null) {
            return null;
        }
        singleton.start = start;
        final int index = Arrays.binarySearch(intervals, singleton, compareIntervals);
        final int ip;
        if (index < 0) {

            // position does not match min exactly, we get the insertion point instead:
            ip = -index - 1;
            return checkOverlap(intervals, ip, start, end);
        } else {
            ip = index;
            return checkOverlap(intervals, ip + 1, start, end);
        }

    }

    private Interval checkOverlap(final Interval[] intervals, final int insertionPoint, final int start, final int end) {
        if (insertionPoint <= 0) {
            return null;
        }
        final Interval previous = intervals[insertionPoint - 1];
        if (previous.start <= start && previous.end >= end) {
            return previous;
        }

        return null;

    }

    private Interval singleton = new Interval();
}
