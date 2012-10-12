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

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.longs.Long2IntArrayMap;
import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;

import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.Collection;
import java.util.Collections;

/**
 * A mode that extracts splicing events (instances when a read cross a splice junction) from alignment files.
 * @author Fabien Campagne
 *         Date: 9/14/12
 *         Time: 10:54 AM
 */
public class ExportSplicingEvents {
    private DoubleIndexedIdentifier reverseIds = null;


    final Object2IntArrayMap<SpliceEvent> eventMap = new Object2IntArrayMap<SpliceEvent>();
    int lastTargetIndex = -1;
    private String sampleFilename = "sample-id";
    final PrintWriter output;
    /**
     * The minimum mapping quality that alignments must have to be considered in a tally of the splicing events. Spliced
     * reads with a mapping quality less than this threshold are ignored. The default value is 255.
     */
    private int qualTheshold=255;

    public ExportSplicingEvents(final Writer output) {
        this.output = new PrintWriter(output);
        eventMap.defaultReturnValue(0);
    }

    public void process(final Collection<Alignments.AlignmentEntry> entries,
                        final DoubleIndexedIdentifier reverseIds) throws IOException {
        this.reverseIds = reverseIds;
        for (final Alignments.AlignmentEntry entry : entries) {

            processEntry(eventMap, entry);
        }
        flushPreviousRef(null, eventMap);
        output.flush();
    }

    public void process(final String filename) throws IOException {
        final AlignmentReader reader = new AlignmentReaderImpl(filename);
        reader.readHeader();
        final IndexedIdentifier targetIds = reader.getTargetIdentifiers();
        reverseIds = new DoubleIndexedIdentifier(targetIds);
        sampleFilename = FilenameUtils.getBaseName(filename);
        lastTargetIndex = -1;

        for (final Alignments.AlignmentEntry entry : reader) {

            processEntry(eventMap, entry);
        }
        flushPreviousRef(null, eventMap);
        output.flush();
    }

    private void processEntry(final Object2IntArrayMap<SpliceEvent> eventMap, final Alignments.AlignmentEntry entry) {
        assert reverseIds != null : " reverse Id cannot be null";
        if ( entry.hasMappingQuality() && entry.getMappingQuality()<qualTheshold) {
            // ignore alignments that have a mapping quality field with a value strictly less than the threshold.
            return;
        }
        if (entry.hasSplicedForwardAlignmentLink()) {
            flushPreviousRef(entry, eventMap);

            final Alignments.RelatedAlignmentEntry link = entry.getSplicedForwardAlignmentLink();
            final int firstBase = entry.getPosition() + entry.getTargetAlignedLength() + 1;
            if (link.getTargetIndex() != entry.getTargetIndex()) {
                // this is not a splicing event, but a fusion, ignore.
            } else {
                final int endBase = link.getPosition()+1;
                // increment the count of the observed splice junction:
                final char strand = entry.getMatchingReverseStrand() ? '-' : '+';
                observeSplice(firstBase, endBase, strand);
            }
        }
    }

    private void observeSplice(final int firstBase, final int endBase, final char strand) {
        final SpliceEvent key = new SpliceEvent(firstBase, endBase, strand);
        final int count = eventMap.getInt(key);

        eventMap.put(key, count + 1);
    }

    public void setMinMappingQuality(int qualThreshold) {
        this.qualTheshold=qualThreshold;
    }


    private class SpliceEvent implements Comparable<SpliceEvent>{
        int firstBase;
        int lastBase;
        char strand;

        public SpliceEvent(final int firstBase, final int endBase, final char strand) {
            this.firstBase = firstBase;
            this.lastBase = endBase;
            this.strand = strand;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            SpliceEvent that = (SpliceEvent) o;

            if (firstBase != that.firstBase) return false;
            if (lastBase != that.lastBase) return false;
            if (strand != that.strand) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = firstBase;
            result = 31 * result + lastBase;
            result = 31 * result + (int) strand;
            return result;
        }

        @Override
        public int compareTo(SpliceEvent spliceEvent) {
         return  firstBase-spliceEvent.firstBase;
        }
    }


    private void flushPreviousRef(final Alignments.AlignmentEntry entry, final Object2IntArrayMap<SpliceEvent> eventMap) {
        if (lastTargetIndex == -1 && entry!=null) {
            lastTargetIndex=entry.getTargetIndex();
            return;
        }
        if (entry == null || entry.getTargetIndex() != lastTargetIndex) {

            final ObjectSet<SpliceEvent> spliceEvents = eventMap.keySet();
            ObjectArrayList<SpliceEvent> list= new ObjectArrayList<SpliceEvent>();
            list.addAll(spliceEvents);
            Collections.sort(list);
            for (final SpliceEvent event : list) {

                final int firstBase = event.firstBase;
                final int lastBase = event.lastBase;
                output.print(sampleFilename);
                output.print("\t");
                output.print(reverseIds.getId(lastTargetIndex));
                output.print("\t");
                output.print(firstBase);
                output.print("\t");
                output.print(lastBase);
                output.print("\t");
                output.print(event.strand);
                output.print("\t");
                output.print("??"); // motif
                output.print("\t");
                output.print(eventMap.getInt(event)); // the number of times the splicing event was observed in the sample
                output.print("\t");
                output.print("0.0");
                output.println();

            }
            eventMap.clear();
            if (entry != null) {
                lastTargetIndex = entry.getTargetIndex();
            }

        }

    }

}



