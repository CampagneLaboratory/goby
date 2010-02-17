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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;

import java.io.Closeable;
import java.io.Reader;
import java.util.NoSuchElementException;

/**
 * Parses the Last MAF format.
 * As per http://genome.ucsc.edu/FAQ/FAQformat#format5
 * ----------------------------------------------------------------------------
 *  Quoting this documentation
 * ----------------------------------------------------------------------------
 * The multiple alignment format stores a series of multiple alignments in a
 * format that is easy to parse and relatively easy to read. This format stores
 * multiple alignments at the DNA level between entire genomes. Previously used
 * formats are suitable for multiple alignments of single proteins or regions
 * of DNA without rearrangements, but would require considerable extension to
 * cope with genomic issues such as forward and reverse strand directions,
 * multiple pieces to the alignment, and so forth.
 *
 * General Structure
 *
 * The .maf format is line-oriented. Each multiple alignment ends with a blank
 * line. Each sequence in an alignment is on a single line, which can get quite
 * long, but there is no length limit. Words in a line are delimited by any
 * white space. Lines starting with # are considered to be comments. Lines
 * starting with ## can be ignored by most programs, but contain meta-data of
 * one form or another.
 *
 * The file is divided into paragraphs that terminate in a blank line. Within a
 * paragraph, the first word of a line indicates its type. Each multiple
 * alignment is in a separate paragraph that begins with an "a" line and
 * contains an "s" line for each sequence in the multiple alignment. Some MAF
 * files may contain other optional line types:
 *
 *      an "i" line containing information about what is in the aligned species
 *      DNA before and after the immediately preceding "s" line
 *
 *      an "e" line containing information about the size of the gap between
 *      the alignments that span the current block
 *
 *      a "q" line indicating the quality of each aligned base for the species
 *
 * Parsers may ignore any other types of paragraphs and other types of lines
 * within an alignment paragraph."
 * ----------------------------------------------------------------------------
 */
public class LastParser implements Closeable {
    private LineIterator iterator;
    private boolean entryFound;
    private float score;
    private final ObjectArrayList<AlignedSequence> alignedSequences;
    private final ObjectArraySet<AlignedSequence> poolOfAlignedSequences;

    public LastParser(final Reader reader) {
        iterator = new LineIterator(new FastBufferedReader(reader));
        alignedSequences = new ObjectArrayList<AlignedSequence>();
        poolOfAlignedSequences = new ObjectArraySet<AlignedSequence>();
    }

    /**
     * Returns true when the MAF file has another entry.
     *
     * @return True or False.
     */
    public boolean hasNext() {
        if (entryFound) {
            return true;
        } else {
            readScore();
            readMatches();
            return entryFound;
        }
    }

    /**
     * Reads MAF file input lines until ALL aligned sequences have been parsed, one per line.
     * Currently only reading lines beginning with "s".
     *
     * According to the MAF format, "e", "i" and "q" lines may also populate the multiple-
     * alignment block.
     */
    private void readMatches() {
        // N.B. MAF entries can have more than 2 lines ... stop reading at blank line
        while (iterator.hasNext()) {
            final MutableString line = iterator.next();
            if (line.startsWith("#")) {
                continue;
            }
            /* TODO (suggestion) // continue reading the alignment block by ignoring these lines
            if (line.startsWith("e") || line.startsWith("i") || line.startsWith("q")) {
                continue;
            }
            */
            if (line.startsWith("s")) {
                final String[] tokens = line.toString().split("[\\ ]+");
                final AlignedSequence seq;

                seq = getAvailableAlignedSequence();
                seq.sequenceIdentifier.setLength(0);
                seq.sequenceIdentifier.append(tokens[1]);
                seq.alignedStart = Integer.parseInt(tokens[2]);
                seq.alignedLength = Integer.parseInt(tokens[3]);
                seq.strand = tokens[4].charAt(0);
                seq.sequenceLength = Integer.parseInt(tokens[5]);
                seq.alignment.setLength(0);
                seq.alignment.append(tokens[6]);

                alignedSequences.add(seq);
            } else {
                // An alignment is terminated by a blank line.
                break;
            }
        }
    }

    private void readScore() {
        MutableString line;
        while (iterator.hasNext()) {
            line = iterator.next();
            if (line.startsWith("#")) {
                continue;
            }
            if (line.startsWith("a")) {
                final String[] tokens = line.toString().split("[\\ ]+");

                this.score = Float.parseFloat(tokens[1].substring(6));
                //return alignedSequences to the pool so that we can reuse them for the next parsed entry:
                poolOfAlignedSequences.addAll(alignedSequences);
                alignedSequences.clear();
                entryFound = true;
                return;
            }
        }
    }

    /**
     * Returns the score of the alignment.
     *
     * @return
     */
    public float getScore() {
        return score;
    }

    public void next() {
        if (!hasNext()) {
            throw new NoSuchElementException("next() cannot be called when hasNext returns false.");
        }
        // do nothing, we parsed the entry already.
        entryFound = false;
    }

    /**
     * Returns the list of aligned sequences associated with this entry.
     *
     * @return
     */
    public ObjectArrayList<AlignedSequence> getAlignedSequences() {
        return alignedSequences;
    }

    public synchronized AlignedSequence getAvailableAlignedSequence() {
        final AlignedSequence seq;
        if (poolOfAlignedSequences.isEmpty()) {
            seq = new AlignedSequence();
            return seq;
        } else {
            seq = poolOfAlignedSequences.iterator().next();
            poolOfAlignedSequences.remove(seq);
            return seq;
        }
    }

    /**
     * {@inheritDoc}
     */
    public void close() {
        iterator = null;
    }
}
