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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Display the sequence variations found in alignments.
 *
 * @author Fabien Campagne
 */
public class DisplaySequenceVariationsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "display-sequence-variations";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Display the sequence variations found in an alignment";

    /**
     * The input filenames.
     */
    private String[] inputFilenames;

    /**
     * The output file.
     */
    private String outputFilename;
    /**
     * The input basenames.
     */
    private String[] basenames;
    private MyIterateAlignments alignmentIterator;
    private FirstPassIterateAlignments firstPassIterator;
    private boolean thresholds;
    private int minimumUniqueReadIndices = 1;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    enum OutputFormat {
        CONCISE,
        TSV,
        TAB_DELIMITED,
        TAB_SINGLE_BASE,
    }

    private OutputFormat outputFormat;

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilenames = jsapResult.getStringArray("input");
        basenames = AlignmentReaderImpl.getBasenames(inputFilenames);
        outputFilename = jsapResult.getString("output");
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());
        minimumUniqueReadIndices = jsapResult.getInt("minimum-read-indices");

        alignmentIterator = new MyIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);

        thresholds = minimumUniqueReadIndices > 0;
        if (thresholds) {
            firstPassIterator = new FirstPassIterateAlignments();
            firstPassIterator.parseIncludeReferenceArgument(jsapResult);
        }
        return this;
    }

    /**
     * Display sequence variations.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintStream stream = null;
        try {
            stream = outputFilename == null ? System.out
                    : new PrintStream(new FileOutputStream(outputFilename));
            switch (outputFormat) {
                case CONCISE:
                    break;
                case TAB_DELIMITED:
                case TAB_SINGLE_BASE:
                case TSV:
                    stream.println("basename\tquery-index\ttarget-id\tposition-on-reference\tread-index\tvar-from\tvar-to\ttype");
                    break;
            }

            if (thresholds) {
                firstPassIterator.iterate(basenames);
            }

            alignmentIterator.setOutputWriter(stream, outputFormat);
            if (thresholds) {
                alignmentIterator.setFirstPass(firstPassIterator);
            }
            // Iterate through each alignment and write sequence variations to output file:
            alignmentIterator.iterate(basenames);
        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
        }
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new DisplaySequenceVariationsMode().configure(args).execute();
    }

    /**
     * Collect the list of read indices where each variation is observed, for a given
     * reference position.
     */
    private static class FirstPassIterateAlignments extends IterateAlignments {
        // reference index -> reference Position -> readIndex list
        private final Int2ObjectOpenHashMap<Int2ObjectOpenHashMap<IntArraySet>> readIndicesForReferencePositions =
                new Int2ObjectOpenHashMap<Int2ObjectOpenHashMap<IntArraySet>>();

        public Int2ObjectOpenHashMap<Int2ObjectOpenHashMap<IntArraySet>> getReadIndicesForReferencePositions() {
            return readIndicesForReferencePositions;
        }

        private FirstPassIterateAlignments() {
            super();
            readIndicesForReferencePositions.defaultReturnValue(new Int2ObjectOpenHashMap<IntArraySet>());
        }

        public IntArraySet getReadIndices(final int referenceIndex, final int referencePosition) {
            final Int2ObjectOpenHashMap<IntArraySet> referencePositionsMap = readIndicesForReferencePositions.get(referenceIndex);
            return referencePositionsMap.get(referencePosition);
        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader, final Alignments.AlignmentEntry alignmentEntry) {
            final int referenceIndex = alignmentEntry.getTargetIndex();
            final int alignmentPositionOnReference = alignmentEntry.getPosition();
            final Int2ObjectOpenHashMap<IntArraySet> referencePositionsMap = readIndicesForReferencePositions.get(referenceIndex);
            for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                final int referencePosition = var.getPosition() + alignmentPositionOnReference;
                IntArraySet readIndexList = referencePositionsMap.get(referencePosition);
                if (readIndexList == null) {
                    readIndexList = new IntArraySet();
                }
                readIndexList.add(var.getReadIndex());
                referencePositionsMap.put(referencePosition, readIndexList);
            }

        }
    }

    private class MyIterateAlignments extends IterateAlignments {
        private PrintStream outputStream;
        private OutputFormat outputFormat;
        private FirstPassIterateAlignments firstPassIterator;

        public void setOutputWriter(final PrintStream stream, final OutputFormat format) {
            this.outputStream = stream;
            this.outputFormat = format;
        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader, final Alignments.AlignmentEntry alignmentEntry) {
            String basename = alignmentReader.basename();
            // remove the path:
            basename = FilenameUtils.getBaseName(basename);
            switch (outputFormat) {

                case CONCISE: {
                    if (alignmentEntry.getSequenceVariationsCount() > 0) {
                        outputStream.print(String.format("%d %s ",

                                alignmentEntry.getQueryIndex(),
                                getReferenceId(alignmentEntry.getTargetIndex())));
                        boolean variations = false;
                        for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                            // convert variation position to position on the reference:
                            final int positionOnReference = alignmentEntry.getPosition() + var.getPosition();
                            final int referenceIndex = alignmentEntry.getTargetIndex();
                            boolean keepVar = true;
                            keepVar = determineKeepVariation(positionOnReference, referenceIndex, keepVar);
                            if (keepVar) {
                                variations = true;
                                outputStream.print(String.format("%d:%d:%s/%s,",


                                        positionOnReference,
                                        var.getReadIndex(),
                                        var.getFrom(),
                                        var.getTo()));
                            }
                        }
                        if (variations) {
                            outputStream.println();
                        }
                    }
                }
                break;
                case TSV:
                case TAB_DELIMITED: {
                    final boolean variations = false;

                    for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {

                        // convert variation position to position on the reference:
                        final int positionOnReference = alignmentEntry.getPosition() + var.getPosition();
                        final int readIndex = var.getReadIndex();
                        final String from = var.getFrom();
                        final String to = var.getTo();
                        final int referenceIndex = alignmentEntry.getTargetIndex();
                        final byte[] qualityScores;
                        if (var.hasToQuality()) {
                            System.out.println("Score");
                            qualityScores = var.getToQuality().toByteArray();
                        } else {
                            System.out.println("No score");
                            qualityScores = null;
                        }

                        boolean keepVar = true;
                        keepVar = determineKeepVariation(positionOnReference, referenceIndex, keepVar);
                        if (keepVar && !isAllNs(to)) {
                            printTab(alignmentEntry, basename, positionOnReference, readIndex, from, to, qualityScores);
                        }
                    }


                }
                break;
                case TAB_SINGLE_BASE: {
                    for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {

                        // convert variation position to position on the reference:
                        final int positionOnReference = alignmentEntry.getPosition() + var.getPosition();
                        final int readIndex = var.getReadIndex();
                        final String from = var.getFrom();
                        final String to = var.getTo();
                        final int fromLength = from.length();
                        final int toLength = to.length();
                        final int referenceIndex = alignmentEntry.getTargetIndex();
                        final byte[] qualityScores;
                        if (var.hasToQuality()) {
                            qualityScores = var.getToQuality().toByteArray();
                        } else {
                            qualityScores = null;
                        }
                        boolean keepVar = true;
                        keepVar = determineKeepVariation(positionOnReference, referenceIndex, keepVar);
                        if (keepVar && !isAllNs(to)) {
                            final int maxLength = Math.max(fromLength, toLength);
                            int fromOffset = 0;
                            int toOffset = 0;
                            int readIndexIncrementValue = (alignmentEntry.getMatchingReverseStrand() ? -1 : 1);                            
                            byte[] toScore;
                            if (qualityScores != null) {
                                toScore = new byte[1];
                            } else {
                                toScore = null;
                            }
                            for (int i = 0; i < maxLength; i++) {
                                final char fromChar = var.getFrom().charAt(i);
                                final char toChar = var.getTo().charAt(i);
                                if (qualityScores != null) {
                                    if (i < qualityScores.length) {
                                        toScore[0] = qualityScores[i];
                                    } else {
                                        toScore = null;
                                    }
                                }
                                printTab(alignmentEntry, basename,
                                        positionOnReference + fromOffset,
                                        readIndex + toOffset,
                                        i < fromLength ? Character.toString(fromChar) : "",
                                        i < toLength ? Character.toString(toChar) : "",
                                        toScore);
                                if (fromChar != '-') {
                                    fromOffset += 1;
                                }
                                if (toChar != '-') {
                                    toOffset += readIndexIncrementValue;
                                }
                            }
                        }
                    }
                }
                break;
            }
        }

        private boolean determineKeepVariation(final int positionOnReference, final int referenceIndex, boolean keepVar) {
            if (thresholds) {
                final IntArraySet indices = firstPassIterator.getReadIndices(referenceIndex, positionOnReference);
                if (indices != null && indices.size() < minimumUniqueReadIndices) {
                    keepVar = false;
                }
            }
            return keepVar;
        }

        private boolean isAllNs(final String to) {
            for (int i = 0; i < to.length(); ++i) {
                if (to.charAt(i) != 'N') {
                    return false;
                }

            }
            return true;
        }

        private void printTab(final Alignments.AlignmentEntry alignmentEntry, final String basename, final int positionOnReference, final int readIndex, final String from, final String to, final byte[] qualityScores) {
            final String type;
            if (from.contains("-") || from.length() == 0) {
                // insertion in read sequence.
                type = "READ_INSERTION";
            } else if (to.contains("-") || to.length() == 0) {
                // deletion in read sequence.
                type = "READ_DELETION";
            } else {
                // one or more bases are mutated. no insertions or deletions.
                type = "MUTATION";
            }
            StringBuffer qualStr = new StringBuffer();
            if (qualityScores != null) {
                for (int i = 0; i < qualityScores.length; i++) {
                    if (i > 0) {
                        qualStr.append(",");
                    }
                    qualStr.append(qualityScores[i]);
                }
            }
            outputStream.println(String.format("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s",
                    basename,
                    alignmentEntry.getQueryIndex(),
                    getReferenceId(alignmentEntry.getTargetIndex()),
                    positionOnReference,
                    readIndex,
                    from,
                    to,
                    type,
                    qualStr.toString()));
        }

        public void setFirstPass(final FirstPassIterateAlignments firstPassIterator) {
            this.firstPassIterator = firstPassIterator;
        }
    }
}
