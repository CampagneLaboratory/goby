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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Converts a compact alignment to plain text.
 *
 * @author Fabien Campagne
 */
public class AlignmentToTextMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-text";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts a compact alignment to text formats.";

    /**
     * The output file.
     */
    private String outputFilename;

    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;
    private AlignmentToTextIterateAlignments alignmentIterator;

    /**
     * The values to use for read lengths if none are found in the alignment entries/header.
     */
    private int defaultReadLength;

    /**
     * If header is written, used in PLAIN output (not SAM).
     */
    private boolean headerWritten = false;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    enum OutputFormat {
        PLAIN,
        SAM
    }

    private OutputFormat outputFormat;

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        final String[] inputFiles = jsapResult.getStringArray("input");
        basenames = AlignmentReaderImpl.getBasenames(inputFiles);
        outputFilename = jsapResult.getString("output");
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());
        if (outputFormat == OutputFormat.SAM) {
            // No output header for SAM format
            headerWritten = true;
        }
        defaultReadLength = jsapResult.getInt("constant-read-length");
        alignmentIterator = new AlignmentToTextIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);

        return this;
    }

    private class AlignmentToTextIterateAlignments extends IterateAlignments {
        private PrintStream outputStream;
        private OutputFormat outputFormat;
        private AlignmentReader cachedReader;
        private boolean hasReadIds;
        private DoubleIndexedIdentifier readIds;
        private int[] referenceLengths;

        public void setOutputWriter(final PrintStream outputStreawm, final OutputFormat outputFormat) {
            this.outputStream = outputStreawm;
            this.outputFormat = outputFormat;
        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader, final Alignments.AlignmentEntry alignmentEntry) {
            final int referenceIndex = alignmentEntry.getTargetIndex();

            if (cachedReader != alignmentReader) {
                hasReadIds = alignmentReader.getQueryIdentifiers().size() > 0;
                referenceLengths = alignmentReader.getTargetLength();
            }

            int startPosition = alignmentEntry.getPosition();
            final int alignmentLength = alignmentEntry.getQueryAlignedLength();

            if (!headerWritten) {
                printHeader(outputStream);
            }


            for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {
                final int queryIndex = alignmentEntry.getQueryIndex();

                // Get the length of the reference (if available)
                final int referenceLength;
                if (referenceLengths != null && ArrayUtils.getLength(referenceLengths) >= referenceIndex) {
                    referenceLength = referenceLengths[referenceIndex];
                } else {
                    referenceLength = -1;
                }
                switch (outputFormat) {
                    case PLAIN:
                        outputStream.printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%g\t%d\t%d\t%s\t%d%n",
                                hasReadIds ? readIds.getId(queryIndex) : queryIndex,
                                alignmentEntry.hasFragmentIndex() ? alignmentEntry.getFragmentIndex() : 0,
                                alignmentEntry.hasPairFlags() ? zeroPad(Integer.toBinaryString(alignmentEntry.getPairFlags()), 9) : 0,
                                alignmentEntry.hasPairAlignmentLink() ? alignmentEntry.getPairAlignmentLink().getFragmentIndex() : "",
                                alignmentEntry.hasPairAlignmentLink() ? getReferenceId(alignmentEntry.getPairAlignmentLink().getTargetIndex()) : "",
                                alignmentEntry.hasPairAlignmentLink() ? alignmentEntry.getPairAlignmentLink().getPosition() : "",
                                getReferenceId(alignmentEntry.getTargetIndex()),
                                referenceLength,
                                alignmentEntry.getNumberOfIndels(),
                                alignmentEntry.getNumberOfMismatches(),
                                alignmentEntry.getScore(),
                                startPosition,
                                alignmentLength,
                                alignmentEntry.getMatchingReverseStrand() ? "-" : "+",
                                alignmentEntry.hasMappingQuality() ? alignmentEntry.getMappingQuality() : 255);
                        break;
                    case SAM:
                        final int flag;
                        startPosition++;  // SAM is 1-based
                        if (alignmentEntry.hasPairFlags()) {
                            flag = alignmentEntry.getPairFlags();
                        } else {
                            // strand is encoded in 000010000 (binary), shift left by 4 bits.
                            flag = (alignmentEntry.getMatchingReverseStrand() ? 1 : 0) << 4;
                        }
                        final int mappingQuality;
                        if (alignmentEntry.hasMappingQuality()) {
                            mappingQuality = alignmentEntry.getMappingQuality();
                        } else {
                            mappingQuality = 255;
                        }
                        final String cigar = calculateCigar(alignmentEntry);
                        final int targetIndex = alignmentEntry.getTargetIndex();

                        String MRNM = "=";
                        int mPos = startPosition;
                        int inferredInsertSize = 0;
                        if (alignmentEntry.hasPairAlignmentLink()) {
                            final int pairTargetIndex = alignmentEntry.getPairAlignmentLink().getTargetIndex();
                            mPos = alignmentEntry.getPairAlignmentLink().getPosition() + 1;
                            if (pairTargetIndex != targetIndex) {
                                MRNM = getReferenceId(pairTargetIndex).toString();
                            }
                        }

                        final int readLength;
                        // check entry then header for the query/read length otherwise use default

                        readLength = alignmentEntry.getQueryLength();


                        final MutableString readSequence = getReadSequence(alignmentEntry, readLength);

                        final String readQualities = "........";  // TODO - should make this the correct length
                        outputStream.printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s%n",
                                hasReadIds ? readIds.getId(queryIndex) : queryIndex,
                                flag,
                                getReferenceId(targetIndex),
                                startPosition,
                                mappingQuality,
                                cigar,
                                MRNM,
                                mPos,
                                inferredInsertSize,
                                readSequence,
                                readQualities,
                                getTags(alignmentEntry, readLength));
                        break;
                }
            }
        }


    }

    private void printHeader(PrintStream outputStream) {
        headerWritten = true;
        outputStream.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n",
                "queryIndex",
                "queryFragmentIndex",
                "pairFlags",
                "pairFragmentIndex",
                "pairTarget",
                "pairPosition",
                "targetIdentifier",
                "referenceLength",
                "numIndels",
                "numMismatches",
                "score",
                "position",
                "alignmentLength",
                "strand",
                "mappingQuality");
    }

    private String zeroPad(final String val, final int length) {
        int addZeros = length - val.length();
        if (addZeros <= 0) {
            return val;
        }
        String format = String.format("%%0%dd%%s", addZeros);
        return String.format(format, 0, val);
    }

    private MutableString getReadSequence(final Alignments.AlignmentEntry alignmentEntry, final int readLength) {
        final MutableString sequence = new MutableString(readLength);
        if (readLength > 0) {
            for (int i = 0; i < readLength; ++i) {
                sequence.append('.');
            }
            for (final Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {
                final String to = var.getTo();
                if (var.getFrom().length() == to.length()) {
                    for (int i = 0; i < to.length(); i++) {
                        sequence.setCharAt(var.getReadIndex() - 1, to.charAt(i));
                    }
                }
            }
        }
        return sequence;
    }

    private String getTags(final Alignments.AlignmentEntry alignmentEntry, final int readLenth) {
        return String.format("NM:i:%d",
                alignmentEntry.getNumberOfMismatches());
        //  getMdAttribute(alignmentEntry, readLenth));
    }

    private String getMdAttribute(final Alignments.AlignmentEntry alignmentEntry, final int readLenth) {
        return "";
    }

    /* private String getMdAttribute(Alignments.AlignmentEntry alignmentEntry, int readLength) {
       MutableString mdAttribute = new MutableString();
       int previousReadIndex = 0;
       int readIndex = 0;
       int alreadyMatched = 0;
      boolean reverseStrand=alignmentEntry.getMatchingReverseStrand();
       for (Alignments.SequenceVariation var : alignmentEntry.getSequenceVariationsList()) {

           final String to = var.getTo();
           final String from = var.getFrom();
           readIndex = var.getReadIndex() - 1;
           alreadyMatched += to.length();
           if (from.length() > to.length()) {

               if (readIndex != previousReadIndex) {
                   mdAttribute.append(Integer.toString(readIndex - previousReadIndex));
                   previousReadIndex = readIndex;
               }
               mdAttribute.append("^" + from);


           } else if (from.length() < to.length()) {

               if (readIndex != previousReadIndex) {
                   mdAttribute.append(Integer.toString(readIndex - previousReadIndex));
                   previousReadIndex = readIndex;
               }
               mdAttribute.append(to);


           } else {
               // point mutation:
               mdAttribute.append(Integer.toString(to.length()));
               mdAttribute.append(to);
               previousReadIndex = readIndex;
           }

       }

       mdAttribute.append(Integer.toString(readLength - alreadyMatched));


       return mdAttribute.toString();
   }
    */
    private String calculateCigar(final Alignments.AlignmentEntry alignmentEntry) {
        return (alignmentEntry.getQueryAlignedLength() - alignmentEntry.getNumberOfIndels()) + "M";
    }

    /**
     * Display the alignments as text files.
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
                case PLAIN:
                    printHeader(stream);
                    break;
                case SAM:
                    stream.printf("@HD\tVN:1.0%n" + "@PG\tGoby\tVN:"
                            + VersionUtils.getImplementationVersion(GobyDriver.class) + "%n");

                    for (final String basename : basenames) {
                        final AlignmentReaderImpl reader = new AlignmentReaderImpl(basename);
                        reader.readHeader();
                        final IndexedIdentifier identifiers = reader.getTargetIdentifiers();
                        for (final MutableString targetId : identifiers.keySet()) {
                            if (targetId != null) {
                                final int[] targetLengths = reader.getTargetLength();
                                if (targetLengths != null) {
                                    stream.printf("@SQ\tSN:%s\tLN:%d%n", targetId,
                                            targetLengths[identifiers.getInt(targetId)]);
                                } else {
                                    stream.printf("@SQ\tSN:%s%n", targetId);
                                }
                            }
                        }
                    }
                    break;
            }

            alignmentIterator.setOutputWriter(stream, outputFormat);
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
     * @throws JSAPException error parsing
     * @throws IOException   error parsing or executing.
     */

    public static void main(final String[] args) throws JSAPException, IOException {
        new AlignmentToTextMode().configure(args).execute();
    }
}
