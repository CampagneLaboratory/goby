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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.IterateAlignments;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

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
    private MyIterateAlignments alignmentIterator;
    private int defaultReadLength;


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
        basenames = AlignmentReader.getBasenames(inputFiles);
        outputFilename = jsapResult.getString("output");
        outputFormat = OutputFormat.valueOf(jsapResult.getString("format").toUpperCase());
        defaultReadLength = jsapResult.getInt("constant-read-length");
        alignmentIterator = new MyIterateAlignments();
        alignmentIterator.parseIncludeReferenceArgument(jsapResult);

        return this;
    }

    private class MyIterateAlignments extends IterateAlignments {
        private PrintWriter outputWriter;
        private OutputFormat outputFormat;
        private AlignmentReader cachedReader;
        private boolean hasReadIds;
        private DoubleIndexedIdentifier readIds;
        private int[] referenceLengths;

        public void setOutputWriter(final PrintWriter outputWriter, final OutputFormat outputFormat) {
            this.outputWriter = outputWriter;
            this.outputFormat = outputFormat;
        }

        @Override
        public void processAlignmentEntry(final AlignmentReader alignmentReader, final Alignments.AlignmentEntry alignmentEntry) {
            final int referenceIndex = alignmentEntry.getTargetIndex();

            if (cachedReader != alignmentReader) {
                hasReadIds = alignmentReader.getQueryIdentifiers().size() > 0;
                final DoubleIndexedIdentifier readIds =
                        new DoubleIndexedIdentifier(alignmentReader.getQueryIdentifiers());
                referenceLengths = alignmentReader.getTargetLength();
            }

            final int startPosition = alignmentEntry.getPosition();
            final int alignmentLength = alignmentEntry.getQueryAlignedLength();
            for (int i = 0; i < alignmentEntry.getMultiplicity(); ++i) {
                final int queryIndex = alignmentEntry.getQueryIndex();

                // Get the length of the reference (if available)
                final int referenceLength;
                if (ArrayUtils.getLength(referenceLengths) >= referenceIndex) {
                    referenceLength = referenceLengths[referenceIndex];
                } else {
                    referenceLength = -1;
                }
                switch (outputFormat) {
                    case PLAIN:
                        outputWriter.write(String.format("%s\t%s\t%d\t%d\t%d\t%g\t%d\t%d\t%b%n",
                                hasReadIds ? readIds.getId(queryIndex) : queryIndex,
                                getReferenceId(alignmentEntry.getTargetIndex()),
                                referenceLength,
                                alignmentEntry.getNumberOfIndels(),
                                alignmentEntry.getNumberOfMismatches(),
                                alignmentEntry.getScore(),
                                startPosition,
                                alignmentLength,
                                alignmentEntry.getMatchingReverseStrand()));
                        break;
                    case SAM:
                        final int flag = (alignmentEntry.getMatchingReverseStrand() ? 1 : 0) << 4;   // strand is encoded in 0x10, shift left by 4 bits.
                        final int mappingQuality = 255;
                        final String cigar = calculateCigar(alignmentEntry);
                        final String MRNM = "=";
                        final int inferredInsertSize = 0;
                        int readLength = defaultReadLength;
                        if (alignmentReader.hasQueryLengths()) {
                            readLength = alignmentReader.getQueryLength(alignmentEntry.getQueryIndex());
                        }

                        final MutableString readSequence = getReadSequence(alignmentEntry, readLength);

                        final String readQualities = "........";
                        outputWriter.write(String.format("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s%n",
                                hasReadIds ? readIds.getId(queryIndex) : queryIndex,
                                flag,
                                getReferenceId(alignmentEntry.getTargetIndex()),
                                startPosition,
                                mappingQuality,
                                cigar,
                                MRNM,
                                startPosition,
                                inferredInsertSize,
                                readSequence,
                                readQualities,
                                getTags(alignmentEntry, readLength)
                        ));
                        break;
                }
            }
        }
    }

    private MutableString getReadSequence(final Alignments.AlignmentEntry alignmentEntry, final int readLength) {
        final MutableString sequence = new MutableString();
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
        final PrintWriter writer = outputFilename == null ? new PrintWriter(System.out) :
                new PrintWriter(new FileWriter(outputFilename));
        switch (outputFormat) {
            case PLAIN:
                writer.printf("queryId\treferenceId\treferenceLength\tnumberOfIndels\tnumberOfMismatches\tscore\tstartPosition\talignmentLength\tmatchesReverseStrand%n");
            case SAM:
                writer.printf("@HD\tVN:1.0%n" +
                        "@PG\tGoby\tVN:" + VersionUtils.getImplementationVersion(GobyDriver.class) + "%n");

                for (final String basename : basenames) {
                    final AlignmentReader reader = new AlignmentReader(basename);
                    reader.readHeader();
                    final IndexedIdentifier identifiers = reader.getTargetIdentifiers();
                    for (final MutableString targetId : identifiers.keySet()) {
                        if (targetId != null) {
                            final int[] targetLengths = reader.getTargetLength();
                            if (targetLengths != null) {
                                writer.printf("@SQ\tSN:%s\tLN:%d%n", targetId, targetLengths[identifiers.getInt(targetId)]);
                            } else {
                                writer.printf("@SQ\tSN:%s%n", targetId);
                            }
                        }
                    }
                }
                break;
        }

        try {
            alignmentIterator.setOutputWriter(writer, outputFormat);
            // Iterate through each alignment and write sequence variations to output file:
            alignmentIterator.iterate(basenames);
        } finally {
            IOUtils.closeQuietly(writer);
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
