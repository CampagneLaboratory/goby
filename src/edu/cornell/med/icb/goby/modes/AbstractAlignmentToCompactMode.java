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

import com.google.protobuf.ByteString;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.alignments.filters.AlignmentQualityFilter;
import edu.cornell.med.icb.goby.alignments.filters.PercentMismatchesQualityFilter;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Abstract class for modes that convert alignment formats to compact format.
 * The JSAP file of the concrete implementation must include the following options in addition to the options specific
 * to the concrete mode.
 * <p/>
 * "input"
 * "output"
 * "query-reads-ids"
 * "target-reference-ids"
 * "propagate-query-ids"
 * "propagate-target-ids"
 * "read-index-filter"
 * "ambiguity-threshold"
 * "quality-filter-parameters"
 * <p/>
 * <p/>
 *
 * @author Stuart Andrews
 *         Date: Sep 11, 2009
 *         Time: 9:30:16 PM
 */
public abstract class AbstractAlignmentToCompactMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AbstractAlignmentToCompactMode.class);

    /**
     * default ambiguity threshold.
     */
    protected static final int DEFAULT_M_PARAM = 2;

    /**
     * Input file.
     */
    protected String inputFile;

    /**
     * Output file.
     */
    protected String outputFile;

    /**
     * Query / target identifiers.
     */
    protected String queryReadIdsFilename;
    protected String targetReferenceIdsFilename;
    protected boolean propagateQueryIds;
    protected boolean propagateTargetIds;

    /**
     * Conversion parameters.
     */
    protected String qualityFilterParameters = "";
    protected AlignmentQualityFilter qualityFilter;
    protected File readIndexFilterFile;
    protected int mParameter = DEFAULT_M_PARAM;
    protected int numberOfReads = -1;

    /**
     * Identifiers for the query sequences in this alignment.
     */
    protected final IndexedIdentifier queryIds = new IndexedIdentifier();
    /**
     * Identifiers for the target sequences in this alignment.
     */
    protected final IndexedIdentifier targetIds = new IndexedIdentifier();
    protected int numberOfReadsFromCommandLine;

    /**
     * Indicate that the file being imported is from a third party. This means that:
     * <ol>
     * <li>queryNames are not integers, but are strings that need to be converted to indices.</li>
     * <li>targetNames should be treated as strings and defined from the input.</li>
     * <li>targetLength information should be read from the input, not from a supplied target
     * file in compact format.</li>
     * </ol>
     * <p/>
     * False by default when constructed, overidden by configure with default
     * configuration=true when run as a mode on the command line, set to false
     * explictly each time another Goby mode needs to import internally the result
     * of a Goby search.
     */
    protected boolean thirdPartyInput = true;
    protected int smallestQueryIndex;
    protected int largestQueryIndex = -1;

    /**
     * This method is deprecated, store read lengths directly into the alignment entry instead.
     *
     * @return
     */
    @Deprecated
    protected int[] createReadLengthArray() {
        return new int[largestQueryIndex - smallestQueryIndex + 1];
    }

    /**
     * Scan.
     *
     * @param readIndexFilter
     * @param writer
     * @param targetIds
     * @param tmhWriter
     * @return number of alignment entries written
     * @throws IOException error parsing
     */
    protected abstract int scan(ReadSet readIndexFilter, IndexedIdentifier targetIds,
                                AlignmentWriter writer,
                                AlignmentTooManyHitsWriter tmhWriter) throws IOException;


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
        //
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFile = jsapResult.getString("input");
        outputFile = jsapResult.getString("output");

        queryReadIdsFilename = jsapResult.getString("query-reads-ids");
        targetReferenceIdsFilename = jsapResult.getString("target-reference-ids");
        propagateQueryIds = jsapResult.getBoolean("propagate-query-ids");
        propagateTargetIds = jsapResult.getBoolean("propagate-target-ids");

        readIndexFilterFile = jsapResult.getFile("read-index-filter");
        mParameter = jsapResult.getInt("ambiguity-threshold");
        qualityFilterParameters = jsapResult.getString("quality-filter-parameters");
        thirdPartyInput = jsapResult.getBoolean("third-party-input");

        return this;
    }

    /**
     * Run the alignment-to-compact mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        // read target/query identifier lookup table, and initialize output alignment
        // file with this information
        final TransferIds transferIds = new TransferIds().invoke();
        final ReadSet readIndexFilter = transferIds.getReadIndexFilter();
        final AlignmentWriter writer = transferIds.getWriter();

        targetIds.clear();
        targetIds.putAll(transferIds.getTargetIds());

        // initialize too-many-hits output file
        final AlignmentTooManyHitsWriter tmhWriter =
                new AlignmentTooManyHitsWriter(outputFile, mParameter);

        try {

            // initialize numberOfReads
            if (numberOfReads < 0 && transferIds.numberOfReads != 0) {
                numberOfReads = transferIds.numberOfReads;
            }
            if (numberOfReads <= 0) {
                numberOfReads = numberOfReadsFromCommandLine;
            }
            if (numberOfReads <= 0) {
                System.err.println("Cannot determine number of reads. Must set property or provide reads file with -q");
                return;
            }


            // initialize quality filter
            qualityFilter = new PercentMismatchesQualityFilter();
            qualityFilter.setParameters(qualityFilterParameters);

            final int numAligns = scan(readIndexFilter, targetIds, writer, tmhWriter);

            System.out.println("Number of alignments written: " + numAligns);
            if (propagateQueryIds && !queryIds.isEmpty()) {
                // we collected query ids, let's write them to the header:
                writer.setQueryIdentifiers(queryIds);
            }
            if (propagateTargetIds && !targetIds.isEmpty()) {
                // we collected target ids, let's write them to the header:
                writer.setTargetIdentifiers(targetIds);
            }


            writer.setSmallestSplitQueryIndex(smallestQueryIndex);
            assert largestQueryIndex > -1 : "largestQueryIndex must be set (set with --number-of-reads when running from command line).";

            writer.setLargestSplitQueryIndex(largestQueryIndex);


        } finally {
            writer.close();
            tmhWriter.close();
        }
    }

    protected void evaluateStatistics(final AlignedSequence reference, final AlignedSequence query, final AlignmentStats stats) {
        final MutableString queryAligned = query.alignment;
        final MutableString targetAligned = reference.alignment;
        //  assert reference.alignedLength == query.alignedLength :"aligned length differ for queryIndex="+query.sequenceIdentifier;
        final int length = Math.max(query.alignedLength, reference.alignedLength);
        int numIndels = 0;
        int numMismatches = 0;
        for (int i = 0; i < length; i++) {
            final char queryBase = queryAligned.charAt(i);
            final char targetBase = targetAligned.charAt(i);
            if (queryBase == '-' && targetBase != '-' || queryBase != '-' && targetBase == '-') {
                numIndels++;
            }
            if (queryBase != '-' && targetBase != '-' && queryBase != targetBase) {
                numMismatches++;
            }
        }
        stats.numberOfIndels = numIndels;
        stats.numberOfMismatches = numMismatches;
    }


    public void setInputFile(final String inputFile) {
        this.inputFile = inputFile;
    }

    public void setOutputFile(final String outputFile) {
        this.outputFile = outputFile;
    }

    public String getOutputFile() {
        return outputFile;
    }

    public void setQueryReadIdsFilename(final String queryReadIdsFilename) {
        this.queryReadIdsFilename = queryReadIdsFilename;
    }

    public void setTargetReferenceIdsFilename(final String targetReferenceIdsFilename) {
        this.targetReferenceIdsFilename = targetReferenceIdsFilename;
    }

    public void setPropagateTargetIds(final boolean propagateTargetIds) {
        this.propagateTargetIds = propagateTargetIds;
    }

    public void setPropagateQueryIds(final boolean propagateQueryIds) {
        this.propagateQueryIds = propagateQueryIds;
    }

    public void setAmbiguityThreshold(final int mParameter) {
        this.mParameter = mParameter;
    }

    public void setQualityFilterParameters(final String qualityFilterParameters) {
        this.qualityFilterParameters = qualityFilterParameters;
    }

    public void setNumberOfReads(final int numberOfReads) {
        this.numberOfReads = numberOfReads;
    }

    public void setThirdPartyInput(final boolean thirdPartyInput) {
        this.thirdPartyInput = thirdPartyInput;
    }

    public void setSmallestQueryIndex(final int smallestQueryIndex) {
        this.smallestQueryIndex = smallestQueryIndex;
    }

    public void setLargestQueryIndex(final int largestQueryIndex) {
        this.largestQueryIndex = largestQueryIndex;
    }

    public void setNumberOfQuerySequences(final int numberOfReads) {
        this.numberOfReadsFromCommandLine = numberOfReads;
    }

    public static int getTargetIndex(final IndexedIdentifier targetIds,
                                     final CharSequence targetIdentifier,
                                     final boolean thirdPartyInput) {
        int targetIndex = -1;
        try {
            if (thirdPartyInput) {
                targetIndex = targetIds.registerIdentifier(new MutableString(targetIdentifier));
            } else {
                targetIndex = Integer.parseInt(targetIdentifier.toString());
            }
        } catch (NumberFormatException e) {
            if (targetIds != null) {
                final Integer object = targetIds.get(targetIdentifier);
                if (object == null) {
                    LOG.warn("Input file contains a target id that is not defined in the target compact reads: " + targetIdentifier);
                    targetIndex = targetIds.registerIdentifier(new MutableString(targetIdentifier));
                } else {
                    targetIndex = object;
                }
                if (targetIndex == -1) {
                    System.out.println("Cannot convert reference identifier to index. " + targetIdentifier);
                    System.exit(1);
                }
            }
        }
        return targetIndex;
    }

    /**
     * Compare read and reference sequences to determine sequence variations. The variations found
     * are appended to the alignment entry builder.
     *
     * @param currentEntry      alignment entry where variations will be stored.
     * @param alignmentLength   length of the sequence alignment (common length of reference and read sequences)
     * @param referenceSequence The reference sequence
     * @param readSequence      The read sequence
     * @param queryLength
     * @param baseQualities     ASCII encoded, remove 33 to get Phred quality score (see http://bioinformatics.oxfordjournals.org/cgi/reprint/btp352v1.pdf)
     */
    public static void extractSequenceVariations(final Alignments.AlignmentEntry.Builder currentEntry, final int alignmentLength,
                                                 final MutableString referenceSequence,
                                                 final MutableString readSequence,
                                                 final int readStartPosition,
                                                 final int queryLength, final boolean reverseStrand,
                                                 byte[] baseQualities) {
        //     System.out.printf("Extracting variations from %n%s%n%s%n",
        //             referenceSequence, readSequence);

        final MutableString from = new MutableString();
        final MutableString to = new MutableString();
        int variationPosition = Integer.MAX_VALUE;
        int minLength = Math.min(referenceSequence.length(), readSequence.length());
        minLength = Math.min(alignmentLength, minLength);
        // will record the number of gaps in the read, up to the variation position. We need to
        // subtract this number from the position to obtain the read index.
        int readIndexAdjustment = 0;
        int newAdjustment = 0;
        for (int position = 0; position < minLength; ++position) {

            final char referenceBase = referenceSequence.charAt(position);
            final char queryBase = readSequence.charAt(position);

            if (referenceBase != queryBase) {
                from.append(referenceBase);
                to.append(queryBase);
                variationPosition = Math.min(variationPosition, position);
                if (queryBase == '-') {

                    ++newAdjustment;
                }

            } else {
                appendNewSequenceVariation(currentEntry, from, to, variationPosition, readStartPosition, queryLength,
                        reverseStrand, readIndexAdjustment, baseQualities);
                variationPosition = Integer.MAX_VALUE;
                from.setLength(0);
                to.setLength(0);
                readIndexAdjustment = newAdjustment;
            }

        }
        appendNewSequenceVariation(currentEntry, from, to, variationPosition, readStartPosition, queryLength, reverseStrand, readIndexAdjustment, baseQualities);
    }

    /**
     * @param currentEntry
     * @param from
     * @param to
     * @param variationPosition
     * @param readStartPosition
     * @param queryLength
     * @param reverseStrand
     * @param readIndexAdjustment
     * @param baseQualities       ASCII encoded, remove 33 to get Phred quality score.
     */
    private static void appendNewSequenceVariation(final Alignments.AlignmentEntry.Builder currentEntry,
                                                   final MutableString from,
                                                   final MutableString to,
                                                   final int variationPosition,
                                                   final int readStartPosition,
                                                   final int queryLength,
                                                   final boolean reverseStrand,
                                                   final int readIndexAdjustment,
                                                   byte[] baseQualities) {
        if (variationPosition != Integer.MAX_VALUE) {
            final Alignments.SequenceVariation.Builder sequenceVariation =
                    Alignments.SequenceVariation.newBuilder();

            sequenceVariation.setFrom(from.toString());
            sequenceVariation.setTo(to.toString());
            sequenceVariation.setPosition(variationPosition + 1); // positions start at 1
            // calculate the readIndex, taking strand and query length into consideration:
            final int readIndex = (reverseStrand ? (queryLength - (variationPosition - readIndexAdjustment) - readStartPosition) :
                    variationPosition - readIndexAdjustment + readStartPosition);
            if (readIndex > queryLength) {
                assert readIndex <= queryLength : String.format(" readIndex %d must be smaller than read length %d .",
                        readIndex,
                        queryLength);
                LOG.warn(String.format(
                        "Ignoring sequence variations for a read since readIndex %d must be smaller than read length %d. query index=%d reference index=%d%n", readIndex,
                        queryLength, currentEntry.getQueryIndex(),
                        currentEntry.getTargetIndex()));
                //System.exit(1);
                return;
            }
            final int correctedReadIndex = readIndex + (reverseStrand ? 0 : 1);     // positions start at 1
            sequenceVariation.setReadIndex(correctedReadIndex);
            if (baseQualities != null) { // transfer read qualities for this sequence variation:
                byte[] toQualities = new byte[to.length()];
                int j = 0;

                for (int i = correctedReadIndex - 1; i < correctedReadIndex - 1 + to.length(); i++) {
                    if (i < baseQualities.length) {
                        toQualities[j++] = (byte) ((int) baseQualities[i] - 33);
                    } else {
                        LOG.warn(String.format("index i=%d too large max=%d", i, baseQualities.length));
                    }
                }


                sequenceVariation.setToQuality(ByteString.copyFrom(toQualities));
            }
            // do not offset if the match is in the reverse strand, since subtracting from the length takes care of offseting already.
            //        System.out.printf("Appending variation: %d %s/%s ", variationPosition, from, to);
            currentEntry.addSequenceVariations(sequenceVariation);
            // reset since they are used:
            from.setLength(0);
            to.setLength(0);
        }
    }

    public class TransferIds {
        private ReadSet readIndexFilter;
        private AlignmentWriter writer;
        private IndexedIdentifier targetIds;
        public int numberOfReads = -1;
        public int numberOfReadsForSplit;
        private int numberOfTargets;

        public ReadSet getReadIndexFilter() {
            return readIndexFilter;
        }

        public AlignmentWriter getWriter() {
            return writer;
        }

        public IndexedIdentifier getTargetIds() {
            return targetIds;
        }

        /**
         * Postcondition: output Ids.size() == maximumSequenceIndex + 1.
         */
        private ObjectArrayList<String> processIds(final String idsFilename) throws FileNotFoundException {
            final ObjectArrayList<String> ids = new ObjectArrayList<String>(500000);
            int minSequenceIndex = Integer.MAX_VALUE;
            int maxSequenceIndex = Integer.MIN_VALUE;
            ReadsReader idsReader = null;
            try {
                idsReader = new ReadsReader(new FileInputStream(idsFilename));

                boolean atLeastOneId = false;
                ids.size(500000);
                while (idsReader.hasNext()) {
                    final Reads.ReadEntry readEntry = idsReader.next();
                    final int readIndex = readEntry.getReadIndex();
                    if (readEntry.hasReadIdentifier()) {
                        // resize as necessary:
                        if (readIndex >= ids.size()) {
                            final double newSize = ids.size() * 1.5;
                            //  System.out.println("resizing to " + newSize);
                            ids.size((int) newSize);
                        }
                        // set element:
                        ids.set(readIndex, readEntry.getReadIdentifier());
                        atLeastOneId = true;
                    }
                    minSequenceIndex = Math.min(minSequenceIndex, readIndex);
                    maxSequenceIndex = Math.max(maxSequenceIndex, readIndex);
                }
            } finally {
                if (idsReader != null) {
                    try {
                        idsReader.close();
                    } catch (IOException e) {
                        LOG.warn("Error closing " + idsFilename, e);
                    }
                }
            }
            ids.size(maxSequenceIndex + 1);
            this.numberOfReadsForSplit = (maxSequenceIndex - minSequenceIndex) + 1;
            ids.trim();
            assert ids.size() == maxSequenceIndex + 1;
            return ids;
        }

        public AbstractAlignmentToCompactMode.TransferIds invoke() throws IOException {
            // setup multiplicity set:
            readIndexFilter = new ReadSet();
            if (readIndexFilterFile == null) {
                readIndexFilter = null;
            } else {
                readIndexFilter.load(readIndexFilterFile);
            }

            writer = new AlignmentWriter(outputFile);
            targetIds = new IndexedIdentifier();

            // first write reference ids to compact header, if these ids are provided on the command line:
            if (targetReferenceIdsFilename != null) {
                System.out.println("Scanning target file..");
                // read reference ids from file
                final ObjectArrayList<String> ids = processIds(targetReferenceIdsFilename);
                this.numberOfTargets = ids.size();
                System.out.println("Target file had " + this.numberOfTargets + " entries.");
                // write ids to header
                writer.setNumTargets(this.numberOfTargets);
                if (this.numberOfTargets > 0 && propagateTargetIds) {
                    writer.setTargetIdentifiersArray(ids.toArray(new String[ids.size()]));
                    System.out.println("Wrote " + ids.size() + " target ids to alignment header.");
                } else {
                    System.out.println("Target ids are not propagated to output header.");
                }

                for (int index = 0; index < ids.size(); index++) {
                    final String id = ids.get(index);
                    if (id != null) {
                        targetIds.put(new MutableString(id), index);
                    }
                }
            }
            // write query/read ids to compact header, if provided as well:
            if (queryReadIdsFilename != null) {
                // read read ids from file
                System.out.println("Scanning query file..");
                final ObjectArrayList<String> ids = processIds(queryReadIdsFilename);
                this.numberOfReads = ids.size();
                //   System.out.println("Query file had " + this.numberOfReads + " entries.");
                for (final String id : ids) {
                    if (id != null) {
                        queryIds.registerIdentifier(new MutableString(id));
                    }
                }
                // write ids to header
                writer.setNumQueries(this.numberOfReadsForSplit);
                if (this.numberOfReads > 0 && propagateQueryIds) {
                    writer.setQueryIdentifiersArray(ids.toArray(new String[ids.size()]));
                    System.out.println("Wrote " + ids.size() + " query ids to alignment header.");
                } else {
                    System.out.println("Query ids are not propagated to output header.");
                }
            }
            return this;
        }
    }
}
