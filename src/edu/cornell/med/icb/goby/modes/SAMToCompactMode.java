/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsWriter;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.perms.ReadNameToIndex;
import edu.cornell.med.icb.goby.reads.QualityEncoding;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.Util;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Converts binary BWA alignments in the SAM format to the compact alignment format.
 *
 * @author Fabien Campagne
 */
public class SAMToCompactMode extends AbstractAlignmentToCompactMode {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SAMToCompactMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sam-to-compact";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts binary BWA alignments in the SAM "
            + "format to the compact alignment format (new version).";

    /**
     * Native reads output from aligner.
     */
    protected String samBinaryFilename;

    private int dummyQueryIndex;

    private boolean bsmap;

    private QualityEncoding qualityEncoding = QualityEncoding.SANGER;
    /**
     * Flag to indicate if log4j was configured.
     */
    private boolean debug;

    public String getSamBinaryFilename() {
        return samBinaryFilename;
    }

    public void setSamBinaryFilename(final String samBinaryFilename) {
        this.samBinaryFilename = samBinaryFilename;
    }


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Get the quality encoding scale used for the input fastq file.
     *
     * @return the quality encoding scale used for the input fastq file
     */
    public QualityEncoding getQualityEncoding() {
        return qualityEncoding;
    }

    /**
     * Set the quality encoding scale to be used for the input fastq file.
     * Acceptable values are "Illumina", "Sanger", and "Solexa".
     *
     * @param qualityEncoding the quality encoding scale to be used for the input fastq file
     */
    public void setQualityEncoding(final QualityEncoding qualityEncoding) {
        this.qualityEncoding = qualityEncoding;
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
        // configure baseclass
        super.configure(args);
        final JSAPResult jsapResult = parseJsapArguments(args);
        bsmap = jsapResult.getBoolean("bsmap");

        numberOfReadsFromCommandLine = jsapResult.getInt("number-of-reads");
        qualityEncoding = QualityEncoding.valueOf(jsapResult.getString("quality-encoding").toUpperCase());
        sortedInput = jsapResult.getBoolean("sorted");
        this.largestQueryIndex = numberOfReadsFromCommandLine;
        this.smallestQueryIndex = 0;
        // don't even dare go through the debugging code if log4j was not configured. The debug code
        // is way too slow to run unintentionally in production!
        debug = Util.log4JIsConfigured();
        return this;
    }

    boolean sortedInput;

    @Override
    protected int scan(final ReadSet readIndexFilter, final IndexedIdentifier targetIds,
                       final AlignmentWriter writer, final AlignmentTooManyHitsWriter tmhWriter)
            throws IOException {
        int numAligns = 0;

        final ProgressLogger progress = new ProgressLogger(LOG);
        // the following is required to set validation to SILENT before loading the header (done in the SAMFileReader constructor)
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        final SAMFileReader parser = new SAMFileReader(new File(inputFile));


        progress.start();

        final SAMRecordIterator samIterator = parser.iterator();

        SamHelper samHelper = new SamHelper();
        samHelper.setQualityEncoding(qualityEncoding);
        numberOfReads = 0;

        // int stopEarly = 0;
        SAMRecord prevRecord = null;

        final SAMFileHeader fileHeader = parser.getFileHeader();
        if (sortedInput) {
            // if the input is sorted, request creation of the index when writing the alignment.
            final int numTargets = fileHeader.getSequenceDictionary().size();
            final int[] targetLengths = new int[numTargets];
            for (int i = 0; i < numTargets; i++) {
                final SAMSequenceRecord seq = fileHeader.getSequence(i);
                final int targetIndex = getTargetIndex(targetIds, seq.getSequenceName(), thirdPartyInput);
                targetLengths[targetIndex] = seq.getSequenceLength();
            }
            writer.setTargetLengths(targetLengths);
            writer.setSorted(true);
        }
        while (samIterator.hasNext()) {
            samHelper.reset();
            numberOfReads++;
            final SAMRecord samRecord = samIterator.next();
            if (samRecord.getReadUnmappedFlag()) {
                if (debug && LOG.isDebugEnabled()) {
                    LOG.debug(String.format("NOT keeping unmapped read %s", samRecord.getReadName()));
                }
                continue;
            }

            if (sortedInput) {
                // check that input entries are indeed in sort order. Abort otherwise.
                if (prevRecord != null) {
                    final int compare = prevRecord.getAlignmentStart() - samRecord.getAlignmentStart();//  samComparator.compare(prevRecord, samRecord);
                    if (compare > 0) {
                        final String message = String.format("record %s has position before previous record: %s",
                                samRecord.toString(), prevRecord.toString()
                        );
                        System.err.println("You cannot specify --sorted when the input file is not sorted. For instance: " + message);

                        LOG.warn(message);
                        System.exit(0);
                    }
                }
            }

            prevRecord = samRecord;

            // try to determine readMaxOccurence, the maximum number of times a read name occurs in a complete alignment.
            int readMaxOccurence = 1;
            final boolean readIsPaired = samRecord.getReadPairedFlag();
            final boolean b = readIsPaired && !samRecord.getReadUnmappedFlag();
            if (b) {
                // if the reads are paired, we expect to see the read name  at least twice.
                readMaxOccurence++;
                // Unfortunately the SAM/BAM format does not provide the exact number of times
                // a read matched the reference sequence. We could find this number in an non sorted BAM file
                // by counting how many times the same read name appears in a continous block of constant read name.
                // However, for sorted input, we need this number to know when to stop
                // keep a given read name in memory with its associated query index.
                // Since we can't keep all the read names in memory continuously (these are strings and
                // consume much memory), it is unclear how to determine  query-index-occurrences in a sorted
                // SAM/BAM file that contains multiple best hits for read or mate. We can handle these cases
                // correctly when working directly in the aligner and writing Goby format because the information
                // is available at the time of alignment, but discarded afterwards.
                //TODO see if we find a solution to better handle ambiguous matches.
            }

            final Object xoString = samRecord.getAttribute("X0");
            final int numTotalHits = xoString == null ? 1 : (Integer) xoString;
            readMaxOccurence = Math.max(numTotalHits * (readIsPaired ? 2 : 1), readMaxOccurence);
            final String readName = samRecord.getReadName();

            final int queryIndex = nameToQueryIndices.getQueryIndex(readName, readMaxOccurence);

            final int targetIndex = getTargetIndex(targetIds, samRecord.getReferenceName(), thirdPartyInput);

            final int fragmentIndex;
            final int mateFragmentIndex;
            if (readIsPaired) {
                if (samRecord.getFirstOfPairFlag()) {
                    fragmentIndex = 0;
                    mateFragmentIndex = 1;
                } else {
                    fragmentIndex = 1;
                    mateFragmentIndex = 0;
                }
            } else {
                fragmentIndex = 0;
                mateFragmentIndex = 1;
            }

            if (bsmap) {
                // reference is provided in attribute XR
                final String specifiedReference = (String) samRecord.getAttribute("XR");
                final String directions = (String) samRecord.getAttribute("XS");
                final boolean reverseStrand = directions.equals("-+") || directions.equals("+-");
                samHelper.setSourceWithReference(queryIndex, samRecord.getReadString(), samRecord.getBaseQualityString(), specifiedReference, samRecord.getAlignmentStart(), reverseStrand);
            } else {
                samHelper.setSource(queryIndex, samRecord.getReadString(), samRecord.getBaseQualityString(), samRecord.getCigarString(), (String) samRecord.getAttribute("MD"), samRecord.getAlignmentStart(), samRecord.getReadNegativeStrandFlag());
            }

            // positions reported by BWA appear to start at 1. We convert to start at zero.
            int multiplicity = 1;

            // we have a multiplicity filter. Use it to determine multiplicity.
            if (readIndexFilter != null) {
                /* Multiplicity of a read is the number of times the (exact) sequence
          of the read is identically repeated across a sample file.  The filter
          removes duplicates to avoid repeating the same alignments.  Once
          aligned, these are recorded multiplicity times. */
                multiplicity = readIndexFilter.getMultiplicity(queryIndex);
            }
            largestQueryIndex = Math.max(queryIndex, largestQueryIndex);
            smallestQueryIndex = Math.min(queryIndex, smallestQueryIndex);

            // the record represents a mapped read..
            final Alignments.AlignmentEntry.Builder currentEntry = Alignments.AlignmentEntry.newBuilder();

            currentEntry.setMultiplicity(multiplicity);
            currentEntry.setQueryIndex(samHelper.getQueryIndex());
            currentEntry.setTargetIndex(targetIndex);
            currentEntry.setPosition(samHelper.getPosition());
            currentEntry.setQueryPosition(samHelper.getQueryPosition());
            currentEntry.setFragmentIndex(fragmentIndex);
            currentEntry.setQueryLength(samHelper.getQueryLength());
            currentEntry.setScore(samHelper.getScore());
            currentEntry.setNumberOfIndels(samHelper.getNumDeletions() + samHelper.getNumInsertions());
            currentEntry.setNumberOfMismatches(samHelper.getNumMisMatches());
            currentEntry.setMatchingReverseStrand(samHelper.isReverseStrand());
            currentEntry.setQueryAlignedLength(samHelper.getQueryAlignedLength());
            currentEntry.setTargetAlignedLength(samHelper.getTargetAlignedLength());
            currentEntry.setPairFlags(samRecord.getFlags());

            if (readIsPaired) {
                if (!samRecord.getMateUnmappedFlag()) {
                    final Alignments.RelatedAlignmentEntry.Builder relatedBuilder =
                            Alignments.RelatedAlignmentEntry.newBuilder();
                    final int mateTargetIndex = getTargetIndex(targetIds, samRecord.getMateReferenceName(), thirdPartyInput);
                    final int mateAlignmentStart = samRecord.getMateAlignmentStart() - 1; // Goby is 0-based
                    relatedBuilder.setFragmentIndex(mateFragmentIndex);
                    relatedBuilder.setPosition(mateAlignmentStart);
                    relatedBuilder.setTargetIndex(mateTargetIndex);
                    currentEntry.setPairAlignmentLink(relatedBuilder);
                }
            }

            for (final SamSequenceVariation var : samHelper.getSequenceVariations()) {
                appendNewSequenceVariation(currentEntry, var, samHelper.getQueryLength());
                if (debug && LOG.isDebugEnabled()) {
                    LOG.debug(String.format("Added seqvar=%s for queryIndex=%d to alignment", var.toString(), queryIndex));
                }
            }

            final Alignments.AlignmentEntry alignmentEntry = currentEntry.build();

            if (qualityFilter.keepEntry(samHelper.getQueryLength(), alignmentEntry)) {
                // only write the entry if it is not ambiguous. i.e. less than or equal to mParameter
                if (numTotalHits <= mParameter) {

                    writer.appendEntry(alignmentEntry);
                    numAligns += multiplicity;
                    if (debug && LOG.isDebugEnabled()) {
                        LOG.debug(String.format("Added queryIdndex=%d to alignment", queryIndex));
                    }
                } else {
                    // TMH writer adds the alignment entry only if hits > thresh
                    tmhWriter.append(queryIndex, numTotalHits, samHelper.getQueryLength());
                    if (debug && LOG.isDebugEnabled()) {
                        LOG.debug(String.format("Added queryIndex=%d to TMH", queryIndex));
                    }
                    // remove the query name from memory since we are not writing these entries anyway
                    while (queryIndex == nameToQueryIndices.getQueryIndex(samRecord.getReadName(), 0)) ;
                }
            } else {
                if (debug && LOG.isDebugEnabled()) {
                    LOG.debug(String.format("NOT keeping queryIndex=%d", queryIndex));
                }
            }

            progress.lightUpdate();

        }
        samIterator.close();

        // stu 090817 - mimicking LastToCompactMode.scan() statistics
        if (readIndexFilter != null) {
            writer.putStatistic("keep-filter-filename", readIndexFilterFile.getName());
        }
        writer.putStatistic("number-of-entries-written", numAligns);
        writer.setNumQueries(Math.max(numberOfReads, numberOfReadsFromCommandLine));
        writer.printStats(System.out);

        // write information from SAM file header
        final SAMFileHeader samHeader = fileHeader;
        final SAMSequenceDictionary samSequenceDictionary = samHeader.getSequenceDictionary();
        final List<SAMSequenceRecord> samSequenceRecords = samSequenceDictionary.getSequences();
        int targetCount = targetIds.size();
        if (targetIds.size() != 0 && (targetIds.size() != samSequenceRecords.size())) {

            LOG.warn("targets: " + targetIds.size() + ", records: " + samSequenceRecords.size());
        }
        targetCount = Math.max(samSequenceRecords.size(), targetCount);
        final int[] targetLengths = new int[targetCount];
        for (final SAMSequenceRecord samSequenceRecord : samSequenceRecords) {
            final int index = samSequenceRecord.getSequenceIndex();
            if (debug && LOG.isDebugEnabled()) {
                LOG.debug("Sam record: " + samSequenceRecord.getSequenceName() + " at " + index);
            }
            targetLengths[index] = samSequenceRecord.getSequenceLength();
        }

        writer.setTargetLengths(targetLengths);
        progress.stop();
        return numAligns;
    }

    static void appendNewSequenceVariation(
            final Alignments.AlignmentEntry.Builder currentEntry,
            final SamSequenceVariation var, final int queryLength) {

        int readIndex = var.getReadIndex();
        if (readIndex > queryLength) {
            assert readIndex <= queryLength : String.format(" readIndex %d must be smaller than read length %d .",
                    readIndex, queryLength);
            LOG.warn(String.format(
                    "Ignoring sequence variations for a read since readIndex %d must be smaller than read length %d. query index=%d reference index=%d%n",
                    readIndex, queryLength, currentEntry.getQueryIndex(), currentEntry.getTargetIndex()));
            return;
        }

        final Alignments.SequenceVariation.Builder sequenceVariation =
                Alignments.SequenceVariation.newBuilder();

        sequenceVariation.setFrom(var.getFromString().toString());
        sequenceVariation.setTo(var.getToString().toString());
        sequenceVariation.setPosition(var.getRefPosition()); // positions start at 1
        sequenceVariation.setReadIndex(readIndex);  // readIndex starts at 1
        if (var.getQual() != null) {
            sequenceVariation.setToQuality(ByteString.copyFrom(var.getQualByteArray()));
        }
        currentEntry.addSequenceVariations(sequenceVariation);
    }

    ReadNameToIndex nameToQueryIndices = new ReadNameToIndex("ignore-this-for-now");

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new SAMToCompactMode().configure(args).execute();
    }
}