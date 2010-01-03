/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import edu.cornell.med.icb.alignments.AlignmentTooManyHitsWriter;
import edu.cornell.med.icb.alignments.AlignmentWriter;
import edu.cornell.med.icb.alignments.Alignments;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.reads.ReadSet;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;

/**
 * Converts binary BWA alignments in the SAM format to the compact alignment format.
 *
 * @author Fabien Campagne
 */
public class SAMToCompactMode extends AbstractAlignmentToCompactMode { // AbstractGobyMode { //LastToCompactMode {
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "sam-to-compact";
    public static final String MODE_DESCRIPTION = "Converts binary BWA alignments in the SAM format to the compact alignment format.";

    public String getSamBinaryFilename() {
        return samBinaryFilename;
    }

    //
    public void setSamBinaryFilename(String samBinaryFilename) {
        this.samBinaryFilename = samBinaryFilename;
    }


    // Native reads output from aligner
    protected String samBinaryFilename;

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SAMToCompactMode.class);


    public String getModeName() {
        return MODE_NAME;
    }

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
        // configure baseclass
        super.configure(args);
        //
        return this;
    }

    protected int scan(ReadSet readIndexFilter, IndexedIdentifier targetIds, AlignmentWriter writer,
                     AlignmentTooManyHitsWriter tmhWriter) throws IOException {
        int numAligns = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        final SAMFileReader parser = new SAMFileReader(new File(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);


        int[] readLengths = new int[numberOfReads];
        final CloseableIterator<SAMRecord> recordCloseableIterator = parser.iterator();

        progress.start();

        while (recordCloseableIterator.hasNext()) {
            SAMRecord samRecord = recordCloseableIterator.next();
            int queryIndex = Integer.parseInt(samRecord.getReadName());

            final int depth = samRecord.getReadLength();
            // save length
            readLengths[queryIndex] = depth;

            // if SAM reports read is unmapped (we don't know how or why), skip record
            if (samRecord.getReadUnmappedFlag()) {
                continue;
            }

            // SAM mismatch info must be available, o.w. skip record
            String mismatches = (String) samRecord.getAttribute("MD");
            // System.out.println("mismatches: " + mismatches);
            if (mismatches == null) {
                continue;
            }

            int targetIndex = Integer.parseInt(samRecord.getReferenceName());
            // positions reported by BWA appear to start at 1. We convert to start at zero.
            int position = samRecord.getAlignmentStart()-1;
            boolean reverseStrand = samRecord.getReadNegativeStrandFlag();
            java.util.List<net.sf.samtools.AlignmentBlock> blocks = samRecord.getAlignmentBlocks();

            float score = 0;
            int targetAlignedLength = 0;
            int numMismatches = 0;
            int numIndels = 0;
            int queryAlignedLength = 0;
            int multiplicity = 1;

            // count the number of mismatches:
            for (char c : mismatches.toCharArray()) {
                if (!Character.isDigit(c)) {
                    numMismatches++;
                }
            }
            for (CigarElement cigar : samRecord.getCigar().getCigarElements()) {
                final int length = cigar.getLength();
                switch (cigar.getOperator()) {
                    case M:
                        // match or mismatch?!   CIGAR cannot differentiate.
                        // This means here that the score rewards mutations as much as matches.  We would
                        // have to parse the sequence to determine what is what.
                        score += length;
                        break;
                    case I:
                        // insertion to the reference :
                        score -= length;
                        numIndels += length;
                        break;
                    case P:
                        //padding, no score impact:
                        break;
                    case D:
                        // deletion from the reference
                        score -= length;
                        numIndels += length;
                        break;
                }
            }
            score -= numMismatches;
            for (AlignmentBlock block : blocks) {
                // System.out.println("block length: " + block.getLength() + " readStart: " + block.getReadStart() + " referenceStart: " + block.getReferenceStart());
                targetAlignedLength += block.getLength();
                queryAlignedLength += block.getLength();
            }

            // we have a multiplicity filter. Use it to determine multiplicity.
            if (readIndexFilter != null) {
                /* Multiplicity of a read is the number of times the (exact) sequence
                  of the read is identically repeated across a sample file.  The filter
                  removes duplicates to avoid repeating the same alignments.  Once
                  aligned, these are recorded multiplicity times. */
                multiplicity = readIndexFilter.getMultiplicity(queryIndex);
            }

            // the record represents a mapped read..
            final Alignments.AlignmentEntry.Builder currentEntry = Alignments.AlignmentEntry.newBuilder();

            //      System.out.println("score: " + score);
            currentEntry.setNumberOfIndels(numIndels);
            currentEntry.setNumberOfMismatches(numMismatches);
            currentEntry.setMatchingReverseStrand(reverseStrand);
            currentEntry.setMultiplicity(multiplicity);
            currentEntry.setPosition(position);
            currentEntry.setQueryAlignedLength(queryAlignedLength);
            currentEntry.setQueryIndex(queryIndex);
            currentEntry.setScore(score);
            currentEntry.setTargetAlignedLength(targetAlignedLength);
            currentEntry.setTargetIndex(targetIndex);
            final Alignments.AlignmentEntry alignmentEntry = currentEntry.build();

            int numTotalHits = (Integer) samRecord.getAttribute("X0");

            if (qualityFilter.keepEntry(depth, alignmentEntry)) {
                // only write the entry if it is not ambiguous. i.e. less than or equal to mParameter
                if (numTotalHits <= mParameter) {
                    writer.appendEntry(alignmentEntry);
                    numAligns += multiplicity;
                }
                // TMH writer adds the alignment entry only if hits > thresh
                tmhWriter.append(queryIndex, numTotalHits, depth);
            }

            progress.lightUpdate();

        }
        recordCloseableIterator.close();

        // TODO write statistics to writer.
        // stu 090817 - mimicking LastToCompactMode.scan() statistics
        if (readIndexFilter != null) {
            writer.putStatistic("keep-filter-filename", readIndexFilterFile.getName());
        }
        writer.putStatistic("number-of-entries-written", numAligns);
        writer.printStats(System.out);
        writer.setQueryLengths(readLengths);
        progress.stop();
        return numAligns;
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
        new SAMToCompactMode().configure(args).execute();
    }


}
