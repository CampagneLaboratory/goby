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
import edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsWriter;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
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
            + "format to the compact alignment format.";

    /**
     * Native reads output from aligner.
     */
    protected String samBinaryFilename;

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
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        // configure baseclass
        super.configure(args);
        return this;
    }

    @Override
    protected int scan(final ReadSet readIndexFilter, final IndexedIdentifier targetIds, final AlignmentWriter writer,
                       final AlignmentTooManyHitsWriter tmhWriter) throws IOException {
        int numAligns = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        final SAMFileReader parser = new SAMFileReader(new File(inputFile));
        parser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);


        final int[] readLengths = new int[numberOfReads];
        final CloseableIterator<SAMRecord> recordCloseableIterator = parser.iterator();

        progress.start();

        while (recordCloseableIterator.hasNext()) {
            final SAMRecord samRecord = recordCloseableIterator.next();
            final int queryIndex = Integer.parseInt(samRecord.getReadName());

            final int depth = samRecord.getReadLength();
            // save length
            readLengths[queryIndex] = depth;

            // if SAM reports read is unmapped (we don't know how or why), skip record
            if (samRecord.getReadUnmappedFlag()) {
                continue;
            }

            // SAM mismatch info must be available, o.w. skip record
            final String mismatches = (String) samRecord.getAttribute("MD");
            // System.out.println("mismatches: " + mismatches);
            if (mismatches == null) {
                continue;
            }

            final int targetIndex = Integer.parseInt(samRecord.getReferenceName());
            // positions reported by BWA appear to start at 1. We convert to start at zero.
            final int position = samRecord.getAlignmentStart() - 1;
            final boolean reverseStrand = samRecord.getReadNegativeStrandFlag();
            final List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();

            float score = 0;
            int targetAlignedLength = 0;
            int numMismatches = 0;
            int numIndels = 0;
            int queryAlignedLength = 0;
            int multiplicity = 1;

            // count the number of mismatches:
            for (final char c : mismatches.toCharArray()) {
                if (!Character.isDigit(c)) {
                    numMismatches++;
                }
            }

            for (final CigarElement cigar : samRecord.getCigar().getCigarElements()) {
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
            for (final AlignmentBlock block : blocks) {
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

            String attributeXA = (String) samRecord.getAttribute("XA");
            String sequence = samRecord.getReadString();
            extractSequenceVariations(attributeXA, sequence, currentEntry);
            final Alignments.AlignmentEntry alignmentEntry = currentEntry.build();

            final int numTotalHits = (Integer) samRecord.getAttribute("X0");

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

    private void extractSequenceVariations(String attributeXA, String sequence, Alignments.AlignmentEntry.Builder currentEntry) {
      //  System.out.println("Found: " + attributeXA + " " + sequence);
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws JSAPException error parsing
     * @throws IOException   error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new SAMToCompactMode().configure(args).execute();
    }
}
