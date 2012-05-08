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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.AlignmentWriterImpl;
import edu.cornell.med.icb.goby.alignments.perms.QueryIndexPermutation;
import edu.cornell.med.icb.goby.compression.MessageChunksWriter;
import edu.cornell.med.icb.goby.readers.sam.SamComparison;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionRegistry;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import it.unimi.dsi.Util;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * Converts alignments in the SAM or BAM format to the compact alignment format.
 *
 * @author Kevin Dorff
 * @author Fabien Campagne
 */
public class SAMComparisonMode extends AbstractGobyMode {

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(SAMComparisonMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "sam-comparison";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Goby's modes sam-to-compact and compact-to-sam are " +
            "capable of taking a source SAM/BAM file, converting it to Goby Compact Alignment, " +
            "then converting that Goby Compact Alignment to SAM/BAM. This mode is written to " +
            "verify the the source SAM/BAM against the resultant SAM/BAM (optionally taking the intermediate " +
            "Goby Compact Alignment as an argument). It is expected that the source and destination SAM/BAM files " +
            "will contain the same data with the exception of lines in the source files that are 'unmapped', which " +
            "will not appear in the destination SAM/BAM file. The fidelity of the two files depends on the " +
            "parameters given to sam-to-compact (keep soft clips, keep quality, etc.). " +
            "This is NOT a general purpose SAM/BAM comparison tool.";

    private File sourceBamFile;
    private File destBamFile;
    private String destGobyBasename;
    private boolean mappedQualitiesPreserved;
    private boolean softClipsPreserved;
    private boolean checkMate;
    private boolean canonicalMdzForComparison = true;

    private boolean runningFromCommandLine = false;
    private boolean debug = false;

    // --------- FOR COMPARISON ---------
    private SamComparison samComparison = new SamComparison();

    @RegisterThis
    public static DynamicOptionClient doc = new DynamicOptionClient(SAMComparisonMode.class,
            "ignore-read-origin:boolean, When this flag is true do not import read groups.:false"
    );

    public static DynamicOptionClient doc() {
        return doc;
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
     * Get the expected/source sam alignment file. This contains the values we EXPECT TO FIND, ie, the
     * source SAM/BAM file.
     * @return the expected/source sam alignment file.
     */
    public File getSourceBamFile() {
        return sourceBamFile;
    }

    /**
     * Set the expected/source sam alignment file. This contains the values we EXPECT TO FIND, ie, the
     * source SAM/BAM file.
     * @param sourceBamFile the expected/source sam alignment file.
     */
    public void setSourceBamFile(final File sourceBamFile) {
        this.sourceBamFile = sourceBamFile;
    }

    /**
     * Get the ACTUAL sam alignment file. This contains the alignments containing the values we will compare against
     * the expected (source).
     * @return the actual/dest sam dest file
     */
    public File getDestBamFile() {
        return destBamFile;
    }

    /**
     * Set the ACTUAL sam alignment file. This contains the alignments containing the values we will compare against
     * the expected (source).
     * @param destBamFile the actual/dest sam dest file
     */
    public void setDestBamFile(final File destBamFile) {
        this.destBamFile = destBamFile;
    }

    /**
     * Get the Goby Compact Alignment basename (optional). This is the source of intermediate alignments that were
     * created from the expected alignments and was used to create the actual alignments. Providing this will generate
     * some additional output when there are comparison failures.
     * @return the goby alignment basename
     */
    public String getDestGobyBasename() {
        return destGobyBasename;
    }

    /**
     * Set the Goby Compact Alignment basename (optional). This is the source of intermediate alignments that were
     * created from the expected alignments and was used to create the actual alignments. Providing this will generate
     * some additional output when there are comparison failures.
     * @param destGobyBasename the goby alignment basename
     */
    public void setDestGobyBasename(final String destGobyBasename) {
        this.destGobyBasename = destGobyBasename;
    }

    /**
     * Get if it is assumed that the compact file created from the BAM/SAM
     * file preserved mapped qualities.
     * @return if it is assumed ...
     */
    public boolean isMappedQualitiesPreserved() {
        return mappedQualitiesPreserved;
    }

    /**
     * Set if it is assumed that the compact file created from the BAM/SAM
     * file preserved mapped qualities.
     * @param mappedQualitiesPreserved if it is assumed...
     */
    public void setMappedQualitiesPreserved(final boolean mappedQualitiesPreserved) {
        this.mappedQualitiesPreserved = mappedQualitiesPreserved;
    }

    /**
     * Get if it is assumed that the compact file created from the BAM/SAM
     * file preserved soft clips.
     * @return if it is assumed ...
     */
    public boolean isSoftClipsPreserved() {
        return softClipsPreserved;
    }

    /**
     * Set if it is assumed that the compact file created from the BAM/SAM
     * file preserved soft clips.
     * @param softClipsPreserved if it is assumed ...
     */
    public void setSoftClipsPreserved(final boolean softClipsPreserved) {
        this.softClipsPreserved = softClipsPreserved;
    }

    /**
     * Get if the details about mate reads will be checked.
     * If the source SAM/BAM file is a complete file you can set this to true,
     * if you are using an incomplete source SAM/BAM file, this should be
     * set to false. Default is false.
     * @return if mates will be checked
     */
    public boolean isCheckMate() {
        return checkMate;
    }

    /**
     * Set if the details about mate reads will be checked.
     * If the source SAM/BAM file is a complete file you can set this to true,
     * if you are using an incomplete source SAM/BAM file, this should be
     * set to false. Default is false.
     * @return if mates will be checked
     */
    public void setCheckMate(final boolean checkMate) {
        this.checkMate = checkMate;
    }

    /**
     * Get if canonical MD:Z comparisons will be made.
     * When true, the source and destination MD:Z values will be passed through an algorithm
     * to make them canonical (place 0's in places where 0's should exist but might not).
     * By default this is enabled.
     * @return if ...
     */
    public boolean isCanonicalMdzForComparison() {
        return canonicalMdzForComparison;
    }

    /**
     * Set if canonical MD:Z comparisons will be made.
     * When true, the source and destination MD:Z values will be passed through an algorithm
     * to make them canonical (place 0's in places where 0's should exist but might not).
     * By default this is enabled.
     * @param canonicalMdzForComparison if ...
     */
    public void setCanonicalMdzForComparison(final boolean canonicalMdzForComparison) {
        this.canonicalMdzForComparison = canonicalMdzForComparison;
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

        final JSAPResult jsapResult = parseJsapArguments(args);
        sourceBamFile = jsapResult.getFile("source-bam");
        destBamFile = jsapResult.getFile("destination-bam");
        destGobyBasename = jsapResult.getString("destination-goby");
        mappedQualitiesPreserved = jsapResult.getBoolean("mapped-qualities-preserved");
        softClipsPreserved = jsapResult.getBoolean("soft-clips-preserved");
        checkMate = jsapResult.getBoolean("check-mate");
        canonicalMdzForComparison = jsapResult.getBoolean("canonical-mdz");

        // is way too slow to run unintentionally in production!
        DynamicOptionRegistry.register(MessageChunksWriter.doc());
        DynamicOptionRegistry.register(AlignmentWriterImpl.doc());
        DynamicOptionRegistry.register(QueryIndexPermutation.doc());
        runningFromCommandLine = true;
        return this;
    }

    @Override
    public void execute() throws IOException {
        // don't even dare go through the debugging code if log4j was not configured. The debug code
        debug = Util.log4JIsConfigured();

        if (!sourceBamFile.exists()) {
            System.err.println("--source-bam SAM/BAM file couldn't be found " + sourceBamFile.toString());
            System.exit(-1);
        }
        if (!destBamFile.exists()) {
            System.err.println("--destination-bam SAM/BAM file couldn't be found " + destBamFile.toString());
            System.exit(-1);
        }

        System.out.println("Comparing source bam and destination bam");
        final SAMFileReader sourceParser = new SAMFileReader(new FileInputStream(sourceBamFile));
        final SAMFileReader destParser = new SAMFileReader(new FileInputStream(destBamFile));
        // We need to set the validation to silent because an incomplete file (if the source isn't the entire file)
        // we can see errors that wouldn't exist in a real conversion.
        sourceParser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        destParser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SAMRecordIterator sourceIterator = sourceParser.iterator();
        final SAMRecordIterator destIterator = destParser.iterator();
        AlignmentReaderImpl gobyReader = null;
        if (destGobyBasename != null && destGobyBasename.length() > 0) {
            try {
                gobyReader = new AlignmentReaderImpl(destGobyBasename);
            } catch (IOException e) {
                LOG.error("Error opening --destination-goby file " + destGobyBasename, e);
                return;
            }
        }

        if (gobyReader == null && !mappedQualitiesPreserved) {
            System.err.println("WARNING: If no --destination-goby file is provided and " +
                    "--mapped-qualities-preserved is not selected, nothing about qualities will be compared.");
        }

        final ProgressLogger progress = new ProgressLogger(LOG);
        progress.displayFreeMemory = true;
        progress.start();
        samComparison.setMappedQualitiesPreserved(mappedQualitiesPreserved);
        samComparison.setSoftClipsPreserved(softClipsPreserved);
        samComparison.setCheckMate(checkMate);
        samComparison.setCanonicalMdzForComparison(canonicalMdzForComparison);
        samComparison.reset();
        int numRecordsSkipped = 0;
        while (sourceIterator.hasNext()) {
            samComparison.setExpectedSamRecord(sourceIterator.next());
            if (samComparison.getExpectedSamRecord().getReadUnmappedFlag()) {
                // We don't store unmapped reads, so skip this source record
                numRecordsSkipped++;
                continue;
            }
            if (!destIterator.hasNext()) {
                LOG.error("Not enough records in --destination-bam SAM/BAM file");
                return;
            }
            samComparison.setActualSamRecord(destIterator.next());
            if (gobyReader != null) {
                if (!gobyReader.hasNext()) {
                    LOG.error("Not enough records in goby compact-alignment file");
                    return;
                }
                samComparison.setGobyAlignment(gobyReader.next());
            }

            samComparison.compareSamRecords();
            progress.lightUpdate();
        }
        progress.stop();
        System.out.println("Number of records compared   : " + samComparison.getReadNum());
        System.out.println("Number of records unmapped   : " + numRecordsSkipped);
        System.out.println("Number of comparison failures: " + samComparison.getComparisonFailureCount());
    }
}
