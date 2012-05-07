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
        samComparison.mappedQualitiesPreserved = mappedQualitiesPreserved;
        samComparison.softClipsPreserved = softClipsPreserved;
        samComparison.checkMate = true;
        samComparison.reset();
        samComparison.readNum = 0;
        while (sourceIterator.hasNext()) {
            samComparison.expectedSamRecord = sourceIterator.next();
            if (samComparison.expectedSamRecord.getReadUnmappedFlag()) {
                // We don't store unmapped reads, so skip this source record
                continue;
            }
            if (!destIterator.hasNext()) {
                LOG.error("Not enough records in --destination-bam SAM/BAM file");
                return;
            }
            samComparison.actualSamRecord = destIterator.next();
            if (gobyReader != null) {
                if (!gobyReader.hasNext()) {
                    LOG.error("Not enough records in goby compact-alignment file");
                    return;
                }
                samComparison.gobyAlignment = gobyReader.next();
            }

            samComparison.compareSamRecords();
            progress.lightUpdate();
        }
        progress.stop();
        System.out.println("Number of records compared   : " + samComparison.readNum);
        System.out.println("Number of comparison failures: " + samComparison.comparisonFailureCount);
    }
}
