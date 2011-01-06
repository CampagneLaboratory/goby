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

package edu.cornell.med.icb.goby.aligners;

import edu.cornell.med.icb.goby.config.GobyConfiguration;
import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.goby.modes.SAMToCompactMode;
import edu.cornell.med.icb.goby.util.LoggingOutputStream;
import edu.mssm.crover.cli.CLI;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.RandomStringUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Level;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Facade to the Burrows-Wheeler Aligner (BWA) aligner. For details about BWA,
 * see <a href="http://bio-bwa.sourceforge.net/">http://bio-bwa.sourceforge.net/</a>
 * <p/>
 *
 * @author Fabien Campagne
 *         Date: Jul 9, 2009
 *         Time: 5:21:08 PM
 */
public class BWAAligner extends AbstractAligner {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(BWAAligner.class);

    private static final String BWA_EXEC = SystemUtils.IS_OS_WINDOWS ? "bwa.exe" : "bwa";
    private static final String BWA_INDEX_ARG = "index";
    private static final String BWA_ALIGN_ARG = "aln";
    private static final String BWA_SAMSE_ARG = "samse";
    private static final int DEFAULT_SEED_LENGTH = 35;
    private static final boolean NATIVE_GOBY_SUPPORTED = false;

    private String databasePrefix;
    private String saiBinaryFilename = "";
    private String samBinaryFilename = "";

    private int seedLength = DEFAULT_SEED_LENGTH;

    public BWAAligner() {
        super();
        extensions = new String[]{"rsa", "rpac", "rbwt", "pac", "bwt", "ann", "amb", "sa"};
    }

    public void setConfiguration(final Configuration configuration) {
        pathToExecutables = configuration.getString(GobyConfiguration.EXECUTABLE_PATH_BWA, "");
    }

    public String getDefaultDbNameForReferenceFile(final File referenceFile) {
        return FilenameUtils.getBaseName(referenceFile.getName()) + ".fasta";
    }

    /**
     * Return the reference file converter used by aligner algorithm.
     */
    @Override
    public CompactToFastaMode getReferenceCompactToFastaConverter() {
        final CompactToFastaMode processor = new CompactToFastaMode();
        processor.setHashOutputFilename(true);
        processor.setIndexToHeader(false);
        // Reference conversion (always from nt-space) *MAY* be needed by alignment algorithm to match colorspace reads platforms
        processor.setOutputColorMode(false); // BWA handles the conversion
        processor.setOutputFakeNtMode(false); // BWA handles the conversion
        processor.setOutputFakeQualityMode(false); // BWA handles the conversion
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet("ACTG"); // BWA handles the conversion
        processor.setOutputFormat(CompactToFastaMode.OutputFormat.FASTA);
        return processor;
    }


    /**
     * Return the reads file converter used by aligner algorithm.
     */
    @Override
    public CompactToFastaMode getReadsCompactToFastaConverter() {
        final CompactToFastaMode processor = new CompactToFastaMode();
        processor.setIndexToHeader(true);
        processor.setHashOutputFilename(true);
        // Since colorSpace is determined by reads platform, colorspace conversion is never needed for reads
        processor.setOutputFakeNtMode(colorSpace); // BWA uses fakeNt to represent colorSpace data
        processor.setOutputFakeQualityMode(colorSpace); // BWA uses fakeNt to represent colorSpace data
        // BWA reads FASTQ format:
        processor.setOutputFormat(CompactToFastaMode.OutputFormat.FASTQ);
        // Filter using the default alphabet for specified output mode
        processor.setAlphabet("ACTG"); // BWA expects fake-nt or nt, either way, "ACTG" is the alphabet
        processor.setTrimAdaptorLength((colorSpace) ? 2 : 0); // BWA's solid2fastq remove 2 bases - prevent aligner from finding false matches with adaptor at the start of the sequence

        return processor;
    }


    /**
     * Returns {@link edu.cornell.med.icb.goby.modes.SAMToCompactMode} processor, initialized
     * with correct input file.
     */
    @Override
    public AbstractAlignmentToCompactMode getNativeAlignmentToCompactMode(final String outputBasename) {
        // can not initialize, unless correct input file exists
        assert samBinaryFilename != null : "Can not initialize, unless SAM Binary input file exists";
        assert new File(samBinaryFilename).exists() : "Can not initialize, unless correct SAM Binary file exists";
        final SAMToCompactMode processor = new SAMToCompactMode();
        processor.setInputFile(samBinaryFilename);
        processor.setNumberOfReads(numberOfReads);
        processor.setThirdPartyInput(false);
        return processor;
    }

    /**
     * Parse and validate aligner specific options.
     * Input options are comma separated, with the syntax key1=value1,key2=value2.
     * "alignerOptions" format should match aligner's command-line specification e.g. "-key1 value1 -key2 value2"
     * This method is declared abstract so that aligners have control over the parsing of their arguments.
     */
    @Override
    public void setAlignerOptions(final String options) {
        seedLength = DEFAULT_SEED_LENGTH;
        alignerOptions = "";
        if (options == null) {
            return;
        }
        final String[] opts = options.split("[,=]");

        // seedLength parameter
        if (CLI.isKeywordGiven(opts, "l")) {
            seedLength = CLI.getIntOption(opts, "l", DEFAULT_SEED_LENGTH);
        }

        alignerOptions += (CLI.isKeywordGiven(opts, "n")) ? " -n " + CLI.getIntOption(opts, "n", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "o")) ? " -o " + CLI.getIntOption(opts, "o", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "e")) ? " -e " + CLI.getIntOption(opts, "e", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "d")) ? " -d " + CLI.getIntOption(opts, "d", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "i")) ? " -i " + CLI.getIntOption(opts, "i", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "k")) ? " -k " + CLI.getIntOption(opts, "k", -999) : "";
        //alignerOptions += (CLI.isKeywordGiven(opts, "l")) ? " -l " + CLI.getIntOption(opts, "l", -999) : ""; // seedLength
    
        alignerOptions += (CLI.isKeywordGiven(opts, "t")) ? " -t " + CLI.getIntOption(opts, "t", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "M")) ? " -M " + CLI.getIntOption(opts, "M", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "O")) ? " -O " + CLI.getIntOption(opts, "O", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "E")) ? " -E " + CLI.getIntOption(opts, "E", -999) : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "c")) ? " -c " : "";
        alignerOptions += (CLI.isKeywordGiven(opts, "N")) ? " -N " : "";
        //alignerOptions += (CLI.isKeywordGiven(opts, "R")) ? " -R " : ""; // R param is only applicable for paired-end read alignment
        assert (!alignerOptions.contains("-999")) : "Parsing error.";
    }


    public File[] indexReference(final File referenceFileOrDbBasename) throws IOException {
        databasePrefix = referenceFileOrDbBasename.toString();
        if (isDatabaseBasename(databasePrefix)) {
            return matchingExtension(referenceFileOrDbBasename.toString(), extensions);
        }
        databasePrefix = getDatabasePath(databaseDirectory, databaseName);
        if (isDatabaseBasename(databasePrefix)) {
            return matchingExtension(databasePrefix, extensions);
        }
        final File fastaReferenceFile = prepareReference(referenceFileOrDbBasename);

        forceMakeParentDir(databasePrefix);

        final String command = FilenameUtils.concat(pathToExecutables, BWA_EXEC);
        final CommandLine commandLine = createCommandLine(command);
        commandLine.addArgument(BWA_INDEX_ARG);
        commandLine.addArgument("-a");
        commandLine.addArgument(databaseIndexingStyle(fastaReferenceFile));
        if (colorSpace) {
            commandLine.addArgument("-c");
        }
        commandLine.addArgument("-p");
        commandLine.addArgument(databasePrefix);
        commandLine.addArgument(fastaReferenceFile.toString());

        LOG.info("About to execute " + commandLine);
        final StopWatch timer = new StopWatch();
        timer.start();
        final DefaultExecutor executor = new DefaultExecutor();
        OutputStream logStream = null;
        try {
            logStream = new LoggingOutputStream(getClass(), Level.INFO, "");
            executor.setStreamHandler(new PumpStreamHandler(logStream));

            final int exitValue = executor.execute(commandLine);
            LOG.info("Exit value = " + exitValue);
        } finally {
            IOUtils.closeQuietly(logStream);
        }
        LOG.info("Command executed in: " + timer.toString());

        listFiles(databaseDirectory);
        return matchingExtension(databasePrefix, extensions);
    }

    private String databaseIndexingStyle(final File fastaReferenceFile) {
        // use bwtsw for genome fasta files larger than 5M, if otherwise.
        return fastaReferenceFile.length() < 5 * 1024 * 1024 ? "is" : "bwtsw";
    }

    /**
     * Return string option of form "-l <seedLength>".
     * <p/>
     * If a positive seedLength has been initialized, use this value. Otherwise,
     * use -l 35 to restrict the seed to the first 35 bp of each read. This
     * considerably speeds up searches with ~100 bp reads.
     * <p/>
     * Only use this option if the seedLength is shorter than the minimum read length,
     * otherwise bwa may fail (despite what the documentation says).
     *
     * @return A string of the form "-l <seedLength>" or an empty string.
     */
    private String getLOption() {
        return (minReadLength >= seedLength) ? ("-l " + seedLength) : "";
    }

    public File[] align(final File referenceFile, final File readsFile, final String outputBasename) throws InterruptedException, IOException {
        assert pathToExecutables != null : "path to executables must be defined.";

        final File nativeReads = prepareReads(readsFile);
        if (nativeReads.toString().contains(" ")) {
            throw new IllegalArgumentException("BWA does not support spaces in filename: "
                    + nativeReads);
        }
        if (databaseName == null) {
            databaseName = getDefaultDbNameForReferenceFile(referenceFile);
        }
        indexReference(referenceFile);
        saiBinaryFilename = FilenameUtils.concat(workDirectory,
                File.createTempFile(RandomStringUtils.randomAlphabetic(10), ".sai").getName());
        samBinaryFilename = FilenameUtils.concat(workDirectory,
                File.createTempFile(RandomStringUtils.randomAlphabetic(10), ".sam").getName());
        if (LOG.isDebugEnabled()) {
            LOG.debug("sai file: " + saiBinaryFilename);
            LOG.debug("sam file: " + samBinaryFilename);
        }
        forceMakeParentDir(saiBinaryFilename);
        forceMakeParentDir(samBinaryFilename);
        LOG.info("Searching.");

        final String alignCommand = FilenameUtils.concat(pathToExecutables, BWA_EXEC);
        final CommandLine alignCommandLine = createCommandLine(alignCommand);
        alignCommandLine.addArgument(BWA_ALIGN_ARG);
        if (colorSpace) {
            alignCommandLine.addArgument("-c");
        }
        alignCommandLine.addArguments(getLOption());
        alignCommandLine.addArguments(alignerOptions);
        alignCommandLine.addArgument(databasePrefix);
        alignCommandLine.addArgument(nativeReads.toString());

        LOG.info("About to execute " + alignCommandLine);
        final StopWatch timer = new StopWatch();
        timer.start();
        final DefaultExecutor alignExecutor = new DefaultExecutor();
        OutputStream saiOutputStream = null;
        OutputStream alignLogStream = null;
        try {
            saiOutputStream = new BufferedOutputStream(new FileOutputStream(saiBinaryFilename));
            alignLogStream = new LoggingOutputStream(getClass(), Level.INFO, "");
            alignExecutor.setStreamHandler(new PumpStreamHandler(saiOutputStream, alignLogStream));

            final int exitValue = alignExecutor.execute(alignCommandLine);
            LOG.info("Exit value = " + exitValue);
        } finally {
            IOUtils.closeQuietly(alignLogStream);
            IOUtils.closeQuietly(saiOutputStream);
        }
        LOG.info("Command executed in: " + timer.toString());

        // convert sai to SAM format
        final String samseCommand = FilenameUtils.concat(pathToExecutables, BWA_EXEC);
        final CommandLine samseCommandLine = createCommandLine(samseCommand);
        samseCommandLine.addArgument(BWA_SAMSE_ARG);
        samseCommandLine.addArgument(databasePrefix);
        samseCommandLine.addArgument(saiBinaryFilename);
        samseCommandLine.addArgument(nativeReads.toString());

        LOG.info("About to execute " + alignCommandLine);
        timer.reset();
        timer.start();
        final DefaultExecutor samseExecutor = new DefaultExecutor();
        OutputStream samseOutputStream = null;
        OutputStream samseLogStream = null;
        try {
            samseOutputStream = new BufferedOutputStream(new FileOutputStream(samBinaryFilename));
            samseLogStream = new LoggingOutputStream(getClass(), Level.INFO, "");
            samseExecutor.setStreamHandler(new PumpStreamHandler(samseOutputStream, samseLogStream));

            final int exitValue = samseExecutor.execute(samseCommandLine);
            LOG.info("Exit value = " + exitValue);
        } finally {
            IOUtils.closeQuietly(samseLogStream);
            IOUtils.closeQuietly(samseOutputStream);
        }
        LOG.info("Command executed in: " + timer.toString());

        // convert native alignment into compact reads
        final File[] buildResults = processAlignment(referenceFile, readsFile, outputBasename);
        if (!keepTemporaryFiles) {
            FileUtils.deleteQuietly(new File(saiBinaryFilename));
            FileUtils.deleteQuietly(new File(samBinaryFilename));
        }
        return buildResults;
    }

    public boolean isNativeGobySupported() {
        return NATIVE_GOBY_SUPPORTED;
    }

}
