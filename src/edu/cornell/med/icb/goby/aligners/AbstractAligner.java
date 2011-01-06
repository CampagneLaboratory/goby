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

import edu.cornell.med.icb.goby.modes.AbstractAlignmentToCompactMode;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Locale;

/**
 * @author Fabien Campagne
 *         Date: Jul 10, 2009
 *         Time: 3:26:45 PM
 */
// TODO - replace mParameter and qualityFilterParameters with filterParams?
public abstract class AbstractAligner implements Aligner {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AbstractAligner.class);

    protected String pathToExecutables;
    protected String databaseDirectory = ".";
    protected String workDirectory = ".";
    protected String[] extensions;
    protected String databaseName;
    protected int numberOfReads;
    protected int numberOfReferences;
    protected File readIndexFilter;
    protected File referenceIndexFilter;
    protected boolean colorSpace;
    protected int minReadLength;
    protected String alignerOptions = "";
    protected String qualityFilterParameters = "threshold=0.05";
    protected int mParameter = 2;
    protected boolean keepTemporaryFiles;
    private int smallestQueryIndex;
    private int largestQueryIndex;
    protected long splitStartPosition = 0;
    protected long splitEndPosition = 0;
    protected boolean pairedEndCompactInput = false;
    protected String pairedEndDirections = "FR";
    protected boolean bisulfiteInput = false;

    /**
     * When true, will not delete the native aligner files generated during alignment. These files are deleted by default.
     *
     * @param keepTemporaryFiles
     */
    public void setKeepTemporaryFiles(final boolean keepTemporaryFiles) {
        this.keepTemporaryFiles = keepTemporaryFiles;
    }

    /**
     * Return the reference file converter used by aligner algorithm.
     */
    public abstract CompactToFastaMode getReferenceCompactToFastaConverter();

    /**
     * Return the reads file converter used by aligner algorithm.
     */
    public abstract CompactToFastaMode getReadsCompactToFastaConverter();

    /**
     * Return the compact reads converter for the native alignment.
     */
    public abstract AbstractAlignmentToCompactMode getNativeAlignmentToCompactMode(final String outputBasename);


    public abstract boolean isNativeGobySupported();

    /**
     * Parse and validate aligner specific options.
     * Input options are comma separated, with the syntax key1=value1,key2=value2.
     * "alignerOptions" format should match aligner's command-line specification e.g. "-key1 value1 -key2 value2"
     * This method is declared abstract so that aligners have control over the parsing of their arguments.
     */
    public abstract void setAlignerOptions(final String alignerOptions);

    public void setDatabaseName(final String databaseName) {
        this.databaseName = databaseName;
    }

    public String getDatabaseName() {
        return databaseName;
    }

    public void setReferenceIndexFilter(final File referenceIndexFilterFile) {
        referenceIndexFilter = referenceIndexFilterFile;
    }

    public void setReadIndexFilter(final File readIndexFilterFile) {
        readIndexFilter = readIndexFilterFile;
    }

    public void setColorSpace(final boolean colorSpace) {
        this.colorSpace = colorSpace;
    }

    public String getAlphabet() {
        return colorSpace ? "0123" : "ACTG";
    }

    public void setQualityFilterParameters(final String qualityFilterParameters) {
        this.qualityFilterParameters = qualityFilterParameters;
    }

    public void setAmbiguityThreshold(final int mParameter) {
        this.mParameter = mParameter;
    }

    public void setDatabaseDirectory(final String path) {
        databaseDirectory = path;
    }

    public void setPathToExecutables(final String path) {
        pathToExecutables = path;
    }

    public void setWorkDirectory(final String path) {
        workDirectory = path;
    }

    /**
     * If the aligner directly supports reading Goby Compact-Reads (check the aligners isNativeGobySupport()),
     * this sets the start position when reading from the input file. The default, 0, indicates start of file.
     * @return the split start position
     */
    public long getSplitStartPosition() {
        return splitStartPosition;
    }

    /**
     * If the aligner directly supports reading Goby Compact-Reads (check the aligners isNativeGobySupport()),
     * this sets the start position when reading from the input file. The default, 0, indicates start of file.
     * @param splitStartPosition the split start position
     */
    public void setSplitStartPosition(long splitStartPosition) {
        this.splitStartPosition = splitStartPosition;
    }

    /**
     * If the aligner directly supports reading Goby Compact-Reads (check the aligners isNativeGobySupport()),
     * this sets the end position when reading from the input file. The default, 0, indicates end of file.
     * @return the split end position
     */
    public long getSplitEndPosition() {
        return splitEndPosition;
    }

    /**
     * If the aligner directly supports reading Goby Compact-Reads (check the aligners isNativeGobySupport()),
     * this sets the end position when reading from the input file. The default, 0, indicates end of file.
     * @param splitEndPosition the split end position
     */
    public void setSplitEndPosition(long splitEndPosition) {
        this.splitEndPosition = splitEndPosition;
    }

    /**
     * If the aligner directly supports reading Goby Compact-Reads, this sets if the input Goby compact-reads
     * file is paired-end and should be processed as such, if the aligners supports it.
     * @return if the input is paired end Goby compact-reads
     */
    public boolean isPairedEndCompactInput() {
        return pairedEndCompactInput;
    }

    /**
     * If the aligner directly supports reading Goby Compact-Reads, this sets if the input Goby compact-reads
     * file is paired-end and should be processed as such, if the aligners supports it.
     * @param pairedEndCompactInput if the input is paired end Goby compact-reads
     */
    public void setPairedEndCompactInput(boolean pairedEndCompactInput) {
        this.pairedEndCompactInput = pairedEndCompactInput;
    }

    /**
     * If pairedEndCompactInput is set to true, this specifies the directions for the sequence pair.
     * "FF", "FR", "RF", "RR" and valid values. Default is FR.
     * @return the paired end directions if the input is paired end Goby compact-reads
     */
    public String getPairedEndDirections() {
        return pairedEndDirections;
    }

    /**
     * If pairedEndCompactInput is set to true, this specifies the directions for the sequence pair.
     * "FF", "FR", "RF", "RR" and valid values. If you try to set to a value other than these four
     * nothing will happen. Default is "FR".
     * @param pairedEndDirections the paired end directions if the input is paired end Goby compact-reads
     */
    public void setPairedEndDirections(String pairedEndDirections) throws IllegalArgumentException {
        if (pairedEndDirections.equals("FF") ||
                pairedEndDirections.equals("FR") ||
                pairedEndDirections.equals("RF") ||
                pairedEndDirections.equals("RR")) {
            this.pairedEndDirections = pairedEndDirections;
        } else {
            throw new IllegalArgumentException();
        }
    }

    /**
     * If the input has be bisulfite processed (currently only supported by gsnap and requires a bisulfite reference).
     * @return if Bisulfite Input
     */
    public boolean isBisulfiteInput() {
        return bisulfiteInput;
    }

    /**
     * If the input has be bisulfite processed (currently only supported by gsnap and requires a bisulfite reference).
     * @param bisulfiteInput if Bisulfite Input
     */
    public void setBisulfiteInput(boolean bisulfiteInput) {
        this.bisulfiteInput = bisulfiteInput;
    }

    /**
     * If databasePrefix is a complete basename to a database, this will return databasePrefix.
     * Otherwise this will assume databasePrefix is actually a database-name
     * and return databaseDirectory + "/" + databasePrefix.
     *
     * @param databaseDirectory the database directory
     * @param databasePrefix    the database basename or database-name.
     * @return a database basename
     */
    protected String getDatabasePath(final String databaseDirectory, final String databasePrefix) {
        final String result;
        if (isDatabaseBasename(databasePrefix)) {
            result = databasePrefix;
        } else {
            if (StringUtils.isBlank(databaseDirectory)) {
                result = databasePrefix;
            } else {
                result = FilenameUtils.concat(databaseDirectory, databasePrefix);
            }
        }
        return result.replaceAll(" ", "\\ ");

    }

    /**
     * Returns true if the reference is a database basename.
     * Searches for databasePathPrefix + "." + [extensions] to make sure all those
     * files exist.
     *
     * @param databasePathPrefix path prefix for database files
     * @return true if databasePathPrefix specifies an existing database
     */
    public boolean isDatabaseBasename(final String databasePathPrefix) {
        assert extensions != null : "aligner must define valid database extensions. ";
        for (final String ext : extensions) {
            final String inspectionPath = databasePathPrefix + "." + ext;
            System.err.println("Looking for database... file exists?:" + inspectionPath);
            if (!new File(inspectionPath).exists()) {
                System.err.println("...No");
                return false;
            }
        }
        System.err.println("...Yes");
        return true;
    }

    protected void listFiles(final String directory) {
        assert StringUtils.isNotBlank(directory) : "Directory must not be blank or null";
        if (LOG.isDebugEnabled()) {
            LOG.debug("Files in directory: " + directory);
            final Collection databaseFiles = FileUtils.listFiles(new File(directory), null, false);
            for (final Object databaseFile : databaseFiles) {
                LOG.debug(databaseFile);
            }
        }
    }

    /**
     * Collect the set of files that make up the alignment output of the aligner.
     *
     * @param basename Basename where the alignment was written.
     * @return an array of files that will contain the output from the aligner
     */
    public static File[] buildResults(final String basename) {
        final ObjectArrayList<File> result = new ObjectArrayList<File>();
        final String[] extensions = {".stats", ".entries", ".header", ".tmh"};

        for (final String extension : extensions) {
            final File file = new File(basename + extension);
            if (file.exists()) {
                result.add(file);
            }
        }

        return result.toArray(new File[result.size()]);
    }

    protected File[] matchingExtension(final String filenamePrefix, final String... extension) {
        final ObjectList<File> result = new ObjectArrayList<File>();
        for (final String ext : extension) {
            final File potential = new File(filenamePrefix + "." + ext);
            if (potential.exists()) {
                result.add(potential);
            }
        }
        return result.toArray(new File[result.size()]);
    }

    /**
     * Ensure *parent* directories exist.
     */
    protected void forceMakeParentDir(final String path) {
        final String parentDir = FilenameUtils.getFullPath(path);
        forceMakeDir(parentDir);
    }

    /**
     * Ensure directories exist.
     */
    private void forceMakeDir(final String path) {
        final File dir = new File(path);
        if (path.length() == 0) {
            // TODO: add unit test with various working directories, including "./" case
            // file is perhaps in current directory?
            return;
        }
        try {
            // once again, forceMkdir will sometimes fail if the directory exists. Do not believe its javadoc.
            if (!dir.exists()) {
                FileUtils.forceMkdir(dir);
            }
        } catch (IOException e) {
            LOG.warn("Error trying to create directory. Trying to recover.", e);
        }
        assert dir.exists() : "Abstract Aligner> can not proceed without directory: '" + path + "'";
    }

    public File prepareReference(final File compactReferenceFile) throws IOException {
        LOG.info("Preparing reference sequences");
        forceMakeDir(databaseDirectory);

        final String extension = FilenameUtils.getExtension(compactReferenceFile.getName());
        if ("fasta".equals(extension)) {
            final File outputFile = new File(FilenameUtils.concat(databaseDirectory,
                    FilenameUtils.getName(compactReferenceFile.toString())));
            if (compactReferenceFile.compareTo(outputFile) != 0) {
                FileUtils.copyFile(compactReferenceFile, outputFile);
            }
            return outputFile;
        } else if ("fa".equals(extension)) {
            final File outputFile = new File(FilenameUtils.concat(databaseDirectory,
                    FilenameUtils.getBaseName(compactReferenceFile.toString()) + ".fasta"));
            FileUtils.copyFile(compactReferenceFile, outputFile);
            return outputFile;
        } else {
            // assume compact-reads format:
            // TODO: assert that file is a compact-reads file
            final String outputFilename = FilenameUtils.concat(databaseDirectory,
                    FilenameUtils.getBaseName(compactReferenceFile.getCanonicalPath()) + ".fasta");
            final CompactToFastaMode processor = getReferenceCompactToFastaConverter();
            processor.setInputFilename(compactReferenceFile.getCanonicalPath());
            processor.setReadIndexFilterFile(referenceIndexFilter);
            processor.setOutputFilename(outputFilename);
            processor.execute();
            numberOfReferences = processor.getNumberOfSequences();
            return new File(processor.getOutputFilename());
        }
    }

    public File prepareReads(final File compactReadsFile) throws IOException {
        LOG.info("Preparing reads..");
        forceMakeDir(workDirectory);

        final String extension = FilenameUtils.getExtension(compactReadsFile.getName());
        if ("fasta".equals(extension)) {
            final File outputFile = new File(FilenameUtils.concat(workDirectory,
                    FilenameUtils.getName(compactReadsFile.toString())));
            if (compactReadsFile.compareTo(outputFile) != 0) {
                FileUtils.copyFile(compactReadsFile, outputFile);
            }
            return outputFile;
        } else if ("fa".equals(extension)) {
            final File outputFile = new File(FilenameUtils.concat(workDirectory,
                    FilenameUtils.getBaseName(compactReadsFile.toString()) + ".fasta"));
            FileUtils.copyFile(compactReadsFile, outputFile);
            return outputFile;
        } else {
            // assume compact-reads format:
            // TODO: assert that file is a compact-reads file
            final CompactToFastaMode processor = getReadsCompactToFastaConverter();
            processor.setInputFilename(compactReadsFile.getCanonicalPath());
            processor.setReadIndexFilterFile(readIndexFilter);

            // output file name is based on the input file and type
            final String outputExtension = "."
                   +  processor.getOutputFormat().toString().toLowerCase(Locale.getDefault());
            final String outputFilename = FilenameUtils.concat(workDirectory,
                    FilenameUtils.getBaseName(compactReadsFile.getCanonicalPath())
                            + outputExtension);
            processor.setOutputFilename(outputFilename);
            processor.execute();
            numberOfReads = processor.getNumberOfSequences();
            minReadLength = processor.getMinSequenceLength();
            smallestQueryIndex = processor.getSmallestQueryIndex();
            largestQueryIndex = processor.getLargestQueryIndex();
            return new File(processor.getOutputFilename());
        }
    }

    public File[] processAlignment(final File referenceFile, final File readsFile, final String outputBasename) throws IOException {
        forceMakeParentDir(outputBasename);
        final AbstractAlignmentToCompactMode processor = getNativeAlignmentToCompactMode(outputBasename);

        /*
        // If match quality is BEST_MATCH, then the quality filter is turned off.
        qualityFilterParameters = (matchQuality == MatchQuality.BEST_MATCH) ? "threshold=0.00" : qualityFilterParameters;
        */
        processor.setOutputFile(outputBasename);
        processor.setPropagateQueryIds(false);
        processor.setPropagateTargetIds(true);
        processor.setTargetReferenceIdsFilename(referenceFile.getPath());
        processor.setQueryReadIdsFilename(readsFile.getPath());
        processor.setQualityFilterParameters(qualityFilterParameters);
        processor.setAmbiguityThreshold(mParameter);
        processor.setSmallestQueryIndex(smallestQueryIndex);
        processor.setLargestQueryIndex(largestQueryIndex);
        processor.setNumberOfQuerySequences(numberOfReads);
        processor.execute();

        return buildResults(outputBasename);
    }

    /**
     * Create a new {@link org.apache.commons.exec.CommandLine} prepending the
     * command with "nice" if running on a Unix system.
     *
     * @param command The path to the command to run
     * @return A command line object potential with "nice" prepended
     */
    protected CommandLine createCommandLine(final String command) {
        final CommandLine commandLine;
        if (SystemUtils.IS_OS_UNIX) {
            commandLine = CommandLine.parse("nice");
            commandLine.addArgument(command);
        } else {
            commandLine = CommandLine.parse(command);
        }
        return commandLine;
    }
}
