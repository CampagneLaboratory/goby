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
import edu.cornell.med.icb.goby.aligners.*;
import edu.cornell.med.icb.goby.config.GobyConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Run third party aligners, taking care of data translations.  Data translations include
 * converting compact reads to the aligner input and converting the aligner output to compact
 * alignment format.
 *
 * @author Fabien Campagne
 *         Date: May 4 2009
 *         Time: 12:28 PM
 */
public class AlignMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(AlignMode.class);


    public void setKeepTemporaryFiles(final boolean keepTemporaryFiles) {
        this.keepTemporaryFiles = keepTemporaryFiles;
    }

    private boolean keepTemporaryFiles;

    /**
     * Set quality filter parameters programmatically.
     *
     * @param qualityFilterParams String of the form threshold=<double>, see
     * {@link edu.cornell.med.icb.goby.alignments.filters.AlignmentQualityFilter} for syntax
     * of optional parameters.
     */
    public void setQualityFilterParameters(final String qualityFilterParams) {
        this.qualityFilterParameters = qualityFilterParams;
    }

    /**
     * Supported native aligner types.
     */
    public enum AlignerTypes {  // TODO - move to aligner package and make factory class
        /**
         * Burrows-Wheeler Aligner (BWA).
         *
         * @see edu.cornell.med.icb.goby.aligners.BWAAligner
         * @see <a href="http://bio-bwa.sourceforge.net/">http://bio-bwa.sourceforge.net/</a>
         */
        bwa,

        /**
         * Last.
         *
         * @see edu.cornell.med.icb.goby.aligners.LastAligner
         * @see <a href="http://last.cbrc.jp/">http://last.cbrc.jp/</a>
         */
        last,

        /**
         * Enhanced version of the last aligner.
         *
         * @see edu.cornell.med.icb.goby.aligners.LastagAligner
         */
        lastag,
        /**
         * GSNAP.
         *
         * @see edu.cornell.med.icb.goby.aligners.GSnapAligner
         * @see <a href="http://research-pub.gene.com/gmap/">http://research-pub.gene.com/gmap//</a>
         */
        gsnap
    }

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "align";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Run third party aligners, taking care of data "
            + "translations.  Data translations include converting compact reads to the aligner "
            + "input and converting the aligner output to compact alignment format.";

    private static final int DEFAULT_M_PARAMETER = 2;

    private File readsFile;
    private File referenceFile;
    private File readIndexFilterFile;
    private File referenceIndexFilterFile;

    private File workDirectory;
    private File databaseDirectory;
    private String databasePath;

    private boolean getDefaultDatabaseNameFlag;
    private boolean indexFlag;
    private boolean searchFlag;

    private String alignerName;
    private AlignerTypes alignerType;
    private String alignerOptions;

    private boolean colorSpace;
    // TODO 090909 replace with setFilterOptions()
    private String qualityFilterParameters;
    private int mParameter = DEFAULT_M_PARAMETER;
    private String outputBasename;

    private long splitStartPosition = 0;
    private long splitEndPosition = 0;
    private boolean pairedEndCompactInput = false;
    private String pairedEndDirections = "FR";
    private boolean bisulfiteInput = false;

    private File[] alignmentFiles;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    public void setAlignerType(final AlignerTypes alignerType) {
        this.alignerType = alignerType;
    }

    public void setColorSpace(final boolean colorSpace) {
        this.colorSpace = colorSpace;
    }

    public void setSearchFlag(final boolean searchFlag) {
        this.searchFlag = searchFlag;
    }

    public void setIndexFlag(final boolean indexFlag) {
        this.indexFlag = indexFlag;
    }

    public void setAlignerOptions(final String alignerOptions) {
        this.alignerOptions = alignerOptions;
    }

    public void setAlignerName(final String alignerName) {
        this.alignerType = AlignerTypes.valueOf(alignerName);
        this.alignerName = alignerName;
    }

    public void setDatabaseDirectory(final File databaseDirectory) {
        this.databaseDirectory = databaseDirectory;
    }

    public void setWorkDirectory(final File workDirectory) {
        this.workDirectory = workDirectory;
    }

    public void setReferenceFile(final File referenceFile) {
        this.referenceFile = referenceFile;
    }

    public void setOutputBasename(final String outputBasename) {
        this.outputBasename = outputBasename;
    }

    public void setReadsFile(final File readsFile) {
        this.readsFile = readsFile;
    }

    public File[] getAlignmentFiles() {
        return alignmentFiles;
    }

    public void setDatabasePath(final String databasePath) {
        this.databasePath = databasePath;
    }

    public void setReadIndexFilterFile(final File readIndexFilterFile) {
        this.readIndexFilterFile = readIndexFilterFile;
    }

    public void setReferenceIndexFilterFile(final File referenceIndexFilterFile) {
        this.referenceIndexFilterFile = referenceIndexFilterFile;
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
        final Map<String, String> alignerTypesHelp = new HashMap<String, String>(1);
        alignerTypesHelp.put("[ALIGNERTYPES]", ArrayUtils.toString(AlignerTypes.values()));
        final JSAPResult jsapResult = parseJsapArguments(args, alignerTypesHelp);

        readsFile = jsapResult.getFile("reads");
        referenceFile = jsapResult.getFile("reference");
        readIndexFilterFile = jsapResult.getFile("read-index-filter");
        referenceIndexFilterFile = jsapResult.getFile("reference-index-filter");

        workDirectory = jsapResult.getFile("work-directory");
        databaseDirectory = jsapResult.getFile("database-directory");
        databasePath = jsapResult.getString("database-name");

        indexFlag = jsapResult.getBoolean("index");
        searchFlag = jsapResult.getBoolean("search");
        getDefaultDatabaseNameFlag = jsapResult.getBoolean("get-default-database-name");

        alignerName = jsapResult.getString("aligner");
        alignerOptions = jsapResult.getString("options");

        colorSpace = jsapResult.getBoolean("color-space");
        keepTemporaryFiles = jsapResult.getBoolean("keep-temporary-files");
        // TODO 090909 replace with setFilterOptions()
        qualityFilterParameters = jsapResult.getString("quality-filter-parameters");
        mParameter = jsapResult.getInt("ambiguity-threshold");

        // Zero should indicate start of file
        splitStartPosition = jsapResult.getLong("split-start-position");
        // Zero should indicate end of file
        splitEndPosition = jsapResult.getLong("split-end-position");
        pairedEndCompactInput = jsapResult.getBoolean("paired-end-compact-input");
        try {
            setPairedEndDirections(jsapResult.getString("paired-end-directions"));
        } catch (IllegalArgumentException e) {
            System.err.println("--paired-end-directions value supported: " +
                    jsapResult.getString("paired-end-directions"));
            System.exit(1);
        }
        bisulfiteInput = jsapResult.getBoolean("bisulfite-input");

        outputBasename = jsapResult.getString("basename");

        try {
            alignerType = AlignerTypes.valueOf(alignerName.toLowerCase());
        } catch (IllegalArgumentException e) {
            System.err.println("Aligner type is not supported: " + alignerName);
            System.exit(1);
        }
        if (!(indexFlag || searchFlag || getDefaultDatabaseNameFlag)) {
            System.err.println("One of --index or --search or --get-default-database-name must be specified.");
            System.exit(1);
        }
        if (outputBasename == null && searchFlag) {
            System.err.println("In --search mode --basename is required.");
            System.exit(1);
        }
        if (databasePath == null && indexFlag) {
            System.err.println("In --index mode --database-name is required.");
            System.exit(1);
        }
        return this;
    }

    @Override
    public void execute() throws IOException {
        Aligner aligner = null;
        switch (alignerType) {
            case last:
                aligner = new LastAligner();
                break;
            case lastag:
                aligner = new LastagAligner();
                break;
            case bwa:
                aligner = new BWAAligner();
                break;

            case gsnap:
                aligner = new GSnapAligner();
                break;
            default:
                System.err.println("Unsupported aligner: " + alignerType);
                System.exit(2);
        }

        if (aligner == null) {
            System.err.println("Could not initialize aligner.");
            System.exit(1);
        }

        try {
            final Configuration configuration = GobyConfiguration.getConfiguration();
            aligner.setConfiguration(configuration);
            aligner.setWorkDirectory(workDirectory != null ? workDirectory.getPath()
                    : configuration.getString(GobyConfiguration.WORK_DIRECTORY));
            aligner.setDatabaseDirectory(databaseDirectory != null ? databaseDirectory.getPath()
                    : configuration.getString(GobyConfiguration.DATABASE_DIRECTORY));
            aligner.setDatabaseName(databasePath);
            aligner.setReadIndexFilter(readIndexFilterFile);
            aligner.setReferenceIndexFilter(referenceIndexFilterFile);
            aligner.setColorSpace(colorSpace);
            aligner.setAlignerOptions(alignerOptions);
            // TODO 090909 replace with setFilterOptions()
            aligner.setQualityFilterParameters(qualityFilterParameters);
            aligner.setAmbiguityThreshold(mParameter);
            aligner.setKeepTemporaryFiles(keepTemporaryFiles);
            aligner.setSplitStartPosition(splitStartPosition);
            aligner.setSplitEndPosition(splitEndPosition);
            aligner.setPairedEndCompactInput(pairedEndCompactInput);
            aligner.setPairedEndDirections(pairedEndDirections);
            aligner.setBisulfiteInput(bisulfiteInput);

            if (indexFlag) {
                aligner.indexReference(referenceFile);
            } else if (searchFlag) {
                alignmentFiles = aligner.align(referenceFile, readsFile, outputBasename);
            } else if (getDefaultDatabaseNameFlag) {
                System.out.println(aligner.getDefaultDbNameForReferenceFile(referenceFile));
            }
        } catch (InterruptedException e) {
            LOG.error("Interrupted", e);
            Thread.currentThread().interrupt();
        }
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new AlignMode().configure(args).execute();
    }
}
