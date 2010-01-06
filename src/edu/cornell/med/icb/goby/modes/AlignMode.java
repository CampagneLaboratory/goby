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
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.aligners.Aligner;
import edu.cornell.med.icb.goby.aligners.BWAAligner;
import edu.cornell.med.icb.goby.aligners.LastagAligner;
import edu.cornell.med.icb.goby.config.ConfigHelper;
import edu.cornell.med.icb.goby.config.GobyPropertyKeys;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.File;
import java.io.IOException;

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
    private static final Log LOG = LogFactory.getLog(AlignMode.class);

    public enum AlignerTypes {
        lastag, bwa
    }

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "align";
    public static final String MODE_DESCRIPTION = "Run third party aligners, taking care of data translations.  Data translations include converting compact reads to the aligner input and converting the aligner output to compact alignment format.";

    public static final int DEFAULT_M_PARAMETER = 2;

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
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        readsFile                  = jsapResult.getFile    ( "reads"                     ) ;
        referenceFile              = jsapResult.getFile    ( "reference"                 ) ;
        readIndexFilterFile        = jsapResult.getFile    ( "read-index-filter"         ) ;
        referenceIndexFilterFile   = jsapResult.getFile    ( "reference-index-filter"    ) ;

        workDirectory              = jsapResult.getFile    ( "work-directory"            ) ;
        databaseDirectory          = jsapResult.getFile    ( "database-directory"        ) ;
        databasePath               = jsapResult.getString  ( "database-name"             ) ;

        indexFlag                  = jsapResult.getBoolean ( "index"                     ) ;
        searchFlag                 = jsapResult.getBoolean ( "search"                    ) ;
        getDefaultDatabaseNameFlag = jsapResult.getBoolean ( "get-default-database-name" ) ;

        alignerName                = jsapResult.getString  ( "aligner"                   ) ;
        alignerOptions             = jsapResult.getString  ( "options"                   ) ;

        colorSpace                 = jsapResult.getBoolean ( "color-space"               ) ;
        // TODO 090909 replace with setFilterOptions()
        qualityFilterParameters    = jsapResult.getString  ( "quality-filter-parameters" ) ;
        mParameter                 = jsapResult.getInt     ( "ambiguity-threshold"       ) ;
        outputBasename             = jsapResult.getString  ( "basename"                  ) ;

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
            case lastag:
                aligner = new LastagAligner();
                break;
            case bwa:
                aligner = new BWAAligner();
                break;
        }

        if (aligner == null) {
            System.out.println("Could not initialize aligner.");
            System.exit(1);
        }

        try {
            final Configuration configuration = ConfigHelper.getConfiguration();
            aligner.setConfiguration(configuration);
            aligner.setWorkDirectory(workDirectory != null ? workDirectory.getPath() :
                    configuration.getString(GobyPropertyKeys.WORK_DIRECTORY));
            aligner.setDatabaseDirectory(databaseDirectory != null ? databaseDirectory.getPath() :
                    configuration.getString(GobyPropertyKeys.DATABASE_DIRECTORY));
            aligner.setDatabaseName(databasePath);
            aligner.setReadIndexFilter(readIndexFilterFile);
            aligner.setReferenceIndexFilter(referenceIndexFilterFile);
            aligner.setColorSpace(colorSpace);
            aligner.setAlignerOptions(alignerOptions);
            // TODO 090909 replace with setFilterOptions()
            aligner.setQualityFilterParameters(qualityFilterParameters);
            aligner.setAmbiguityThreshold(mParameter);

            if (indexFlag) {
                aligner.indexReference(referenceFile);
            } else if (searchFlag) {
                alignmentFiles = aligner.align(referenceFile, readsFile, outputBasename);
            } else if (getDefaultDatabaseNameFlag) {
                System.out.println(aligner.getDefaultDbNameForReferenceFile(referenceFile));

            }
        } catch (InterruptedException e) {
            LOG.error(e);
            Thread.currentThread().interrupt();
        }
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new AlignMode().configure(args).execute();
    }
}
