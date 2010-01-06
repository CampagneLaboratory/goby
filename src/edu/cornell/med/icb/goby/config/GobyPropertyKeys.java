package edu.cornell.med.icb.goby.config;

/**
 * Defines property keys for local configuration entries used by Goby. See
 * {@link edu.cornell.med.icb.goby.config.ConfigHelper} to retrieve property values.
 *
 * @author Fabien Campagne
 *         Date: Aug 13, 2009
 *         Time: 12:07:23 PM
 */
public class GobyPropertyKeys {
    /**
     * Empty constructor for utility class.
     */
    private GobyPropertyKeys() {
        super();
    }

    /**
     * Path to the directory that contains the lastag executable.
     */
    public static final String EXECUTABLE_PATH_LASTAG = "executables.path.lastag";

    /**
     * Path to the directory that contains the BWA executable.
     */
    public static final String EXECUTABLE_PATH_BWA = "executables.path.bwa";

    /**
     * Path to the work directory. This should be a large scratch location, where results of
     * intermediate calculations will be stored.
     */
    public static final String WORK_DIRECTORY = "work.directory";

    /**
     * Path to the database directory. This is the location where the indexed database files may
     * reside.
     */
    public static final String DATABASE_DIRECTORY = "database.directory";
}
