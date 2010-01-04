package edu.cornell.med.icb.goby.util;

/**
 * Defines property keys for local configuration entries used by Goby. See ConfigHelper to retrieve
 * property values.
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
     * Path to the lastag executable.
     */
    public static final String EXECUTABLE_PATH_LASTAG = "executables.path.lastag";
    /**
     * Path to the BWA executable.
     */
    public static final String EXECUTABLE_PATH_BWA = "executables.path.bwa";
    /**
     * Path to the work directory. This should be a large scratch location, where results of intermediate
     * calculations will be stored.
     */
    public static final String WORK_DIRECTORY = "work.directory";
    /**
     * Path to the database directory. This is the location where the indexed database files may reside.
     */
    public static final String DATABASE_DIRECTORY = "database.directory";
    /**
     * The Amazon AWS Secret Access Key. Required to access S3 resources.
     */
    public static final String AWS_SECRET_ACCESS_KEY = "aws.secret.access.key";
    /**
     * The Amazon AWS Access Key. Required to access S3 resources.
     */
    public static final String AWS_ACCESS_KEY = "aws.access.key";
    /**
     * The name of the gridgain grid to use for computations.
     */
    public static final String GRID_NAME = "grid.name";
    /**
     * The maximum number of jobs to run at any given time on the gridgain node.
     */
    public static final String GRID_MAX_NODE = "grid.max.nodes";
    /**
     * The template to use to call the scp command to copy from remote node to local node.
     */
    public static final String SCP_COMMAND_TEMPLATE_REMOTE_TO_LOCAL = "scp.command.template.remote.local";
    /**
     * The template to use to call the scp command to copy from local node to remote node.
     */
    public static final String SCP_COMMAND_TEMPLATE_LOCAL_TO_REMOTE = "scp.command.template.local.remote";

    /** The key for the default s3 directory for storing s3 queue overflow files. */
    public static final String S3_QUEUE_OVERFLOW_PATH = "s3.queue.overflow.path";

       /**
     * Host name of this node. Should be publicly visible.
     */
    public static final String HOST_NAME = "host.name";

    /** The name of the input queue, in the Global/LocalConfig.groovy. */
    public static final String INPUT_QUEUE_NAME = "input.queue.name";

    /** The name of the output queue, in the Global/LocalConfig.groovy. */
    public static final String OUTPUT_QUEUE_NAME = "output.queue.name";

    /** The S3 bucket name to use. */
    public static final String S3_BUCKET_NAME = "s3.bucket.name";

    public static final String STAT_NUMBER_OF_ALIGNED_READS="number.aligned.reads";

    /** if defined in LocalConfig.groovy this specifies the JMS server, otherwise JGroups will be used. */
    public static final String JMS_SERVER="jms.server";

    /** Cache for local versions of S3 files to make things faster for "local" runs. */
    public static final String S3_LOCAL_CACHE_PATH = "s3.local.cache.path";
}
