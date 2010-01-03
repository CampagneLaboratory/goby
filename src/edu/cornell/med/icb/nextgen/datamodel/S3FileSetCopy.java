package edu.cornell.med.icb.nextgen.datamodel;

import edu.cornell.med.icb.config.ConfigHelper;
import edu.cornell.med.icb.util.GobyPropertyKeys;
import edu.cornell.med.icb.util.GroovyProperties;
import edu.cornell.med.icb.util.ICBFilenameUtils;
import edu.cornell.med.icb.util.LoggingEqualsBuilder;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.jets3t.service.S3ServiceException;
import org.jets3t.service.impl.rest.httpclient.RestS3Service;
import org.jets3t.service.model.S3Bucket;
import org.jets3t.service.model.S3Object;
import org.jets3t.service.multithread.DownloadPackage;
import org.jets3t.service.multithread.S3ServiceSimpleMulti;
import org.jets3t.service.security.AWSCredentials;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Fabien Campagne
 *         Date: Jul 30, 2009
 *         Time: 6:01:05 PM
 */
public class S3FileSetCopy implements FileSetCopy, Serializable {

    private static final long serialVersionUID = 276879213023981460L;

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(S3FileSetCopy.class);

    private String url;
    private String fileType;
    private String bucket;
    private String key;

    /**
     * Size is unknown initially.
     */
    long size = -1;
    /**
     * Timestamp is unknown initially.
     */
    long timeStamp = -1;

    /**
     * The names of the files stored under this S3 location/directory.
     */
    final ObjectArrayList<String> filenames = new ObjectArrayList<String>();
    private boolean resourceIsFileList;

    /**
     * Pattern to parse bucket, key from s3 url.
     */
    private static final Pattern S3_URL_PATTHERN = Pattern.compile("^s3:/{1,2}(.*?)/(.*)$");

    /**
     * Empty constructor is required by serialization mechanism.
     */
    public S3FileSetCopy() {
    }

    public static boolean isS3Url(final String url) {
        final Matcher matcher = S3_URL_PATTHERN.matcher(url);
        return matcher.matches();
    }

    public S3FileSetCopy(final String url, final String fileType) {
        assert url.startsWith("s3:/");
        final Matcher matcher = S3_URL_PATTHERN.matcher(url);
        final boolean matches = matcher.matches();
        assert matches;
        bucket = matcher.group(1);
        key = matcher.group(2);
        this.url = url;
        this.fileType = fileType;
        if (fileType.equals(Resource.FILE_LIST)) {
            resourceIsFileList = true;
        }
    }

    public boolean isLocal() {
        return false;
    }

    public void delete() {
    }

    /**
     * Downloads files from S3 and returns the set of local files.
     *
     * @param fileType
     * @return
     */
    public File[] getFiles(final String... fileType) {
        try {
            final GroovyProperties config = ConfigHelper.loadConfiguration();
            final String s3LocalCacheDir = config.get(GobyPropertyKeys.S3_LOCAL_CACHE_PATH);

            final RestS3Service service = initializeS3Connect();

            final S3Object[] content;
            /**
             * Directory that will hold the local file set:
             */
            final File dir = new File(createLocalTempDirectory());

            final S3ServiceSimpleMulti multi = new S3ServiceSimpleMulti(service);

            final S3Bucket s3Bucket = new S3Bucket(bucket);

            if (resourceIsFileList) {
                // first, download the list to the local directory :
                final File localListFile = newLocalFile(dir, FilenameUtils.getName(key));
                final DownloadPackage[] dpArray = {
                        new DownloadPackage(new S3Object(s3Bucket, key), localListFile) };

                // copy the list to local file:
                final String transferTag = Tagged.createStringTag();
                LOG.info("Starting s3 transfer " + dpArray.length + " files, transfer tag = "
                        + transferTag);
                multi.downloadObjects(s3Bucket, dpArray);
                LOG.info("Completed s3 transfer " + dpArray.length + " files, transfer tag = "
                        + transferTag);

                // Get the parent key of the list file. This is where the filenames defined in the list exist:
                final String parentKey = new File(key).getParent();
                final ObjectArrayList<S3Object> listObjects = new ObjectArrayList<S3Object>();
                try {
                    final LineIterator lines = new LineIterator(new FastBufferedReader(new FileReader(localListFile)));
                    for (final MutableString line : lines.allLines()) {
                        final String filename = line.trim().toString();
                        LOG.trace(String.format("adding file %s to resource %s", filename, url));
                        // build the filename path:
                        listObjects.add(new S3Object(s3Bucket,
                                ICBFilenameUtils.concatPathParts(parentKey, filename)));
                    }
                    content = listObjects.toArray(new S3Object[listObjects.size()]);
                } catch (IOException e) {
                    throw new RuntimeException("Cannot open remote S3 file list at path " + url);
                }

            } else {
                // list the directory or file to obtain S3Object details:
                content = service.listObjects(s3Bucket, key, null);
            }

            final ArrayList<File> files = new ArrayList<File>();

            final ObjectArrayList<DownloadPackage> dp = new ObjectArrayList<DownloadPackage>();
            for (final S3Object object : content) {
                final long s3FileLength = s3ObjectSize(s3Bucket, object.getKey());
                final String cacheFilename = ICBFilenameUtils.concatPathParts(
                        s3LocalCacheDir, object.getKey());
                boolean fileUsedFromCache = false;
                LOG.debug("Checking in local s3 cache for " + cacheFilename);
                if (cacheFilename != null) {
                    final File cacheFile = new File(cacheFilename);
                    if (cacheFile.exists()) {
                        if (cacheFile.length() == s3FileLength) {
                            files.add(cacheFile);
                            fileUsedFromCache = true;
                            LOG.debug("... Retrieved from s3 cache");
                        } else {
                            LOG.info("... cache file size didn't match " +
                                    "cache=" + cacheFile.length() + " s3=" + s3FileLength);
                        }
                    } else {
                        LOG.debug("... cache file didn't exist.");
                    }
                }

                if (!fileUsedFromCache) {
                    final File localFilename = newLocalFile(dir, FilenameUtils.getName(object.getKey()));
                    dp.add(new DownloadPackage(object, localFilename));
                }
            }

            if (!dp.isEmpty()) {
                // Files to be fetched from S3 not in the local s3 cache
                final DownloadPackage[] dpArray = dp.toArray(new DownloadPackage[dp.size()]);
                final String transferTag = Tagged.createStringTag();
                LOG.info("Starting s3 transfer " + dp.size() + " files, transfer tag = "
                        + transferTag);
                multi.downloadObjects(s3Bucket, dpArray);
                LOG.info("Completed s3 transfer " + dp.size() + " files, transfer tag = "
                        + transferTag);

                for (final DownloadPackage d : dpArray) {
                    final File dataFile = d.getDataFile();
                    files.add(dataFile);
                }
            }

            // now, let's cleanup (not obvious from the Jets3t API, and not clear it has any effect).
            for (final S3Object object : content) {
                try {
                    object.closeDataInputStream();
                } catch (IOException e) {
                    LOG.error(e);
                }
            }

            return files.toArray(new File[files.size()]);
        } catch (S3ServiceException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    private File newLocalFile(final File dir, final String source) {
        return new File(ICBFilenameUtils.concatPathParts(dir.getPath(), FilenameUtils.getName(source)));
    }

    public long s3ObjectSize(final S3Bucket bucket, final String key) throws S3ServiceException {
        final RestS3Service service = initializeS3Connect();
        final S3Object s3Object = service.getObject(bucket, key);
        final long length = s3Object.getContentLength();
        LOG.debug(String.format("s3://%s/%s has size=%d", bucket.getName(), key, length));
        try {
            s3Object.closeDataInputStream();
        } catch (IOException e) {
            // ignore, we are just trying to close resources.
        }
        return length;
    }

    public long size() {
        if (size != -1) {
            return size;
        }
        try {
            size = s3ObjectSize(new S3Bucket(bucket), key);
        } catch (S3ServiceException e) {
            // cannot get the actual size for the resource. Perhaps it does not exist yet. Avoid
            // querying S3 again and again in that case.
            size = -1;
            throw new RuntimeException(e);
        }
        return size;
    }

    /**
     * Calls size but squishes RuntimeError - notable useful for files that don't yet exist.
     * @return the size or -1 if the file doesn't exist.
     */
    public long safeSize() {
        try {
            return size();
        } catch (RuntimeException e) {
            return -1;
        }
    }

    public long getTimeStamp() {
        if (timeStamp != -1) {
            return timeStamp;
        }
        final RestS3Service service = initializeS3Connect();

        try {
            final S3Object s3Object = service.getObject(new S3Bucket(bucket), key);
            timeStamp = s3Object.getLastModifiedDate().getTime();
            try {
                s3Object.closeDataInputStream();
            } catch (IOException e) {
                // ignore, we are just trying to close resources.
            }
            LOG.debug(String.format("s3://%s/%s has timeStamp=%d", bucket, key, timeStamp));
            return timeStamp;
        } catch (S3ServiceException e) {
            // cannot get the actual time stamp for the resource. Perhaps it does not exist yet. Avoid
            // querying S3 again and again in that case.
            timeStamp = 0;
            return timeStamp;
        }
    }

    /**
     * Write the files in the local copy to S3, under the specified path.
     *
     * @param localCopy
     */
    public void write(final LocalFileSetCopy localCopy) {
        final RestS3Service service = initializeS3Connect();
        try {
            final S3ServiceSimpleMulti multi = new S3ServiceSimpleMulti(service);

            final File[] files = localCopy.getFiles();
            final S3Object[] s3Objects = new S3Object[files.length];
            int i = 0;
            final S3Bucket s3Bucket = new S3Bucket(bucket);
            for (final File file : files) {

                final S3Object s3Object = new S3Object(s3Bucket, file);
                final String filename = ICBFilenameUtils.concatPathParts(key, file.getName());
                s3Object.setKey(filename);
                filenames.add(filename);
                s3Objects[i++] = s3Object;
            }
            multi.putObjects(s3Bucket, s3Objects);
            timeStamp = -1;
            size = -1;
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("Could not write to resource.", e);
        } catch (S3ServiceException e) {
            throw new RuntimeException("Could not write to resource.", e);
        } catch (IOException e) {
            throw new RuntimeException("Could not write to resource.", e);
        }
    }


    private RestS3Service initializeS3Connect() {
        final GroovyProperties config = ConfigHelper.loadConfiguration();

        final String awsSecretAccessKey = config.assertGet(GobyPropertyKeys.AWS_SECRET_ACCESS_KEY);
        final String awsAccessKey = config.assertGet(GobyPropertyKeys.AWS_ACCESS_KEY);

        final AWSCredentials credentials = new AWSCredentials(awsAccessKey, awsSecretAccessKey);

        // validate the credentials by listing buckets:
        try {
            final RestS3Service service = new RestS3Service(credentials);
            final S3Bucket[] myBuckets = service.listAllBuckets();

            /* System.out.println("How many buckets to I have in S3? " + myBuckets.length);
           for (S3Bucket b: myBuckets) {
               System.out.println("bucket: "+b.getName());
           } */
            return service;
        } catch (S3ServiceException e) {
            throw new RuntimeException("AWS credentials could not be used to list buckets", e);
        }
    }

    public static String createLocalTempDirectory() throws IOException {
        final String workDirectory = ConfigHelper.loadConfiguration().assertGet(GobyPropertyKeys.WORK_DIRECTORY);
        final int tagLength = 10;
        final String newTag = Tagged.createStringTag(tagLength);

        final String path = ICBFilenameUtils.concatPathParts(workDirectory, "local-tmp", newTag);
        final File localDir = new File(path);
        if (!localDir.exists()) {
            FileUtils.forceMkdir(localDir);
        }
        return path;
    }

    @Override
    public String toString() {
        return String.format("[ s3 replicate copy  url: %s timeStamp: %d size: %d ]",
                url, getTimeStamp(), safeSize());
    }

    @Override
    public int hashCode() {
        // you pick a hard-coded, randomly chosen, non-zero, odd number
        // ideally different for each class
        return new HashCodeBuilder(13, 17)
            .append(url)
            .append(fileType)
            .append(bucket)
            .append(key)
            .append(filenames)
            .append(resourceIsFileList)
            .append(getTimeStamp())
            .append(safeSize())
            .toHashCode();
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null) {
            return false;
        }
        if (obj == this) {
            return true;
        }
        if (obj.getClass() != getClass()) {
            return false;
        }
        final S3FileSetCopy rhs = (S3FileSetCopy) obj;
        return new LoggingEqualsBuilder(this.getClass(), Level.INFO)
            .append("url", url, rhs.url)
            .append("fileType", fileType, rhs.fileType)
            .append("bucket", bucket, rhs.bucket)
            .append("key", key, rhs.key)
            .append("filenames", filenames, rhs.filenames)
            .append("resourceIsFileList", resourceIsFileList, rhs.resourceIsFileList)
            .append("timeStamp", getTimeStamp(), rhs.getTimeStamp())
            .append("size", safeSize(), rhs.safeSize())
            .isEquals();
    }
}
