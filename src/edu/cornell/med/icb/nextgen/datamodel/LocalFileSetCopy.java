package edu.cornell.med.icb.nextgen.datamodel;

import edu.cornell.med.icb.config.ConfigHelper;
import edu.cornell.med.icb.util.ExecuteProgram;
import edu.cornell.med.icb.util.GobyPropertyKeys;
import edu.cornell.med.icb.util.ICBFilenameUtils;
import edu.cornell.med.icb.util.LoggingEqualsBuilder;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Fabien Campagne
 *         Date: Jul 24, 2009
 *         Time: 5:40:44 PM
 */
public class LocalFileSetCopy implements FileSetCopy, Serializable {

    private static final long serialVersionUID = -5743342898912058018L;

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(LocalFileSetCopy.class);

    List<String> nodeIpAddresses;

    List<FileTypeAssociation> files = new ArrayList<FileTypeAssociation>();

    /**
     * location must be a cannonical path.
     */
    protected String location;

    private String url;

    boolean checkLocallyAccessible;

    /**
     * Empty constructor is required by serialization mechanism.
     */
    public LocalFileSetCopy() {
        nodeIpAddresses = new ArrayList<String>();
        nodeIpAddresses.add(getLocalNodeAddress());
    }

    @Override
    public String toString() {
        return String.format(
                "[ file replicate copy: num-files: %d, files: %s %s, node-ips: %s ]",
                files.size(), location, ArrayUtils.toString(files), nodeIpAddresses);
    }

    /**
     * Construct a local file set copy associated with the given hostname.
     *
     * @param localHostname
     */
    public LocalFileSetCopy(final String localHostname) {
        nodeIpAddresses = new ArrayList<String>();
        if (localHostname == null) {
            nodeIpAddresses.add(getLocalNodeAddress());
        } else {

            nodeIpAddresses.add(localHostname);

        }
    }

    public LocalFileSetCopy(final String url, final String fileType) {
        this(url, fileType, null);
    }

    public LocalFileSetCopy(String url, final String fileType, final String localHostname) {
        this(localHostname);
        if (url.startsWith("file:/")) {
            // strip out prefix from URL:
            url = url.substring(6);
        }
        this.url = url;
        final File location = new File(url);
        if (fileType.equals(Resource.FILE_LIST)) {
            try {
                this.location = location.getParentFile().getCanonicalPath();
                LOG.warn(String.format("file-list: location=%s", this.location));
            } catch (IOException e) {
                throw new RuntimeException(
                        "Cannot obtain cannonical path for file-list location " +
                        location.getPath(), e);
            }
            defineFileList(location);
        } else {
            // add each file in directory:
            if (!location.isDirectory()) {

                //     a single file:
                try {
                    final File parentFile = location.getParentFile();
                    if (parentFile != null) {
                        this.location = parentFile.getCanonicalPath();
                        LOG.warn(String.format("new-file: location=%s", this.location));

                    } else {
                        // no parent directory. The file is relative to the user directory.
                        final File currentDirectory = new File(System.getProperty("user.dir"));
                        this.location = currentDirectory.getCanonicalPath();
                        LOG.warn(String.format("new-file: location=%s", this.location));
                    }
                } catch (IOException e) {
                    throw new RuntimeException(
                            "Cannot obtain cannonical path of parent directory for location " +
                                    location.getPath(), e);
                }

            } else {
                try {
                    this.location = location.getCanonicalPath();
                    LOG.warn(String.format("new-file: location=%s", this.location));
                } catch (IOException e) {
                    throw new RuntimeException(
                            "Cannot obtain cannonical path for directory location " +
                                    location.getPath(), e);
                }
            }

            final FileTypeAssociation assoc = new FileTypeAssociation(location, fileType);
            files.add(assoc);
        }
    }

    /**
     * Retrieves the content of the file list at URL, and define one file association for each
     * filename defined in the list (the extension of the file is used as file type.
     */
    private void defineFileList(final File filePath) {
        FastBufferedReader reader = null;
        try {
            final String parentPath = cannonical(filePath.getParentFile());
            reader = new FastBufferedReader(new FileReader(filePath));
            final LineIterator lines = new LineIterator(reader);
            for (final MutableString line : lines.allLines()) {
                final String filename = line.trim().toString();
                LOG.trace(String.format("adding file %s to resource %s", filename, url));
                final String concatResult = ICBFilenameUtils.concatPathParts(parentPath, filename);
                assert concatResult != null : String.format(
                        "concat result cannot be null, parent path: %s filename: %s",
                        parentPath, filename);
                final File dataFile = new File(concatResult);

                files.add(new FileTypeAssociation(dataFile,
                        FilenameUtils.getExtension(filename)));
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Cannot define resource form file list at " + url, e);
        } finally {
            IOUtils.closeQuietly(reader);
        }
    }

    private String cannonical(File file) {
        try {
            return file.getCanonicalPath();
        } catch (IOException e) {
            LOG.error(e);
            return null;
        }
    }

    protected void setCheckLocallyAccessible(final boolean checkLocallyAccessible) {
        this.checkLocallyAccessible = checkLocallyAccessible;
    }

    public boolean isLocal() {
        final String localNodeAddress = getLocalNodeAddress();
        if (nodeIpAddresses.contains(localNodeAddress)) {
            return true;
        } else {
            if (checkLocallyAccessible) {
                // check if resource is locally accessible (i.e., accessible through some
                // filesystem shared with the other node.
                boolean accessible = true;
                for (final FileTypeAssociation assoc : files) {
                    if (!assoc.getFile(this).canRead()) {
                        // if a single file from files cannot be read, then the resource is
                        // not accessible.
                        accessible = false;
                    }
                }
                if (accessible) {
                    // add node uuid to list of nodes that can see this resource:
                    nodeIpAddresses.add(localNodeAddress);
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
    }

    private String getLocalNodeAddress() {
        // Determine local host dynamically at the time of execution.
        final String hostName = ConfigHelper.loadConfiguration().assertGet(
                GobyPropertyKeys.HOST_NAME);
        return hostName;

    }

    public void delete() {
        for (final FileTypeAssociation assoc : files) {
            FileUtils.deleteQuietly(assoc.getFile(this));
        }
    }

    public synchronized File[] getFiles(final String... fileTypes) {
        final ObjectArrayList<File> result = new ObjectArrayList<File>();
        if (isLocal()) {
            final ObjectArrayList<String> fileTypeList = new ObjectArrayList<String>(fileTypes);


            for (final FileTypeAssociation assoc : files) {
                if (fileTypeList.size() == 0 || fileTypeList.contains(assoc.fileType)) {
                    result.add(assoc.getFile(this));
                }
            }
        } else {
            // all the files belonging to this set must be copied to the same directory.
            // This is important when they all belong to a 'basename' set of files.
            final String localDestDir;
            try {
                localDestDir = S3FileSetCopy.createLocalTempDirectory();
            } catch (IOException e) {
                LOG.error("Directory could not be created", e);
                throw new RuntimeException("Directory could not be created", e);
            }

            for (final FileTypeAssociation assoc : files) {
                final String sourceFilePath = ICBFilenameUtils.concatPathParts(
                        location, assoc.filename);
                final String sourceFileNodeIp = getRandomHostname(nodeIpAddresses);


                final String localDestinationFilePath = ICBFilenameUtils.concatPathParts(
                        localDestDir, assoc.filename);

                final String privateKeyPath = System.getenv("PRIVATE_KEY_PATH");
                final String scpRemoteToLocalTemplate = ConfigHelper.loadConfiguration().
                        assertGetWithParameter(
                                GobyPropertyKeys.SCP_COMMAND_TEMPLATE_REMOTE_TO_LOCAL,
                                privateKeyPath);

                final String scpCommandLine = String.format(scpRemoteToLocalTemplate,
                        sourceFileNodeIp, sourceFilePath,
                        localDestinationFilePath);
                final ExecuteProgram scpProcess = new ExecuteProgram();
                try {
                    scpProcess.executeToLog(
                            scpCommandLine, LocalFileSetCopy.class, Level.INFO, "BWA[__OUTPUT_TAG__]: ");
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                } catch (IOException e) {
                    throw new RuntimeException(
                            "Received IO Exception while copying file from remote host", e);
                }
                result.add(new File(localDestinationFilePath));
            }
        }
        return result.toArray(new File[result.size()]);
    }

    private String getRandomHostname(final List<String> nodeAddresses) {
        Collections.shuffle(nodeAddresses);
        return nodeAddresses.get(0);
    }

    public long size() {
        long size = 0;
        for (final FileTypeAssociation assoc : files) {
            size += assoc.getFile(this).length();
        }
        return size;
    }

    /**
     * The most recent last modification file of any file in the set.
     *
     * @return the fileset timestamp.
     */
    public long getTimeStamp() {
        long maxTimeStamp = Long.MIN_VALUE;
        for (final FileTypeAssociation assoc : files) {
            maxTimeStamp = Math.max(assoc.getFile(this).lastModified(), maxTimeStamp);
        }
        return maxTimeStamp;
    }

    public void write(final LocalFileSetCopy localCopy) {
        if (localCopy.isLocal()) {
            for (final File src : localCopy.getFiles()) {
                final File dest = new File(ICBFilenameUtils.concatPathParts(
                        location, src.getName()));

                try {
                    if (!src.equals(dest)) FileUtils.copyFile(src, dest);
                } catch (IOException e) {
                    throw new RuntimeException(String.format(
                            "Error copying file from src %s to dest %s",
                            src.getPath(), dest.getPath()), e);
                }
            }
            files.addAll(localCopy.files);
        } else {
            throw new InternalError(
                    "Copying to remote node has not been implemented. Resources can be " +
                            "only be pulled from remote hosts at this time.");
        }
    }

    public void add(final File dataFile, final String fileType) {
        files.add(new FileTypeAssociation(dataFile, fileType));
    }

    @Override
    public int hashCode() {
        // you pick a hard-coded, randomly chosen, non-zero, odd number
        // ideally different for each class
        return new HashCodeBuilder(17, 19)
            .append(location)
            .append(url)
            .append(nodeIpAddresses)
            .append(files)
            .append(checkLocallyAccessible)
            .append(getTimeStamp())
            .append(size())
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
        final LocalFileSetCopy rhs = (LocalFileSetCopy) obj;
        return new LoggingEqualsBuilder(this.getClass(), Level.INFO)
            .append("location", location, rhs.location)
            .append("url", url, rhs.url)
            .append("nodeIpAddresses", nodeIpAddresses, rhs.nodeIpAddresses)
            .append("files", files, rhs.files)
            .append("checkLocallyAccessible", checkLocallyAccessible, rhs.checkLocallyAccessible)
            .append("timeStamp", getTimeStamp(), rhs.getTimeStamp())
            .append("size", size(), rhs.size())
            .isEquals();
    }
}
