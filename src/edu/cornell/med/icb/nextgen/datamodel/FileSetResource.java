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

package edu.cornell.med.icb.nextgen.datamodel;

import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.Serializable;
import java.util.HashSet;

/**
 * @author Fabien Campagne
 *         Date: Jul 24, 2009
 *         Time: 5:26:18 PM
 */
public class FileSetResource implements Resource, Serializable {

    private static final long serialVersionUID = 8280944712943951564L;

    private String id;
    private String fileType;
    private static final Logger LOG = Logger.getLogger(FileSetResource.class);

    final int delay = 1000;
    final int timeout = 1000 * 60 * 2; // 2 minute timeout.

    HashSet<FileSetCopy> replicateCopies;

    public FileSetResource() {
        replicateCopies = new HashSet<FileSetCopy>();
    }

    @Override
    public int hashCode() {
        return id.hashCode();
    }

    @Override
    public boolean equals(final Object o) {
        if (!(o instanceof FileSetResource)) {
            return false;
        }

        final String o_id = ((FileSetResource) o).id;
        return id.equals(o_id);
    }

    @Override
    public synchronized String toString() {
        final StringBuilder buffer = new StringBuilder();
        buffer.append(String.format("[resource id:%s type:%s #replicates:%d", id, fileType, replicateCopies.size()));
        for (final FileSetCopy copy : replicateCopies) {
            buffer.append(' ');
            buffer.append(copy.toString());
            buffer.append("\t");
        }
        buffer.append(" ]");
        return buffer.toString();
    }


    public FileSetResource(final String id, final String url, final String fileType) {
        this.id = id;
        replicateCopies = new HashSet<FileSetCopy>();
        replicateCopies.add(defineByLocation(url, fileType));
        this.fileType = fileType;
    }

    /**
     * Define a resource that will be written to. A URL will be created when the file is written.
     *
     * @param id
     * @param fileType
     */
    public FileSetResource(final String id, final String fileType) {
        this.id = id;
        replicateCopies = new HashSet<FileSetCopy>();
        this.fileType = fileType;
    }


    private FileSetCopy defineByLocation(final String url, final String fileType) {
        if (S3FileSetCopy.isS3Url(url)) {
            return new S3FileSetCopy(url, fileType);
        } else if (url.startsWith("file:/")) {
            final String filenameWithoutPrefix = url.substring(6);
            return new LocalFileSetCopy(filenameWithoutPrefix, fileType);
        } else {
            return new LocalFileSetCopy(url, fileType);
        }
    }

    public void makeLocal() {
        LOG.info("Making local copy of resource: " + id);
        if (!isLocal()) {
            waitForReplicates(delay, timeout);
            synchronized (this) {
                for (final FileSetCopy copy : replicateCopies) {
                    if (!copy.isLocal()) {
                        // getFiles on a remote copy causes the download of the file content onto the local node:
                        final File[] files = copy.getFiles();
                        final LocalFileSetCopy localCopy = newLocalCopy(getParentDirectory(files));
                        for (final File file : files) {
                            localCopy.add(file, FilenameUtils.getExtension(file.getName()));
                        }
                        replicateCopies.add(localCopy);
                        return;
                    }
                }
            }
        }
        LOG.info("Done making local copy of resource: " + id);
    }

    private void waitForReplicates(final int delay, final int timeout) {
        int waitDuration = 0;
        if (replicateCopies.size() == 0) {
            do {
                try {
                    // the cache may need time to propagate the new state of this resource. The state will change
                    // through a call to synchronizeSate()
                    Thread.sleep(delay);
                    waitDuration += delay;
                    if (waitDuration > timeout) {
                        throw new RuntimeException(String.format("Waited for resource %s to update, but no replicate could be found to makeLocal().",
                                getId()));
                    } else {
                        LOG.warn(String.format("Waiting until resource %s has at least one replicate and can be made local",
                                getId()));
                    }
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            } while (replicateCopies.size() == 0);
        }
    }


    public synchronized void clearLocal() {
        for (final FileSetCopy replicate : replicateCopies) {
            if (replicate.isLocal()) {
                replicate.delete();
                replicateCopies.remove(replicate);
            }
        }
    }


    public synchronized boolean isLocal() {
        for (final FileSetCopy replicate : replicateCopies) {
            if (replicate.isLocal()) {
                return true;
            }
        }
        return false;
    }

    public synchronized File[] getLocalFiles(final String... fileType) {
        for (final FileSetCopy replicate : replicateCopies) {
            if (replicate.isLocal()) {
                return replicate.getFiles(fileType);
            }
        }
        return null;
    }

    /**
     * Size is unknown initially.
     */
    long size = -1;

    /**
     * Find the size of the resource without making it local.
     *
     * @return
     */
    public synchronized long peekAtSize() {
        for (final FileSetCopy replicate : replicateCopies) {
            size = replicate.size();
            if (size != 0) {
                return size;
            }
        }
        return -1;
    }

    public synchronized long getSize() {
        if (size != -1) {
            return size;
        }
        makeLocal();
        // get the size of any replicate in the set:
        for (final FileSetCopy replicate : replicateCopies) {
            if (replicate.isLocal()) {
                size = replicate.size();
                // we cache the size in the resource.
                return size;
            }
        }
        return 0;

    }

    public void write(final LocalFileSetCopy localCopy) {
        synchronized (this) {
            if (replicateCopies.size() == 0) {
                // Use the local copy as only copy.
                replicateCopies.add(localCopy);
                LOG.info("Added localCopy to Resource ID: " + id + " fileType: " + fileType + "  " + this);
                // force setting the value of size from the new replicate.
                forceEvaluateSize();
                return;
            }
            assert replicateCopies.size() == 1 : "Cannot write to a resource that has more than one replicate copy.";

            final FileSetCopy destination = replicateCopies.iterator().next();
            destination.write(localCopy);
            forceEvaluateSize();
        }
    }

    private void forceEvaluateSize() {
        // force setting the value of size from the new replicate.
        size = -1;
        getSize();
    }

    public String getId() {
        return id;
    }

    public synchronized long getTimeStamp() {
        long maxTimeStamp = Long.MIN_VALUE;

        for (final FileSetCopy replicate : replicateCopies) {
            maxTimeStamp = Math.max(replicate.getTimeStamp(), maxTimeStamp);
        }
        return maxTimeStamp;
    }


    public void synchronizeState(final Resource incomingResource) {
        if (!(incomingResource instanceof FileSetResource)) {
            return;
        }
        final FileSetResource remoteResource = (FileSetResource) incomingResource;
        /*
       Synchronize all accesses to replicateCopies on this to avoid the following exception:
   java.util.ConcurrentModificationException
       at java.util.HashMap$HashIterator.nextEntry(HashMap.java:793)
       at java.util.HashMap$KeyIterator.next(HashMap.java:828)
       at edu.cornell.med.icb.nextgen.datamodel.FileSetResource.getTimeStamp(FileSetResource.java:250)
       at edu.cornell.med.icb.nextgen.datamodel.Job.updateResource(Job.java:64)
       at edu.cornell.med.icb.gridgain.apps.AlignJob.synchronizeState(AlignJob.java:207)
        */
        synchronized (this) {
            final long remoteTime = remoteResource.getTimeStamp();
            final long localTime = getTimeStamp();
            LOG.debug(String.format(
                    "Doing a remote replicate synchronizeState%n   local=%s%n   remote=%s%n",
                    this, remoteResource));
            if (remoteTime > localTime) {
                // accept the new state of this local resource:
                // TODO could the replicateCopies set grow too large?
                replicateCopies.addAll(remoteResource.replicateCopies);
                LOG.info(String.format("Pulled replicates from remote resource %s, " +
                        " remote timeStamp is %d, local timeStamp was=%d, now is=%d",
                        remoteResource, remoteTime, localTime, getTimeStamp()));
                assert remoteTime == getTimeStamp();
            } else {
                LOG.debug(String.format(
                        "Not adding replicate because of time (remote) !< %d (local) %d",
                        remoteTime, localTime));
            }
        }
    }

    public LocalFileSetCopy newLocalCopy(final String location) {
        final LocalFileSetCopy localCopy = new LocalFileSetCopy();
        localCopy.location = location;
        synchronized (this) {
            if (replicateCopies == null) {
                replicateCopies = new HashSet<FileSetCopy>();
            }
            replicateCopies.add(localCopy);
        }
        return localCopy;
    }

    public static String getParentDirectory(final File[] alignFiles) {
        final ObjectSet<String> result = new ObjectOpenHashSet<String>();
        for (final File file : alignFiles) {
            assert file != null : "File cannot be null to determine parent directory.";
            result.add(file.getParent());
        }
        assert result.size() == 1 : "result directory must contain only one element for all the filenames in the alignment result. " + result.toString();
        return result.iterator().next();
    }
}
