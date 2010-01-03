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

import com.xerox.amazonws.sqs2.SQSException;
import edu.cornell.med.icb.goby.cache.GobyCacheManager;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.jets3t.service.S3ServiceException;

import java.io.IOException;
import java.io.Serializable;
import java.security.NoSuchAlgorithmException;
import java.util.HashSet;

/**
 * A job which transforms inputs into outputs. Jobs can be executed by calling their execute method, or run on a grid
 * by submitting the job to the grid.
 *
 * @author Fabien Campagne
 *         Date: May 22, 2009
 *         Time: 12:47:02 PM
 * @see JobSubmission
 */
public class Job extends Tagged implements Serializable {

    private static final long serialVersionUID = -8920461190471775456L;

    private static final Log LOG = LogFactory.getLog(Job.class);

    /** The id of the parent JobSubmission. */
    public String jobSubmissionId;

    /**
     * The jobs that must finish before this job can start.
     */
    private final HashSet<Job> predecessors;

    /**
     * The input of this job, if any.
     */
    protected Serializable input;

    /**
     * The status of this job. Status indicates the progress of the job in the execution pipeline.
     */
    public JobStatus status = JobStatus.NEW;

    /**
     * The throwable that failed the job, if any.
     */
    private Throwable throwable;

    /**
     * A description of the reason why this job failed, if any.
     */
    private String reasonFailed;

    private void setStatus(final JobStatus newStatus) {
        final boolean statusChanged = status != newStatus;
        status = newStatus;

        // Log the new status for this job
        String failure = "";
        if (reasonFailed != null) {
            failure += " Reason failed: " + reasonFailed;
        }
        final String statusMessage;
        if (throwable != null) {
            failure += " Failure class=" + throwable.getClass().getName()  + " message=" + throwable.getMessage();
            statusMessage = String.format("Job status change: status=%s, class=%s, tag=%s %s",
                    status.toString(), getClass().getSimpleName(),
                    getTag(), failure);
            LOG.info(statusMessage, throwable);
        } else {
            statusMessage = String.format("Job status change: status=%s, class=%s, tag=%s %s",
                    status.toString(), getClass().getSimpleName(),
                    getTag(), failure);
            LOG.info(statusMessage);
        }

        if (statusChanged) {
            final GobyCacheManager cache = GobyCacheManager.getInstance();
            final JobSubmission baseSubmission = cache.getJobSubmission(jobSubmissionId);
            final JobSubmission statusSubmission = StatusJob.createStatusJobSubmission(
                    baseSubmission, statusMessage);
            if (statusSubmission != null) {
                try {
                    statusSubmission.getOutputQueueHelper().addJobSubmissionToSQSQueue(
                        statusSubmission);
                } catch (IOException e) {
                    LOG.error(e);
                } catch (S3ServiceException e) {
                    LOG.error(e);
                } catch (SQSException e) {
                    LOG.error(e);
                } catch (NoSuchAlgorithmException e) {
                    LOG.error(e);
                }
            }
        }
    }

    /**
     * Override this method to implement job specific processing.
     *
     * @return The job itself. Job output should be stored in the cache.
     * @throws InterruptedException interrupted during execution
     */
    public Serializable execute() throws InterruptedException {
        System.out.println("Executed job " + getTag());
        return this;
    }

    /**
     * Update the state of a resource to/from the distributed cache.
     *
     * @param resource       The resource to update.
     * @return the updated resource
     */
    protected Resource updateResource(final Resource resource) {

        if (resource == null) {
            return null;
        }
        final GobyCacheManager cache = GobyCacheManager.getInstance();

        final Resource cacheResourceCopy = cache.getResourceFromCache(jobSubmissionId, resource);

        if (cacheResourceCopy == null) {
            LOG.info("Pushing my new resource " + resource);
            cache.pushResourceToCache(jobSubmissionId, resource);
        } else {
            final long remoteTime = cacheResourceCopy.getTimeStamp();
            final long localTime = resource.getTimeStamp();
            if (cacheResourceCopy.getTimeStamp() > resource.getTimeStamp()) {
                // the cache time stamp is more recent.
                resource.synchronizeState(cacheResourceCopy);
            } else if (cacheResourceCopy.getTimeStamp() < resource.getTimeStamp()) {
                LOG.info(String.format("Pushing my more current resource %s / " +
                        "remote timeStamp=%d, local timeStamp=%d",
                        resource, remoteTime, localTime));
                cache.pushResourceToCache(jobSubmissionId, resource);
            }
            // As write can only occur with zero replicates we shouldn't have two timestamps?
        }
        return resource;
    }

    /**
     * Synchronize the state of this job with copies of the job on the cluster.
     *
     * @return Return the job that has the most current state, or null, if 'this' is more current.
     */
    public Job synchronizeState() {
        return synchronizeState(GobyCacheManager.getInstance().getJobFromCache(getTag()));
    }

    /**
     * Synchronize the state of this job with copies of the job on the cluster.
     *
     * @param remoteJob the remote job to synchronize with
     * @return Return the job that has the most current state, or null, if 'this' is more current.
     */
    public Job synchronizeState(final Job remoteJob) {
        if (remoteJob == null) {
            LOG.warn("Local job is: " + getTag() + ", Remote job was null");
            return null;
        }

        if (remoteJob == this) {
            LOG.debug("Remote job is the same as the local job (" + getTag() + ")");
            return null;
        }

        if (otherJobMoreAdvanced(remoteJob)) {
            LOG.debug("Updating status of job " + getTag() + " to " + remoteJob.status);
            setStatus(remoteJob.status);
            return remoteJob;
        }
        return this;
    }

    /**
     * Check if a job is further along than this job.
     * @param other the other job to compare against
     * @return true if other.status > this.status
     */
    private boolean otherJobMoreAdvanced(final Job other) {
        return other.status.getEncodedValue() > this.status.getEncodedValue();
    }

    /**
     * Construct a new job. The job is assigned a random tag.
     */
    public Job() {
        newTag();
        predecessors = new HashSet<Job>();

    }

    @Override
    public int hashCode() {
        return getIntTag();
    }

    @Override
    public boolean equals(final Object other) {
        if (!(other instanceof Job)) {
            return false;
        }
        final Job otherJob = (Job) other;
        return otherJob.getIntTag() == getIntTag();
    }


    /**
     * The job has been submitted to the execution pipeline.
     */
    public void submit() {
        setStatus(JobStatus.SUBMITTED);
    }

    /**
     * The job has started execution.
     */
    public void executing() {
        setStatus(JobStatus.STARTED);
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append("[").append(getTag()).append('/');
        sb.append(getClass().getSimpleName());
        sb.append(" status=").append(status);
        sb.append(" depends on ").append(listPredecessorTags());
        if (reasonFailed != null) {
            sb.append(" reasonFailed=").append(reasonFailed);
        }
        if (throwable != null) {
            sb.append(" throwable.message=").append(throwable.getMessage());
        }
        sb.append("]");
        return sb.toString();
    }

    /**
     * Returns a String of the comma separated list of predecessor tags.
     *
     * @return the predecessor tags
     */
    private String listPredecessorTags() {
        final StringBuilder sb = new StringBuilder();
        int i = 0;
        for (final Job p : predecessors) {
            if (i++ > 0) {
                sb.append(",");
            }
            sb.append(p.getTag());
        }
        return sb.toString();
    }

    /**
     * The job hgas completed successfully.
     */
    public void completed() {
        setStatus(JobStatus.COMPLETED);
    }

    /**
     * This job has failed with this exception.
     *
     * @param error The exception that caused the failure.
     */
    public void failed(final Throwable error) {
        throwable = error;
        setStatus(JobStatus.FAILED);
    }

    /**
     * This job has failed for the specified reason.
     *
     * @param reasonFailed A message which describes the reason the job failed.
     */
    public void failed(final String reasonFailed) {
        this.reasonFailed = reasonFailed;
        setStatus(JobStatus.FAILED);
    }

    /**
     * Add predecessors to this job.
     *
     * @param predecessors List of jobs that must be completed before this job can start.
     */
    public void addPredecessors(final Job... predecessors) {
        for (final Job p : predecessors) {
            if (JobSubmission.getJobByTag(p, this.predecessors) == null) {
                this.predecessors.add(p);
            }
        }
    }

    /**
     * Return true if the job is ready to start. A job is ready to start when it has no
     * predecessors, or when all its predecessors have successfuly completed.
     * @return true if this job is ready
     */
    public boolean isReady() {
        int completedCount = 0;
        final GobyCacheManager cache = GobyCacheManager.getInstance();
        for (Job p : predecessors) {
            p = cache.getJobFromCache(p.getTag());
            if (p.status == JobStatus.COMPLETED) {
                completedCount++;
            }
        }
        final boolean jobIsReady = (completedCount == predecessors.size());
        if (jobIsReady) {
            setStatus(JobStatus.READY);
        }
        return jobIsReady;
    }

    public HashSet<Job> getPredecessors() {
        return predecessors;
    }

    public Job setInput(final Serializable input) {
        this.input = input;
        return this;
    }

    /**
     * Returns a description of the reason why this job failed, if any.
     *
     * @return the failure reason string
       */
    public String getFailureReason() {
        return reasonFailed;
    }

    /**
     * Shortcut method for getting the cache.
     * @return the cache
     */
    public GobyCacheManager cache() {
        return GobyCacheManager.getInstance();
    }

}
