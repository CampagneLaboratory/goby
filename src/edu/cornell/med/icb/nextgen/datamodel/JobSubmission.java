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

import edu.cornell.med.icb.gridgain.GobyCacheManager;
import it.unimi.dsi.fastutil.ints.Int2BooleanMap;
import it.unimi.dsi.fastutil.ints.Int2BooleanOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import javax.management.InstanceAlreadyExistsException;
import javax.management.MBeanRegistrationException;
import javax.management.MBeanServer;
import javax.management.MalformedObjectNameException;
import javax.management.NotCompliantMBeanException;
import javax.management.ObjectName;
import java.io.Serializable;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

/**
 * Describes a set of jobs and their execution status on the cluster.
 *
 * @author Fabien Campagne
 *         Date: May 22, 2009
 *         Time: 1:42:28 PM
 */

public class JobSubmission implements Serializable, JobSubmissionMBean {

    private static final long serialVersionUID = 5793414087779431012L;

    private static final Log LOG = LogFactory.getLog(JobSubmission.class);

    /**
     * Completed job outputs are stored in this map, keyed by int job tag:
     */
    private Int2ObjectMap outputs;
    /**
     * Indicates whether cyclic dependency was detected in this submission.
     */
    public boolean cyclicDependency = false;
    /**
     * The number of jobs in this submission.
     */
    private int totalNumJobs = 0;
    /**
     * The unique identifier of this submission.
     */
    protected String id;

    /** Job to execute upon JobSubmission failure. */
    public Job uponFailureJob;

    /** Job to execute upon JobSubmission success. */
    public Job uponSuccessJob;

    /** Job to execute upon JobSubmission completion. */
    public Job uponCompletionJob;

    /**
     * All of the resources for this JobSubmission.
     */
    protected ObjectList<Resource> resources = new ObjectArrayList<Resource>();

    /**
     * The helper for writing to the SQS output queue.
     */
    protected transient JobSubmissionQueueHelper outputQueueHelper = null;

    /*
     * Jobs that have been submitted for execution as part of this schedule. Jobs are ordered
     * by start time, on the machine where they should be executed.
     */
    private ArrayList<Job> submittedJobs;

    /*
     * Jobs that have been submitted and are ready to start. They either have no predecessor,
     * or their predecessors have successfully completed.
     */
    private ArrayList<Job> readyJobs;

    /**
     * Jobs that have started executing.
     */
    private ArrayList<Job> startedJobs;

    /**
     * Jobs that have completed execution without error or exception.
     */
    private ArrayList<Job> completedJobs;

    /**
     * Jobs that have failed.
     */
    private ArrayList<Job> failedJobs;

    /**
     * Define a resource accessible this job submission.
     *
     * @param resourceId And id for the resource.
     * @param url        URL of the resource.
     * @param fileType   Type of the resource.
     * @return a proxy to the resource.
     */
    public Resource defineResource(final String resourceId, final String url, String fileType) {
        final FileSetResource r = new FileSetResource(resourceId, url, fileType);
        resources.add(r);
        return r;
    }

    /**
     * Define a resource that this job submission can write to (and subsequently read from).
     *
     * @param resourceId And id for the resource.
     * @param fileType   Type of the resource.
     * @return a proxy to the resource.
     */
    public Resource defineResource(final String resourceId, final String fileType) {
        final FileSetResource r = new FileSetResource(resourceId, fileType);
        resources.add(r);
        return r;
    }

    /**
     * Create a new job submission.
     */
    public JobSubmission() {
        submittedJobs = new ArrayList<Job>();
        readyJobs = new ArrayList<Job>();
        startedJobs = new ArrayList<Job>();
        completedJobs = new ArrayList<Job>();
        failedJobs = new ArrayList<Job>();
        outputs = new Int2ObjectOpenHashMap();
        id = "submission-" + Tagged.createStringTag();
        final MBeanServer server = ManagementFactory.getPlatformMBeanServer();

        try {
            final ObjectName objectName = new ObjectName(
                    String.format("goby:type=JobSubmission,name=%s", id));

            if (!server.isRegistered(objectName)) {
                server.registerMBean(this, objectName);
            }
        } catch (InstanceAlreadyExistsException e) {
            LOG.error(e);
        } catch (MBeanRegistrationException e) {
            LOG.error(e);
        } catch (NotCompliantMBeanException e) {
            LOG.error(e);
        } catch (MalformedObjectNameException e) {
            LOG.error(e);
        }

    }

    /**
     * Returns true when this schedule contains a job that is ready to start,
     * false when all jobs are either completed or waiting on a
     * predecessor to finish.
     *
     * @return True or False
     */
    public boolean hasJobReadyToStart() {
        return readyJobs.size() > 0;

    }

    /**
     * Returns the next job ready to start. This method returns null if hasJobReadyToStart return false.
     *
     * @return gets the next "ready" job.
     */
    public Job nextReadyJob() {
        determineReadyJobs();
        if (hasJobReadyToStart()) {
            return readyJobs.get(0);
        } else {
            // inspect
            return null;
        }
    }

    /**
     * Inspect submitted jobs and promote to ready status when appropriate.
     */
    public synchronized void determineReadyJobs() {
        assertJobCount();
        final List<Job> toUpdate = copy(submittedJobs);
        for (Job j : toUpdate) {
            if (j.isReady()) {
                // promote to ready status:
                j = GobyCacheManager.getInstance().getJobFromCache(j.getTag());
                removeFromAnyList(j);
                readyJobs.add(j);

            }
        }

        assertJobCount();
    }

    /**
     * Reorganize jobs according to status.
     */
    public synchronized void reconsiderAllJobs() {
        assertJobCount();
        final List<Job> toUpdate = copy(allJobs());
        for (Job j : toUpdate) {
            j = GobyCacheManager.getInstance().getJobFromCache(j.getTag());
            update(j);
        }

        assertJobCount();
    }

    /**
     * Retrieve the status of finally jobs from the cache if needed.
     */
    public void updateFinallyJobStatus() {
        final GobyCacheManager cache = cache();
        if (uponCompletionJob != null) {
            uponCompletionJob = cache.getJobFromCache(uponCompletionJob.getTag());
        }
        if (uponFailureJob != null) {
            uponFailureJob = cache.getJobFromCache(uponFailureJob.getTag());
        }
        if (uponSuccessJob != null) {
            uponSuccessJob = cache.getJobFromCache(uponSuccessJob.getTag());
        }
    }

    private List<Job> allJobs() {
        ArrayList<Job> allJobs = new ArrayList<Job>();
        allJobs.addAll(submittedJobs);
        allJobs.addAll(readyJobs);
        allJobs.addAll(startedJobs);
        allJobs.addAll(completedJobs);
        allJobs.addAll(failedJobs);
        return allJobs;
    }

    private List<Job> copy(List<Job> toBeCopied) {
        final List<Job> result = new ArrayList<Job>();
        result.addAll(toBeCopied);
        return result;
    }

    /**
     * Start a job on the cluster.
     *
     * @param job the job to start
     */
    public synchronized void startJob(final Job job) {
        assertJobCount();
        removeFromAnyList(job);
        startedJobs.add(job);
        job.executing();
        assertJobCount();
    }

    /**
     * Update this schedule with the new status of the job.
     *
     * @param job the job to update
     */
    public void update(final Job job) {
        synchronized (this) {
            assertJobCount();
            switch (job.status) {
                case STARTED: {
                    removeFromAnyList(job);
                    startedJobs.add(job);
                    break;
                }
                case READY: {
                    removeFromAnyList(job);
                    readyJobs.add(job);
                    break;
                }
                case COMPLETED: {
                    removeFromAnyList(job);
                    completedJobs.add(job);
                    break;
                }
                case FAILED: {
                    removeFromAnyList(job);
                    failedJobs.add(job);
                    someFailure(job);
                    break;
                }
                default:
            }
            assertJobCount();
        }
        determineReadyJobs();
        cache().pushJobToCache(job);
    }

    private void assertJobCount() {
        int currentNumJobs = (submittedJobs.size() + readyJobs.size() +
                startedJobs.size() + completedJobs.size() +
                failedJobs.size());
        if (currentNumJobs != this.totalNumJobs) {
            throw new InternalError("The number of jobs in various state lists " +
                    "must match total number of jobs in submission.");
        }
    }

    private void removeFromAnyList(Job taggedJob) {
        submittedJobs.remove(taggedJob);
        startedJobs.remove(taggedJob);
        readyJobs.remove(taggedJob);
        completedJobs.remove(taggedJob);
        failedJobs.remove(taggedJob);
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append("JobSubmission id=").append(id).append("\n");
        if (cyclicDependency) {
            sb.append("Detected cyclic dependency in this job submission..");
        }
        printJobList(sb, submittedJobs);
        printJobList(sb, readyJobs);
        printJobList(sb, startedJobs);
        printJobList(sb, completedJobs);
        printJobList(sb, failedJobs);
        return sb.toString();
    }

    private void printJobList(final StringBuffer sb, final List<Job> list) {
        for (final Job j : list) {
            sb.append(j);
            sb.append('\n');
        }
    }

    /**
     * Return a job in the collection jobList that matches the argument job tag.
     *
     * @param job     Job which provides the tag to look for
     * @param jobList Collection where the job will be located based on tag
     * @return The job from jobList that matches the specified tag, or null.
     */
    public static Job getJobByTag(final Job job, final Collection<Job> jobList) {
        for (final Job j : jobList) {
            if (j != null && j.getIntTag() == job.getIntTag()) {
                return j;
            }
        }
        return null;
    }

    /**
     * Return true when they are no more jobs in submitted stage.
     *
     * @return True or False.
     */
    public boolean doneSubmitting() {
        final boolean result = (submittedJobs.size() == 0);
        determineReadyJobs();
        return result;
    }

    /**
     * Return true when all the jobs either completed successfully or failed, or early termination was triggered.
     *
     * @return True or False.
     */
    public boolean done() {
        determineReadyJobs();
        final boolean result = (completedJobs.size() + failedJobs.size()) == totalNumJobs;
        return result || earlyTermination();
    }


    /**
     * Return true when they are some jobs currently executing on the cluster.
     *
     * @return True or False.
     */
    public boolean hasJobExecuting() {
        return this.startedJobs.size() > 0;
    }

    /**
     * Returns the output of a successful job previously executed:
     *
     * @param tag Int tag of the job for which output is sought.
     * @return Null if the output is not ready.
     */
    public Object getOutput(int tag) {
        return outputs.get(tag);
    }

    /**
     * Add a job to this submission.
     *
     * @param jobs Jobs to execute as part of this submission.
     */
    public void add(final Job... jobs) {
        for (final Job job : jobs) {
            submittedJobs.add(job);
            totalNumJobs++;
            job.jobSubmissionId = this.getId();
        }
    }

    /**
     * Returns true if there are no failed jobs and all jobs have completed.
     *
     * @return True or False.
     */
    public boolean success() {
        return failedJobs.size() == 0 && completedJobs.size() == totalNumJobs;
    }

    /**
     * Validate the submission. Check that the submission does not contain cyclic dependency.
     *
     * @return true if this submission validated
     */
    public boolean validate() {

        final ObjectArrayList<Job> orderedList = new ObjectArrayList<Job>();
        final Int2BooleanMap visitMap = new Int2BooleanOpenHashMap();
        for (Job job : submittedJobs) {
            visitMap.clear();
            visitJob(job, orderedList, visitMap);
        }
        submittedJobs.clear();
        submittedJobs.addAll(orderedList);
        return !cyclicDependency;
    }

    private void visitJob(Job job, ObjectList<Job> orderedList, Int2BooleanMap visitMap) {
        final int tag = job.getIntTag();
        if (!visitMap.get(tag)) {
            visitMap.put(tag, true);
            for (final Job p : job.getPredecessors()) {
                if (!orderedList.contains(p)) {
                    visitJob(p, orderedList, visitMap);
                }
            }
            if (!orderedList.contains(job)) orderedList.add(job);
        } else {
            cyclicDependency = true;
            throw new UnsupportedOperationException(
                    "The submission graph cannot contain cyclic dependencies.");
        }
    }

    /**
     * Return job submission id.
     * @return the job submission id.
     */
    public String getId() {
        return id;
    }

    /**
     * Return a copy of the list of jobs ready to a start.
     *
     * @return the next ready jobs list
     */
    public List<Job> nextReadyJobs() {
        final List<Job> result = new ArrayList<Job>();
        for (final Job ready : readyJobs) {
            result.add(ready);
        }
        return result;
    }

    /**
     * Configure this submission to run job upon failure of the submission.
     *
     * @param job This job will be run immediately after the submission has failed.
     */
    public void uponFailure(final Job job) {
        this.uponFailureJob = job;
    }

    /**
     * Configure this submission to run job upon success of the submission.
     *
     * @param job This job will be run immediately after the submission has succeeded.
     */
    public void uponSuccess(final Job job) {
        this.uponSuccessJob = job;
    }

    /**
     * Configure this submission to run job upon completion of the submission. A job configured to run upon completion
     * always runs after uponFailure or uponSuccess jobs;
     *
     * @param job This job will be run as the very last part of the job submission, irrespective of submission execution status.
     */
    public void uponCompletion(final Job job) {
        this.uponCompletionJob = job;
    }

    /**
     * Return true if this job submission has been finalized, false otherwise. A job submission
     * is finalized when all the uponFailure/Success/Completion jobs have been run
     * (either failed or succeeded).
     *
     * @return true if finalized
     */
    public boolean finalized() {
        final Boolean bool = (Boolean) cache().get(
                "/gridgain/submissions/" + getId() + "/finalized");
        return bool != null && bool;
    }

    /**
     * Set finalized for this job.
     * @param finalized if it is finalized
     */
    public void setFinalized(boolean finalized) {
        cache().put("/gridgain/submissions/" + getId() + "/finalized", finalized);
    }

    /**
     * Let this submission know that some job failed. The submission can decide to abort early,
     * or can try and execute the job again.
     *
     * @param job the job that failed
     * @return True when the job should be re-executed, or false otherwise.
     */
    public boolean someFailure(final Job job) {
        cache().put("/gridgain/submissions/" + getId() + "/early-termination", true);
        // earlyTermination = true;
        return false;
    }

    /**
     * True when the job submission termination policy calls for early termination and a
     * job has failed.
     * @return true if failure calls for early termination
     */
    public boolean earlyTermination() {
        final Boolean bool = (Boolean) cache().get(
                "/gridgain/submissions/" + getId() + "/early-termination");
        return bool != null && bool;
    }

    /**
     * Clear all cached outputs for this submission.
     */
    public void clearOutputs() {
        ObjectSet<Job> allJobs = new ObjectArraySet<Job>();
        allJobs.addAll(submittedJobs);
        allJobs.addAll(startedJobs);
        allJobs.addAll(completedJobs);
        allJobs.addAll(failedJobs);

        /** TODO: WE NEED TO CLEAR THIS SUBMISSION. */
        cache().clearJobOutput(allJobs.toArray(new Job[allJobs.size()]));
    }

    @Override
    /**
     * Clear the cache from this submission job outputs upon finalize.
     */
    protected void finalize() throws Throwable {
        super.finalize();
        clearOutputs();
        final MBeanServer server = ManagementFactory.getPlatformMBeanServer();
        server.unregisterMBean(new ObjectName(id));
    }

    public String getStatus() {
        final int percentSubmitted = divide(getSubmittedJobs(), totalNumJobs);
        final int percentReady = divide(getReadyJobs(), totalNumJobs);
        final int percentStarted = divide(getExecutingJobs(), totalNumJobs);
        final int percentCompleted = divide(getCompletedJobs(), totalNumJobs);
        final int percentFailed = divide(getFailedJobs(), totalNumJobs);
        return String.format(
                " submitted: %d%% ready: %d%% started: %d %% completed %d %% failed %d %% ",
                percentSubmitted,
                percentReady,
                percentStarted, percentCompleted, percentFailed);
    }

    private int divide(int a, int b) {
        return (int) (((double) a * 100) / ((double) b));
    }

    public List<Job> getSubmittedJobsList() {
        return submittedJobs;
    }

    public int getSubmittedJobs() {
        return submittedJobs.size();
    }

    public int getExecutingJobs() {
        return startedJobs.size();
    }

    public int getReadyJobs() {
        return readyJobs.size();
    }


    public int getCompletedJobs() {
        return completedJobs.size();
    }

    public int getFailedJobs() {
        return failedJobs.size();
    }

    public ObjectList<Resource> getResources() {
        return resources;
    }

    public synchronized String[] getReadyToStartIds() {
        reconsiderAllJobs();
        determineReadyJobs();
        final List<String> tags = new LinkedList<String>();
        GobyCacheManager cache = GobyCacheManager.getInstance();
        for (Job job : readyJobs) {
            job.executing();
            cache.pushJobToCache(job);
            tags.add(job.getTag());
        }
        return tags.toArray(new String[tags.size()]);
    }

    public Job getUponSuccessJob() {
        updateFinallyJobStatus();
        return uponSuccessJob;
    }

    public Job getUponFailureJob() {
        updateFinallyJobStatus();
        return uponFailureJob;
    }

    public Job getUponCompletionJob() {
        updateFinallyJobStatus();
        return uponCompletionJob;
    }

    public void setOutputQueueHelper(final JobSubmissionQueueHelper outputQueueHelper) {
        this.outputQueueHelper = outputQueueHelper;
    }

    public JobSubmissionQueueHelper getOutputQueueHelper() {
        return this.outputQueueHelper;
    }

    /**
     * Shortcut method for getting the cache.
     * @return the cache
     */
    public GobyCacheManager cache() {
        return GobyCacheManager.getInstance();
    }
}