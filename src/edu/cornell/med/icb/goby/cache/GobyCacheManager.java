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

package edu.cornell.med.icb.goby.cache;

import edu.cornell.med.icb.nextgen.datamodel.Job;
import edu.cornell.med.icb.nextgen.datamodel.JobSubmission;
import edu.cornell.med.icb.nextgen.datamodel.JobWithCounter;
import edu.cornell.med.icb.nextgen.datamodel.Resource;
import edu.cornell.med.icb.util.ICBFilenameUtils;
import net.sf.ehcache.Cache;
import net.sf.ehcache.CacheManager;
import net.sf.ehcache.Element;
import net.sf.ehcache.Status;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * This data manager handles getting job data from and putting job data into
 * the Distributed EH Cache.
 */
public final class GobyCacheManager {
    private static final Log LOG = LogFactory.getLog(GobyCacheManager.class);

    /**
     * Singleton instance.
     */
    private static final GobyCacheManager INSTANCE = new GobyCacheManager();

    private Cache gobyCache;

    /**
     * Gets singleton.
     *
     * @return Singleton.
     */
    public static GobyCacheManager getInstance() {
        return INSTANCE;
    }

    /**
     * Ensure singleton.
     */
    private GobyCacheManager() {
        super();
    }

    /**
     * Gets the status of the job associated with the int tag.
     *
     * @param jobSubmissionId The ID of the job submission that contains the job referenced by intTag.
     * @return Cached job status (possibly <tt>null</tt>).
     */
    public JobSubmission getJobSubmission(final String jobSubmissionId) {
        return (JobSubmission) safeGet("/gridgain/job-submissions/" + jobSubmissionId);

    }

    private Serializable safeGet(final String key) {
        final Element e = gobyCache.get(key);
        if (e == null) {
            return null;
        } else {
            return e.getValue();
        }
    }

    /**
     * Push the submission and resources to the distributed cache.
     *
     * @param jobSubmission The job submission to which the job belongs.
     */
    public void pushJobSubmissionToCache(final JobSubmission jobSubmission) {
        if (getJobSubmission(jobSubmission.getId()) == null) {
            gobyCache.put(new Element("/gridgain/job-submissions/" + jobSubmission.getId(), jobSubmission));

            for (final Resource r : jobSubmission.getResources()) {
                pushResourceToCache(jobSubmission.getId(), r);
            }
        }
    }

    public Resource getResourceFromCache(final String jobSubmissionId, final Resource resource) {
        return (Resource) safeGet(ICBFilenameUtils.concatPathParts(
                "resource", jobSubmissionId, resource.getId()));
    }

    public void pushResourceToCache(final String jobSubmissionId, final Resource resource) {
        LOG.info("Pushing resource to cache " + resource);
        gobyCache.put(new Element(
                ICBFilenameUtils.concatPathParts("resource", jobSubmissionId, resource.getId()),
                resource));
    }

    public void storeJobOutput(final Job job, final Serializable jobOutput) {
        storeJobOutput(job, "default", jobOutput);
    }

    /**
     * Store a named job output.  Associate an output to a job and a name.
     *
     * @param job       The job that produced the output.
     * @param name      The name of the output being stored.  The name 'default' is reserved.
     * @param jobOutput The job output.
     */
    public void storeJobOutput(final Job job, final String name, final Serializable jobOutput) {
        final String newOutputKey = job.getTag() + "-output-" + name;
        if (safeGet(newOutputKey) == null) {
            gobyCache.put(new Element(newOutputKey, jobOutput));
        }
        final ArrayList<String> listOfOutputs = getOutputList(job);
        listOfOutputs.add(newOutputKey);
    }

    private ArrayList<String> getOutputList(final Job job) {
        final String outputListKey = job.getTag() + "-output-list";
        ArrayList<String> listOfOutputs = (ArrayList<String>) safeGet(outputListKey);
        if (listOfOutputs == null) {
            listOfOutputs = new ArrayList<String>();
            gobyCache.put(new Element(outputListKey, listOfOutputs));
        }

        return listOfOutputs;
    }

    /**
     * Return the output associated with the job. The job output named 'default' is returned.
     *
     * @param job Job that produced the output.
     * @return The job output when previously stored, or null if the output has not been produced.
     */
    public Serializable getJobOutput(final Job job) {
        return getJobOutput(job, "default");
    }

    public Serializable getJobOutput(final Job job, final String name) {
        return safeGet(job.getTag() + "-output-" + name);
    }

    /**
     * Remove the cached outputs from the cache.
     *
     * @param jobs set of jobs for which output should be removed.
     */
    public void clearJobOutput(final Job... jobs) {
        for (final Job job : jobs) {
            final ArrayList<String> listOfOutputs = getOutputList(job);
            for (final String outputKey : listOfOutputs) {
                // detach each output from the cache:
                gobyCache.remove(outputKey);
            }
        }
    }

    /**
     * Starts the distributed cache instance.
     * @throws GobyCacheException If start failed.
     */
    public void start() throws GobyCacheException {
        System.setProperty(CacheManager.ENABLE_SHUTDOWN_HOOK_PROPERTY, "true");
        final CacheManager cacheManager = CacheManager.getInstance();
        if (!cacheManager.getStatus().equals(Status.STATUS_ALIVE)) {
            throw new GobyCacheException("Cannot start gobyCache");
        }
        gobyCache = (Cache) cacheManager.getEhcache("goby-cache");
        LOG.info("EHCache started.");
    }

    /**
     * Stops JBoss Cache instance.
     */
    public void stop() {
        if (gobyCache != null) {
            CacheManager.getInstance().shutdown();
            System.out.println("gobyCache data manager stopped.");
            gobyCache = null;
        }
    }

    public void pushJobToCache(final Job job) {
        if (LOG.isDebugEnabled()) {
            LOG.debug("Storing job: " + job);
        }
        gobyCache.put(new Element(job.getTag(), job));
    }

    public Job getJobFromCache(final String jobTag) {
        return (Job) safeGet(jobTag);
    }

    public Serializable get(final String key) {
        return safeGet(key);
    }

    public void put(final String key, final Serializable value) {
        gobyCache.put(new Element(key, value));
    }

    /**
     * Returns the value of the grid active counter. This value is synchronized across all the nodes of the grid, and
     * increases by one every second or so.
     *
     * @return Value of the counter or zero if the counter cannot be retrieved.
     */
    public long getActiveCounterValue() {
        final long value;
        final Object activeCoutner = get(JobWithCounter.ACTIVE_COUNTER_NAME);
        if (activeCoutner != null) {
            value = ((JobWithCounter) activeCoutner).counter;
        } else {
            value = 0L;
        }
        return value;
    }
}
