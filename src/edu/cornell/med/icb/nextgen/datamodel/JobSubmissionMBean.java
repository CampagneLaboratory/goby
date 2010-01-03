package edu.cornell.med.icb.nextgen.datamodel;

/**
 * @author Fabien Campagne
 *         Date: Jul 26, 2009
 *         Time: 4:45:30 PM
 */
public interface JobSubmissionMBean {

    String getStatus();

    int getSubmittedJobs();

    int getReadyJobs();

    int getExecutingJobs();

    int getCompletedJobs();

    int getFailedJobs();

}
