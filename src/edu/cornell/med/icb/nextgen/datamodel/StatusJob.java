/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.nextgen.datamodel;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.Serializable;

/**
 * A status job so the web site can update.
 * Use the factory method reateStatusJobSubmission(...) to create of of these.
 *
 * @author Kevin Dorff
 */
public class StatusJob extends Job implements Serializable {

    private static final Log LOG = LogFactory.getLog(JobSubmission.class);

    private static final long serialVersionUID = -1461394355386240582L;

    final String statusText;

    public StatusJob(final String statusText) {
        newTag();
        this.statusText = statusText;
    }

    public String getStatusText() {
        return statusText;
    }

    /**
     * Create a JobSubmission that contains a JobStatus.
     * @param baseJobSubmission the JobSubmission we are creating status for.
     * @param statusText the status text to send
     * @return the JobSubmission containing a StatusJob.
     */
    public static JobSubmission createStatusJobSubmission(
            final JobSubmission baseJobSubmission, final String statusText) {
        if (baseJobSubmission == null) {
            LOG.error("Cannot create StatusJob, baseJobSubmission was null");
            return null;
        }
        if (baseJobSubmission.getId() == null) {
            LOG.error("Cannot create StatusJob, baseJobSubmission.id was null");
            return null;
        }
        if (baseJobSubmission.getOutputQueueHelper() == null) {
            LOG.error("Cannot create StatusJob, baseJobSubmission.outputQueueHelper was null");
            return null;
        }
        final JobSubmission jobSubmission = new JobSubmission();
        jobSubmission.id = baseJobSubmission.getId();
        jobSubmission.setOutputQueueHelper(baseJobSubmission.getOutputQueueHelper());
        final StatusJob statusJob = new StatusJob(statusText);
        jobSubmission.add(statusJob);
        return jobSubmission;
    }
}
