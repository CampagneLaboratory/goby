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

import java.io.PrintStream;
import java.io.Serializable;

/**
 * A job which transforms inputs into outputs. Jobs can be executed by calling their execute method, or run on a grid
 * by submitting the job to the grid. Used primarily by TestJobSubmissionQueueHelper.
 *
 * @author Kevin Dorff
 * @see edu.cornell.med.icb.nextgen.datamodel.JobSubmission
 */
public class EchoJob extends Job implements Serializable{

    private static final long serialVersionUID = 7829536947150591204L;

    /**
     * The PrintStream to write output to in execute(). This is Transient because jobs
     * WILL be serialized and you cannot serialize a PrintStream. output is only relevant
     * when we actually go to execute().
     */
    private transient PrintStream output = System.out;

    /**
     * Construct a new job. The job is assigned a random tag.
     */
    public EchoJob() {
        super();
    }

    /**
     * Set the output (effects where execute() echo's the value to.
     * @param output the new PrintStream to write output to in execute().
     */
    public void setOutput(final PrintStream output) {
        this.output = output;
    }

    /**
     * Check if this Job equals another object (probably a Job).
     * @param other the other object to compare this job to
     * @return true if they are the same job (same intTag).
     */
    @Override
    public boolean equals(Object other) {
        if (!(other instanceof EchoJob)) return false;
        EchoJob otherJob = (EchoJob) other;
        return otherJob.getIntTag() == getIntTag();
    }

    /**
     * Override this method to implement job specific processing.
     *
     * @return The job itself. Job output should be stored in the cache.
     * @throws InterruptedException not used in this class
     */
    public Serializable execute() throws InterruptedException {
        PrintStream out;
        if (output != null) {
            out = output;
        } else {
            out = System.out;
        }
        out.print(input.toString() + " ");
        return this;
    }
}
