package edu.cornell.med.icb.nextgen.datamodel;

/**
 * The various types of status.
 */
public enum JobStatus {
    /**
     * The Job has been instanciated, but not submitted for execution.
     */
    NEW(0),
    /**
     * The job has been submitted for execution.
     */
    SUBMITTED(0),
    /**
     * The job is ready for execution, all its predecessors have completed execution successfully.
     */
    READY(1),
    /**
     * The job has started execution.
     */
    STARTED(2),
    /**
     * The job has failed execution. Either exception, error, or reasonFailed should be available.
     */
    FAILED(3),
    /**
     * The job has successfully completed execution.
     */
    COMPLETED(4);

    private final int encodedValue;

    private JobStatus(int encodedValue) {
        this.encodedValue = encodedValue;
    }
    
    public int getEncodedValue() {
        return encodedValue;
    }
}
