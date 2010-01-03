package edu.cornell.med.icb.nextgen.datamodel;

import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: Aug 19, 2009
 *         Time: 10:31:44 AM
 */
public class JobWithCounter extends Job implements Serializable {

    private static final long serialVersionUID = -4125249045667613350L;

    public static final String ACTIVE_COUNTER_NAME = "/active-counter";
    public long counter;
}
