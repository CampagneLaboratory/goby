package edu.cornell.med.icb.goby.methylation;

import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: Oct 24, 2010
 *         Time: 11:39:49 AM
 */
public class MethylationSite implements Serializable {
    public int chromosome;
    public char strand;
    public int position;
    public int methylatedReadCount;
    public int totalCount;


    public float getMethylationRate() {
        return ((float)methylatedReadCount/ (float)totalCount);
    }
   // public static final long serialVersionUID = 5664745795898488209L;

}
