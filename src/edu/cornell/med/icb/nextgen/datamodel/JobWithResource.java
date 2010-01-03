package edu.cornell.med.icb.nextgen.datamodel;

import java.io.Serializable;

/**
 * @author Fabien Campagne
*         Date: Aug 10, 2009
*         Time: 2:41:24 PM
*/
public class JobWithResource extends Job implements Serializable {

    private static final long serialVersionUID = -7453749703999790492L;

    public Resource resource;

    public JobWithResource(Resource resource) {
        super();
        this.resource = resource;
    }
}
