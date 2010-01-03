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

import java.io.File;

/**
 * @author Fabien Campagne
 *         Date: May 21, 2009
 *         Time: 12:20:09 PM
 */
public class LocalFile extends StoredFile {

    @Override
    public long getSize() {
        return (new File(url)).length();
    }
}
