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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;

import java.io.IOException;

/**
 * @author Fabien Campagne
 * Date: Apr 7, 2009
 * Time: 5:18:28 PM
 */
public class GobyDriver extends GenericToolsDriver {
    private static final String DRIVER_JAR_NAME = "goby.jar";

    public GobyDriver() {
        super(DRIVER_JAR_NAME);
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        new GobyDriver().configure(args).execute();
        System.exit(0);
    }
}
