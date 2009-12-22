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

/**
 * Specifies the maqsupport jar file for all goby modes.
 *
 * @author Fabien Campagne
 * Date: Apr 7, 2009
 * Time: 3:14:02 PM
 */
public abstract class AbstractGobyMode extends AbstractCommandLineMode {
    protected AbstractGobyMode(final String jarFilename) {
        super(jarFilename);
    }

    protected AbstractGobyMode() {
        super("goby.jar");
    }

}
