/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionRegistry;
import edu.cornell.med.icb.util.VersionUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

import java.io.IOException;

/**
 * @author Fabien Campagne
 *         Date: Apr 7, 2009
 *         Time: 5:18:28 PM
 */
public class GobyDriver extends GenericToolsDriver {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(GobyDriver.class);

    private static final String DRIVER_JAR_NAME = "goby.jar";

    public GobyDriver() {
        super(DRIVER_JAR_NAME);
    }

    public static void main(final String[] args) throws IOException, JSAPException {
        final String version = VersionUtils.getImplementationVersion(GobyDriver.class);
        DynamicOptionRegistry.autoRegister();
        if (LOG.isDebugEnabled()) {
            LOG.debug(GobyDriver.class.getName() + " Implementation-Version: " + version);
            LOG.debug("Running with: " + ArrayUtils.toString(args));
        }
        int status = 0;
        try {
            new GobyDriver().configure(args).execute();
            status = 0;
        } catch (Exception e) {
            e.printStackTrace();
            status = 1;
        }
        System.exit(status);


    }
}
