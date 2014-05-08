/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.util.motifs.MotifMatcher;
import edu.cornell.med.icb.goby.util.motifs.SubSequenceMotifMatcher;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Test that Goby can connect to R via JRI/RJava.
 *
 * @author Fabien Campagne
 *         Date: May 8 2014
 *         Time: 11:36 AM
 */
public class TestRConnectionMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "test-r-connection";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Test that Goby can connect to R via JRI/RJava.";


    /**
     * {@inheritDoc}
     */
    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    private static final Logger LOG = Logger.getLogger(TestRConnectionMode.class);

    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException                    error parsing
     * @throws com.martiansoftware.jsap.JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);
        return this;
    }

    @Override
    public void execute() throws IOException {
        final Rengine rEngine = GobyRengine.getInstance().getRengine();
        boolean connectionWorks = rEngine != null && rEngine.isAlive();
        if (connectionWorks) {
            System.out.println("R connection is working properly");
            System.exit(0);
        } else {
            System.out.println("Could not establish connection. Test failed.");
            System.exit(1);
        }

    }


    public static void main
            (
                    final String[] args) throws IOException, JSAPException {
        new TestRConnectionMode().configure(args).execute();
    }
}