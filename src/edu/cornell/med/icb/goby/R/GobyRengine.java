/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.R;

import org.apache.commons.lang.SystemUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.Rengine;

/**
 * Main interface to R. This requires native R libraries and rJava to be installed. From the R
 * console enter:
 * <p>&nbsp;&nbsp;<em>install.packages('rJava')</em></p>
 * <p>When running the Java code, you need to add the R and rJava libraries to the library path.
 * For example on windows, add
 * <em>-Djava.library.path="C:\Program Files (x86)\R\R-2.10.1\library\rJava\jri"</em>.  On
 * Unix, add the R and JRI paths to the <em>LD_LIBRARY_PATH</em> environment variable.
 * <p>
 * See <a href="http://www.r-project.org/">The R Project for Statistical Computing</a> and
 * <a href="http://www.rforge.net/rJava/">rJava</a> for reference.
 */
public final class GobyRengine {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(GobyRengine.class);

    /**
     * Thread that runs the R engine.
     */
    private static final GobyRengine INSTANCE = new GobyRengine();

    /**
     * Rengine for the application.
     */
    private Rengine rengine;

    /**
     * Create a new goby.
     */
    private GobyRengine() {
        super();

        try {
            // tell REngine not to shutdown the jvm if the native R library cannot be loaded
            System.setProperty("jri.ignore.ule", "yes");

            // Tell R to be verbose if we are debugging
            if (LOG.isDebugEnabled()) {
                Rengine.DEBUG = 42;
                LOG.debug("java.library.path: " + SystemUtils.JAVA_LIBRARY_PATH);
            }

            // just making sure we have the right version of everything
            if (!Rengine.versionCheck()) {
                LOG.warn("Rengine cannot be initialized - java files don't match library version.");
                rengine = null;
                return;
            }

            rengine = Rengine.getMainEngine();
            if (rengine == null) {
                // NOTE: Do not use the default Rengine constructor
                rengine = new Rengine(new String[] {"--no-save"},
                        false, new RLoggerMainLoopCallback());
                if (!rengine.waitForR()) {       // will return false if R is dead
                    LOG.warn("Cannot load R");
                    rengine = null;
                }
            }
        } catch (UnsatisfiedLinkError e) {
            LOG.warn("Rengine libraries can not be found", e);
            rengine = null;
        }

        addShutdownHook();
    }

    /**
     * Add a shutdown hook so that the R thread is terminated cleanly on JVM exit.
     */
    private void addShutdownHook() {
        LOG.debug("Adding shutdown hook");
        Runtime.getRuntime().addShutdownHook(
                new Thread(GobyRengine.class.getSimpleName() + "-ShutdownHook") {  // NOPMD
                    @Override
                    public void run() {
                        LOG.info("Shutdown hook is terminating R");
                        if (rengine != null) {
                            rengine.end();
                        }
                    }
                });
    }

    /**
     * Return the main R "engine".
     * @return The interface to R which may be null if R is not available.
     */
    public Rengine getRengine() {
        return Rengine.getMainEngine();
    }

    /**
     * Get the singleton instance of the Rengine for Goby.
     * @return The instance.
     */
    public static GobyRengine getInstance() {
        return INSTANCE;
    }
}
