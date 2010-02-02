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

import org.apache.log4j.Logger;
import org.rosuda.JRI.Rengine;

/**
 * Main interface to R.
 */
public final class GobyRengine {
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(GobyRengine.class);

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

            if (LOG.isDebugEnabled()) {
                Rengine.DEBUG = 42;
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
                rengine = new Rengine(null, false, new RConsoleMainLoopCallback());
                if (!rengine.waitForR()) {
                    LOG.warn("Cannot load R");
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

    public static GobyRengine getInstance() {
        return INSTANCE;
    }
}
