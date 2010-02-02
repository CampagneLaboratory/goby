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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

/**
 * Callback handler for R to log information using {@link org.apache.commons.logging.Log}.
 */
class RLoggerMainLoopCallback implements RMainLoopCallbacks {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(RLoggerMainLoopCallback.class);

    /**
     * Called when R prints output to the console.
     *
     * @param rengine calling engine
     * @param text text to display in the console
     * @param type output type (0=regular, 1=error/warning)
     */
    public void rWriteConsole(final Rengine rengine, final String text, final int type) {
        if (type == 0) {
            LOG.info(text);
        } else {
            LOG.warn(text);
        }
    }

    /**
     * Called when R enters or exits a longer evaluation. It is usually a good idea to signal
     * this state to the user, e.g. by changing the cursor to a "hourglass" and back.
     *
     * @param rengine calling engine
     * @param which identifies whether R enters (1) or exits (0) the busy state
     */
    public void rBusy(final Rengine rengine, final int which) {
        LOG.info(rengine.getName()
                + (which == 0 ? " entering " : " exiting ") + "busy state");
    }

    /**
     * Called when R waits for user input. During the duration of this callback it is safe to
     * re-enter R, and very often it is also the only time. The implementation is free to block
     * on this call until the user hits Enter, but in JRI it is a good idea to call
     * {@link org.rosuda.JRI.Rengine#rniIdle()} occasionally to allow other event handlers
     * (e.g graphics device UIs) to run. Implementations should NEVER return immediately even
     * if there is no input - such behavior will result in a fast cycling event loop which makes
     * the use of R pretty much impossible.
     *
     * @param rengine calling engine
     * @param prompt prompt to be displayed at the console prior to user's input
     * @param addToHistory flag telling the handler whether the input should be considered
     * for adding to history (!=0) or not (0)
     * @return user's input to be passed to R for evaluation
     */
    public String rReadConsole(final Rengine rengine, final String prompt, final int addToHistory) {
        return null;
    }

    /**
     * Called when R want to show a warning/error message (not to be confused with messages
     * displayed in the console output).
     *
     * @param rengine calling engine
     * @param message message to display
     */
    public void rShowMessage(final Rengine rengine, final String message) {
        LOG.error(message);
    }

    /**
     * Called when R expects the user to choose a file.
     *
     * @param rengine calling engine
     * @param newFile flag determining whether an existing or new file is to be selecteed
     * @return path/name of the selected file
     */
    public String rChooseFile(final Rengine rengine, final int newFile) {
        return null;
    }

    /**
     * Called when R requests the console to flush any buffered output.
     *
     * @param rengine calling engine
     */
    public void rFlushConsole(final Rengine rengine) {
    }

    /**
     * Called to save the contents of the history (the implementation is responsible of keeping
     * track of the history).
     *
     * @param rengine calling engine
     * @param filename name of the history file
     */
    public void rSaveHistory(final Rengine rengine, final String filename) {
    }

    /**
     * Called to load the contents of the history.
     *
     * @param re calling engine
     * @param filename name of the history file
     */
    public void rLoadHistory(final Rengine re, final String filename) {
    }
}
