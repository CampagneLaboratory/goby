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

package edu.cornell.med.icb.goby.util;

import org.apache.commons.logging.Log;

import java.util.logging.Logger;

/**
 * Counter to print a specific message up to maxWarnings.
 *
 * @author Fabien Campagne
 *         Date: Apr 16, 2011
 *         Time: 10:07:28 PM
 */
public class WarningCounter {
    int counter;
    private int maxWarnings;

    /**
     * Initialize this warning counter to warn at most maxWarnings times.
     *
     * @param maxWarnings maximum number of times the warning will be issued.
     */
    public WarningCounter(int maxWarnings) {
        this.maxWarnings = maxWarnings;
    }

    public WarningCounter() {
        this.maxWarnings = 10;
    }

    public boolean warnAgain() {
        return ++counter < maxWarnings;
    }

    public void warn(org.apache.commons.logging.Log log, String format, Object... option) {
        if (warnAgain()) {
            log.warn(String.format(format, option));
        }
    }

    public void warn(org.apache.log4j.Logger log, String format, Object... option) {
        if (warnAgain()) {
            log.warn(String.format(format, option));
        }
    }
}
