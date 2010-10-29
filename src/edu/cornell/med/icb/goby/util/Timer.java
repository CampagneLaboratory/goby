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

package edu.cornell.med.icb.goby.util;

import java.util.Date;

/**
 * @author Fabien Campagne
 *         Date: Oct 28, 2010
 *         Time: 12:16:30 PM
 */
public class Timer {
    long startTime;
    long stopTime;

    public void start() {
        startTime = new Date().getTime();
    }

    public void stop() {
        stopTime = new Date().getTime();
    }

    public String toString() {
        long tmp = stopTime - startTime;
        long milliseconds = tmp / 1000;
        long secs = milliseconds / 60;
        long mins = secs / 60;
        long hrs = mins / 60;
        return String.format("Time elapsed %d hrs %d mins %d sec %d ms",
                hrs,
                mins,
                secs,
                milliseconds);
    }

    public long millis() {
        long tmp = stopTime - startTime;
        long milliseconds = tmp / 1000;
        return milliseconds;
    }

    public long seconds() {
        long tmp = stopTime - startTime;
        long milliseconds = tmp / 1000;
        long secs = milliseconds / 60;
        return secs;
    }

    public long minutes() {
        long tmp = stopTime - startTime;
        long milliseconds = tmp / 1000;
        long secs = milliseconds / 60;
        long mins = secs / 60;
        return mins;
    }
}
