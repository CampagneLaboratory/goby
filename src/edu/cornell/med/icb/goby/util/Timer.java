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
 * A simple timer to time operations.
 *
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
        Elapsed elapsed = new Elapsed(startTime, stopTime);


        long hrs = elapsed.removeHours();
        long mins = elapsed.removeMinutes();
        long secs = elapsed.removeSeconds();
        long milliseconds = elapsed.removeMilliseconds();
        return String.format("Time elapsed %d hrs %d mins %d sec %d ms",
                hrs,
                mins,
                secs,
                milliseconds);
    }

    public long millis() {

        long milliseconds = stopTime - startTime;
        return milliseconds;
    }

    public long seconds() {
        long milliseconds = stopTime - startTime;
        long secs = milliseconds / 60;
        return secs;
    }

    public long minutes() {
        long milliseconds = stopTime - startTime;
        long secs = milliseconds / 60;
        long mins = secs / 60;
        return mins;
    }

    private class Elapsed {
        long duration;

        private static final int SECONDS_IN_MS = 1000;
        private static final int MINUTES_IN_MS = SECONDS_IN_MS * 60;
        private static final int HOURS_IN_MS = MINUTES_IN_MS * 60;

        public Elapsed(long startTime, long stopTime) {
            duration = stopTime - startTime;
        }

        int removeAmount(long amountInMilliSeconds) {

            long amount = duration / amountInMilliSeconds;
            if (amount > 0) {
                duration -= amount * amountInMilliSeconds;
            }
            return (int) amount;
        }

        int removeHours() {
            return removeAmount(HOURS_IN_MS);
        }

        int removeMinutes() {
            return removeAmount(MINUTES_IN_MS);
        }

        int removeSeconds() {
            return removeAmount(SECONDS_IN_MS);

        }

        int removeMilliseconds() {
            return removeAmount(1);

        }
    }
}
