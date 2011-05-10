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

package edu.cornell.med.icb.goby;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Helper class to compare Goby versions.
 *
 * @author Fabien Campagne
 *         Date: May 6, 2011
 *         Time: 6:29:16 PM
 */
public class GobyVersion {
    private static String[] versionPairs = {
            "1.9.5-", "20110101000000",
            "goby_1.9.6", "20110510111716"
    };

    private static String now() {
        Date d = new Date();
        SimpleDateFormat formatter = new SimpleDateFormat("yyyymmddHHmmss");

        return formatter.format(d);
    }

    /**
     * Return true when version is more recent or as recent as toBeComparedTo, false otherwise.
     *
     * @param version       The goby version under test.
     * @param toBeCompareTo A specific version we wish to compare version to.
     * @return rrue if version is more or as recent as toBeComparedTo.
     */
    public static boolean isMoreRecent(String version, String toBeCompareTo) {
        return !isOlder(version, toBeCompareTo);
    }

    /**
     * Return true when version is strictly older than toBeComparedTo, false otherwise.
     *
     * @param version       The goby version under test.
     * @param toBeCompareTo A specific version we wish to compare version to.
     * @return Return true when version is strictly older than toBeComparedTo, false otherwise.
     */
    public static boolean isOlder(String version, String toBeCompareTo) {

        if (version.equals(toBeCompareTo)) {
            // same version, the second one is not older.
            return false;
        }

        long time1 = Long.parseLong(reduce(version));
        long time2 = Long.parseLong(reduce(toBeCompareTo));
        return time1 < time2;

    }

    /**
     * Reduce a version string to a date in the format yyyymmddHHmmss.
     * @param version String obtained from the Goby jar manifest.
     * @return Date the version was packaged.
     */
    private static String reduce(String version) {
        if (version.indexOf('(')>=0 && version.indexOf(')')>=0)  {

            // format "development (DATE)"  OR "goby_1.9.6 (DATE)
            String tokens[] = version.split("[()]");
            final String token = tokens[1];
            return token;
        } else {
            for (int i = 0; i < versionPairs.length; i++) {
                if (version.equals(versionPairs[i])) return versionPairs[i + 1];
            }
        } // the version number was not recognized, assume we are dealing with a more recent version.
        LOG.warn(String.format("Version number %s not recognized. Assuming this version is the most recent.", version));
        return now();
    }

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(GobyVersion.class);

}
