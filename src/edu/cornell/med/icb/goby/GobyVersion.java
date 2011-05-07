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

import java.util.Date;
import java.text.SimpleDateFormat;

/**
 * Helper class to compare Goby versions.
 *
 * @author Fabien Campagne
 *         Date: May 6, 2011
 *         Time: 6:29:16 PM
 */
public class GobyVersion {
    private static String[] versionPairs = {
            "1.9.5", "20110501000000",
            "1.9.6", "20110506200142", // TODO update to correct compile time.
            "1.9.6+", now()
    };

    private static String now() {
            Date d=new Date();
       SimpleDateFormat formatter=new SimpleDateFormat("yyyymmddHHmmss");

        return formatter.format(d);
    }

    /**
     * Return true when version is strictly older than toBeComparedl, false otherwise.
     *
     * @param version
     * @param toBeCompared
     * @return
     */
    public static boolean isOlder(String version, String toBeCompared) {

        if (version.equals(toBeCompared)) {
            // same version, the second one is not older.
            return false;
        }
        if (version.indexOf(".") >= 0 && toBeCompared.indexOf(".") > 0) {
            if (version.equals("1.9.5-") && version.equals("1.9.6")) return true;
        } else {

            long time1 = Long.parseLong(reduce(version));
            long time2 = Long.parseLong(reduce(toBeCompared));
            return time1 < time2;
        }
        return false;
    }

    private static String reduce(String version) {
        if (version.startsWith("development")) {
            // format "development (version)"
            String tokens[] = version.split("[()]");
            final String token = tokens[1];
            return token;
        } else {
            for (int i = 0; i < versionPairs.length; i++) {
                if (version.equals(versionPairs[i])) return versionPairs[i + 1];
            }
        }
        return "unknown";
    }
}
