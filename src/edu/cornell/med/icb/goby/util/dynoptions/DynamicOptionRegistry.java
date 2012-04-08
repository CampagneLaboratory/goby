/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.util.dynoptions;

import com.martiansoftware.jsap.JSAPResult;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;

/**
 * A global registry for dynamic options. These options are set on the command line with options of the
 * form -x class-name:key=value
 *
 * @author Fabien Campagne
 *         Date: 3/11/12
 *         Time: 12:05 PM
 */
public class DynamicOptionRegistry {
    private DynamicOptionRegistry() {
    }

    public static void register(final DynamicOptionClient client) {
        if (!registeredDOClients.contains(client)) {
            registeredDOClients.add(client);
        }
    }

    static private ObjectAVLTreeSet<DynamicOptionClient> registeredDOClients = new ObjectAVLTreeSet<DynamicOptionClient>();
    private static String[] dymamicOptions;

    public static void parseCommandLineOptions(final JSAPResult jsapResult) {
        parseCommandLineOptions(jsapResult.getStringArray("dynamic-options"));

    }

    public static void parseCommandLineOptions(final String[] dynamicOptions) {
        // parse dynamic options:
        dymamicOptions = dynamicOptions;
        for (final String dymamicOption : dynamicOptions) {
            boolean parsed = false;
            for (final DynamicOptionClient doc : registeredDOClients) {

                if (doc.acceptsOption(dymamicOption)) {

                    parsed = true;
                    break;
                }
            }
            if (!parsed) {
                System.err.println("Error: none of the installed tools could parse dynamic option: " + dymamicOption);
                System.exit(1);
            }
        }
    }

    public static void printHelp() {
        System.out.println("The following dynamic options have been defined:");
        for (final DynamicOptionClient doc : registeredDOClients) {
            String[] keys = doc.getKeys();
            String[] helpMessages = doc.getHelpMessages();
            String[] defaultValues = doc.getDefaultValues();
            int i = 0;
            System.out.println("-x "+doc.getClassname()+":");
            for (String key : keys) {

                System.out.printf("  - %s: %s default: %s%n",  key, helpMessages[i], defaultValues[i]);
                i++;
            }


        }
    }
}
