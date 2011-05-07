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
import com.google.protobuf.CodedInputStream;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.GobyVersion;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * Converts a compact alignment to plain text.
 *
 * @author Fabien Campagne
 */
public class UpgradeMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "upgrade";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Upgrade goby files to a new version of Goby. We try to devise Goby format to avoid upgrade steps, but sometimes upgrading the data structures cannot be avoided. This tool converts data structures to the latest Goby format.";


    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;
    private boolean silent;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }


    /**
     * Configure.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        final String[] inputFiles = jsapResult.getStringArray("input");
        basenames = AlignmentReaderImpl.getBasenames(inputFiles);
        return this;

    }

    public void execute() throws IOException {
        for (String basename : basenames) {
            upgrade(basename);
        }
    }

    /**
     * Upgrade a Goby alignment as needed.
     * @param basename Basename of the alignment.
     */
    public void upgrade(String basename) {
        try {
            AlignmentReaderImpl reader = new AlignmentReaderImpl(basename, false);
            reader.readHeader();
            String version = reader.getGobyVersion();
           if (!silent) {
               System.out.printf("processing %s with version %s %n", basename, version);
           }
            if (GobyVersion.isOlder(version, "1.9.6")) {
                if (reader.isIndexed()) {
                    // we need to upgrade 1.9.5- alignment indices to the new indexing scheme implemented in 1.9.6+:
                    UpgradeTo1_9_6 upgrader = new UpgradeTo1_9_6();
                    upgrader.setSilent(silent);
                    upgrader.upgrade(basename, reader);
                }
            }
        } catch (IOException e) {
            System.err.println("Could not read alignment " + basename);
            e.printStackTrace();
        }
    }


    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */

    public static void main(final String[] args) throws JSAPException, IOException {
        new UpgradeMode().configure(args).execute();
    }

    public void setSilent(boolean silent) {
        this.silent=silent;
    }
}