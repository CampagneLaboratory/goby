/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This file is part of the Goby IO API.
 *
 *     The Goby IO API is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     The Goby IO API is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with the Goby IO API.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.GobyVersion;
import edu.cornell.med.icb.goby.alignments.*;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.logging.ProgressLogger;

import java.io.IOException;

/**
 * Upgrade goby files to a new version of Goby. We try to devise Goby format to avoid upgrade steps, but sometimes
 * upgrading the data structures cannot be avoided (e.g., when we fix bugs that existed in earlier versions).
 * This tool converts data structures to the latest Goby format.
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
    private static final String MODE_DESCRIPTION = "Upgrade goby files to a new version of Goby. We try to devise Goby format to avoid upgrade steps, but sometimes upgrading the data structures cannot be avoided (e.g., when we fix bugs that existed in earlier versions). This tool converts data structures to the latest Goby format.";

    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;
    private boolean silent;
    private boolean check;


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
        //   check = jsapResult.getBoolean("check");
        return this;


    }

    public void execute() throws IOException {
        for (String basename : basenames) {
            upgrade(basename);
            if (check) {
                check(basename);
            }
        }
    }

    /**
     * Upgrade a Goby alignment as needed.
     *
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
            if (GobyVersion.isOlder(version, "goby_1.9.6")) {
                if (reader.isIndexed()) {
                    // we need to upgrade 1.9.5- alignment indices to the new indexing scheme implemented in 1.9.6+:
                    UpgradeTo1_9_6 upgrader = new UpgradeTo1_9_6();
                    upgrader.setSilent(silent);
                    upgrader.upgrade(basename, reader);
                }
            }
            if (GobyVersion.isOlder(version, "goby_1.9.8.2")) {
                if (reader.isIndexed()) {
                    // we need to upgrade 1.9.5- alignment indices to the new indexing scheme implemented in 1.9.6+:
                    UpgradeTo1_9_8_2 upgrader = new UpgradeTo1_9_8_2();
                    upgrader.setSilent(silent);
                    upgrader.upgrade(basename, reader);
                }
            }
        } catch (IOException e) {
            System.err.println("Could not read alignment " + basename);
            e.printStackTrace();
        }
    }

    public void check(String basename) {
        try {
            AlignmentReaderImpl reader = new AlignmentReaderImpl(basename, false);
            reader.readHeader();
            String version = reader.getGobyVersion();
            if (!silent) {
                System.out.printf("processing %s with version %s %n", basename, version);
            }
            if (GobyVersion.isMoreRecent(version, "1.9.6")) {
                if (reader.isIndexed()) {
                    ObjectList<ReferenceLocation> locations = reader.getLocations(1000);
                    System.out.println("Checking..");
                    ProgressLogger progress = new ProgressLogger();
                    progress.expectedUpdates = locations.size();
                    //  progress.priority = Level.INFO;
                    progress.start();
                    for (ReferenceLocation location : locations) {
                        Alignments.AlignmentEntry entry = reader.skipTo(location.targetIndex, location.position);
                        if (entry == null) {
                            System.err.printf("Entry must be found at position (t=%d,p=%d) %n", location.targetIndex,
                                    location.position);
                            System.exit(1);
                        }
                        if (entry.getTargetIndex() < location.targetIndex) {
                            System.err.printf("Entry must be found on reference >%d for position (t=%d,p=%d) %n",
                                    location.targetIndex, location.targetIndex,
                                    location.position);
                            System.exit(1);
                        }
                        if (entry.getPosition() < location.position) {
                            System.err.printf("Entry must be found at position >=%d for position (t=%d,p=%d) %n",
                                    location.position, entry.getTargetIndex()
                                    ,
                                    entry.getPosition());
                            System.exit(1);
                        }
                        progress.lightUpdate();
                    }
                    progress.stop();
                    System.out.printf("Checked %d skipTo calls", locations.size());
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
        this.silent = silent;
    }
}