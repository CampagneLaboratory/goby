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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import edu.cornell.med.icb.util.VersionUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Output the goby.jar version number to stdout.
 *
 * @author Kevin Dorff
 *         Date: Jan 27 2010
 */
public class VersionMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "version";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Display the version of Goby.";

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
     * @throws IOException   error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        return this;
    }

    @Override
    public void execute() throws IOException {
        String versionPrefix = "development";
        InputStream versionFileStream = this.getClass().getClassLoader().getResourceAsStream("Version.txt");
        if (versionFileStream != null) {
            BufferedReader reader = new BufferedReader(new InputStreamReader(versionFileStream));
            versionPrefix = reader.readLine();  // e.g., 2.0 for instance
        }
        final String version = VersionUtils.getImplementationVersion(GobyDriver.class);
        System.out.printf("Goby Version: %s %s%n", versionPrefix, version.replace("development ", ""));
    }
}
