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

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.algorithm.ExportSplicingEvents;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;

/**
 * @author Fabien Campagne
 *         Date: 9/14/12
 *         Time: 10:48 AM
 */
public class ExtractSplicingEventsMode extends AbstractGobyMode {
    /**
     * Used to log informational and debug messages.
     */
    private static final Log LOG = LogFactory.getLog(ExtractSplicingEventsMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "extract-splicing-events";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Extract splicing events from an alignment file.";
    private String[] inputFilenames;
    private String outputFilename;
    private int qualThreshold;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    @Override
    public AbstractCommandLineMode configure(String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilenames = jsapResult.getStringArray("input");
        outputFilename = jsapResult.getString("output");
        qualThreshold = jsapResult.getInt("min-mapping-quality");

        return this;
    }

    @Override
    public void execute() throws IOException {
        Writer output = null;
        final boolean consoleOutput = outputFilename.equals("-");
        if (consoleOutput) {
            output = new OutputStreamWriter(System.out);
        } else {
            output = new FileWriter(outputFilename);
        }

        for (String filename : inputFilenames) {
            ExportSplicingEvents processor = new ExportSplicingEvents(output);
            processor.setMinMappingQuality(qualThreshold);
            processor.process(filename);

        }
        if (!consoleOutput) {
            output.close();
        }
    }
}
