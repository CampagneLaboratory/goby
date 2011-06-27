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
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;

/**
 * A simple mode to convert a set of filename in the format attr1 delimiter attr2 delimiter ... .ext to a tab
 * delimited attribute file suitable for loading in IGV.
 * This mode is useful if you follow the following convention when naming a set of files"
 * TAG-phenotype-sampleId.counts  you can then extract a sample attribute table with these arguments: -a ignore -a phenotype -a sampleId -d - -o sample-attributes.txt *.counts
 *
 * @author Fabien Campagne
 *         Date: 6/11/11
 *         Time: 5:03 PM
 */
public class FilesToAttributesMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "files-to-attributes";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Scan a list of files and construct an IGV attribute file with attributes about each file.";

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(FilesToAttributesMode.class);
    private String delimiter;
    private String[] attributeNames;
    private String[] filenames;
    private String outputFilename;
    private String suffixTransform;
    private String suffixLookup;
    private String suffixReplace;


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
        delimiter = jsapResult.getString("delimiter");
        attributeNames = jsapResult.getStringArray("attribute");
        filenames = jsapResult.getStringArray("file");
        outputFilename = jsapResult.getString("output");
        suffixTransform = jsapResult.getString("suffix");
        if (suffixTransform != null) {
            String toks[] = suffixTransform.split("/");
            suffixLookup=toks[0];
           if (toks.length>1) {
               suffixReplace=toks[1];
           }   else {
               suffixReplace="";
           }
        }
        return this;
    }

    @Override
    public void execute() throws IOException {
        final PrintWriter out = outputFilename == null ? new PrintWriter(System.out) :
                new PrintWriter(outputFilename);
        try {
            int index = 0;
            out.print("trackName\t");
            for (String attributeName : attributeNames) {
                if (!"ignore".equals(attributeName)) {
                    out.print(attributeName);

                    if (index != attributeNames.length) {
                        out.print('\t');
                    }
                }
                index++;
            }
            out.println();
            for (String file : filenames) {
                final String file1 = file;
                String name = FilenameUtils.getName(file1);
                String filename = FilenameUtils.getBaseName(name);
                String[] tokens = filename.split(delimiter);
                out.printf("%s\t", adjustSuffix(name));
                for (int i = 0; i < attributeNames.length; i++) {
                    if (!"ignore".equals(attributeNames[i])) {
                        out.print(tokens[i]);
                        if (i != attributeNames.length) {
                            out.print('\t');
                        }
                    }

                }
                out.println();
            }
            out.flush();
        } finally {
            out.close();

        }
    }

    private String adjustSuffix(final String name) {

        if (name.endsWith(suffixLookup)) {
            return name.replace(suffixLookup, suffixReplace);
        }
        return name;
    }
}
