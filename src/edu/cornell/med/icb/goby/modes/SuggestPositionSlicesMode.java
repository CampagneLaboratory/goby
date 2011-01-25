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
import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang.ArrayUtils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Converts a compact alignment to plain text.
 *
 * @author Fabien Campagne
 */
public class SuggestPositionSlicesMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "suggest-position-slices";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Suggest how to slice an alignment by position to yield roughly equally sized slices.";

    /**
     * The output file.
     */
    private String outputFilename;

    /**
     * The basename of the compact alignment.
     */
    private String[] basenames;
    private int modulo;
    private int numberOfSlices;
    private int bufferLength;


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
        basenames = AlignmentReader.getBasenames(inputFiles);
        outputFilename = jsapResult.getString("output");
        modulo = jsapResult.getInt("modulo");
        numberOfSlices = jsapResult.getInt("number-of-slices");



        return this;
    }


    /**
     * Suggests slices to process a large alignment file in parallel.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintStream stream = null;
        try {
            stream = outputFilename == null ? System.out
                    : new PrintStream(new FileOutputStream(outputFilename));

            ConcatSortedAlignmentReader input = new ConcatSortedAlignmentReader(basenames);
            input.readHeader();
            DoubleIndexedIdentifier ids = new DoubleIndexedIdentifier(input.getTargetIdentifiers());
            ObjectList<ReferenceLocation> locations = input.getLocations(modulo);
            stream.println("targetIdStart\t%positionStart\tstart:(ref,pos)\ttargetIdEnd\t%positionEnd\tend:(ref,pos)");
            ReferenceLocation[] breakpoints = new ReferenceLocation[numberOfSlices +1];
            for (int i = 0; i < numberOfSlices - 1; i++) {
                breakpoints[i+1] = locations.get(locations.size() / (numberOfSlices - 1) * i);
            }
            breakpoints[0]=new ReferenceLocation(0,0);
            // largest position in the last reference sequence:
            final int lastTargetIndex = ids.size() - 1;
            breakpoints[breakpoints.length-1]=new ReferenceLocation(lastTargetIndex,input.getTargetLength(lastTargetIndex));

            for (int i = 0; i < numberOfSlices ; i++) {

                stream.printf(String.format("%s\t%d\t%s,%d\t%s\t%d\t%s,%d%n",
                        ids.getId(breakpoints[i].targetIndex),
                        breakpoints[i].position,
                        ids.getId(breakpoints[i].targetIndex),
                        breakpoints[i].position,
                        ids.getId(breakpoints[i + 1].targetIndex),
                        breakpoints[i+1].position,
                        ids.getId(breakpoints[i + 1].targetIndex),
                        breakpoints[i+1].position));

            }

        } finally {
            if (stream != System.out) {
                IOUtils.closeQuietly(stream);
            }
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
        new SuggestPositionSlicesMode().configure(args).execute();
    }
}