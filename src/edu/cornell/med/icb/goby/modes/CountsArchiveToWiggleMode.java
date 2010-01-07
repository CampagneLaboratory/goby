/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.counts.CountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.zip.GZIPOutputStream;

/**
 * Converts a full genome counts archive to the Wiggle format.
 * The  <a href="http://genome.ucsc.edu/goldenPath/help/wiggle.html">wiggle format</a> can
 * be imported in the UCSC genome browser to visualize counts in the context of genome annotations.
 *
 * @author Fabien Campagne
 */
public class CountsArchiveToWiggleMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "counts-to-wiggle";
    public static final String MODE_DESCRIPTION = "Converts a full genome counts archive to the Wiggle format.  The wiggle format can be imported in the UCSC genome browser to visualize counts in the context of genome annotations. (http://genome.ucsc.edu/goldenPath/help/wiggle.html)";


    /**
     * The input file.
     */
    private String inputBasename;

    /**
     * The output file.
     */
    private String outputFile;

    boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();
    private String label;
    /**
     * Use to switch from the default "counts" archive to another count archive within the same basename.
     */
    private String alternativeCountArchiveExtension;
    private int resolution;

    public String getModeName() {
        return MODE_NAME;
    }

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
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputBasename = AlignmentReader.getBasename(jsapResult.getString("input"));
        outputFile = jsapResult.getString("output");
        resolution = jsapResult.getInt("resolution");
        alternativeCountArchiveExtension = jsapResult.getString("alternative-count-archive");
        label = jsapResult.getString("label");
        if (label == null) {
            label = new File(inputBasename + "-" + alternativeCountArchiveExtension).getName();
            System.out.println("Setting label from basename: " + label);
        }
        final String includeReferenceNameComas = jsapResult.getString("include-reference-names");
        if (includeReferenceNameComas != null) {
            includeReferenceNames = new ObjectOpenHashSet<String>();
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameComas.split("[,]")));
            System.out.println("Will write wiggles for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }
        if (outputFile == null)
            outputFile = FilenameUtils.getBaseName(inputBasename) +".wiggle"+(includeReferenceNameComas != null ? ("-"+includeReferenceNameComas) :
                    "-all") ;

        return this;
    }

    /**
     * Run the map2text mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(new GZIPOutputStream(new FastBufferedOutputStream(new FileOutputStream(outputFile + ".gz"))));
            writer.write("track type=wiggle_0 name=" + label + " visibility=full viewLimits=1:200\n");
            final AlignmentReader alignment = new AlignmentReader(inputBasename);
            alignment.readHeader();
            alignment.close();
            final IndexedIdentifier referenceIds = alignment.getTargetIdentifiers();
            final DoubleIndexedIdentifier backwards = new DoubleIndexedIdentifier(referenceIds);
            final CountsArchiveReader reader = new CountsArchiveReader(inputBasename, alternativeCountArchiveExtension);
            final int resolution = this.resolution;

            for (int referenceIndex = 0; referenceIndex < reader.getNumberOfIndices(); referenceIndex++) {
                String referenceId = backwards.getId(referenceIndex).toString();
                boolean processThisSequence = true;
                // prepare reference ID for UCSC genome browser.
                if ("MT".equalsIgnoreCase(referenceId)) {
                    // patch chromosome name for UCSC genome browser:
                    referenceId = "M";
                }
                if (referenceId.startsWith("c")) {
                    // ignore c22_H2, c5_H2, and other contigs..
                    processThisSequence = false;
                }

                if (filterByReferenceNames && !includeReferenceNames.contains(referenceId)) {
                    processThisSequence = false;
                }
                if (processThisSequence) {
                    writer.printf("variableStep chrom=chr%s span=%d\n", referenceId, resolution);
                    long sumCount = 0;
                    int numCounts = 0;
                    final CountsReader counts = reader.getCountReader(referenceIndex);
                    int position = 0;
                    while (counts.hasNextTransition()) {
                        counts.nextTransition();
                        final int length = counts.getLength();

                        final int count = counts.getCount();
                        position = counts.getPosition();
                        if (count != 0) {
                            for (int i = 0; i < length; i += resolution) {

                                // positions start at 1 for UCSC genome browser:
                                writer.printf("%d %d\n", position + i + 1, count);


                            }
                        }
                        sumCount += count;
                        numCounts++;
                        //    position += length;
                    }
                    final double averageCount = sumCount / (double) numCounts;
                    System.out.println("average count for sequence " + referenceId + " " + averageCount);
                }
            }
        } finally {
            IOUtils.closeQuietly(writer);
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
        new CountsArchiveToWiggleMode().configure(args).execute();
    }
}
