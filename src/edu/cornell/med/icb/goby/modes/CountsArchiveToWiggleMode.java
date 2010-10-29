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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.counts.CountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.goby.counts.WiggleWindow;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
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
    private static final String MODE_NAME = "counts-to-wiggle";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Converts a full genome counts archive to "
            + "the Wiggle format.  The wiggle format can be imported in the UCSC genome "
            + "browser to visualize counts in the context of genome annotations."
            + "(See http://genome.ucsc.edu/goldenPath/help/wiggle.html)";

    /**
     * The input file.
     */
    private String inputBasename;

    /**
     * The output file.
     */
    private String outputFile;

    /**
     * If true only a subset of references will be processed defined by the set
     * {@link #includeReferenceNames}.
     */
    private boolean filterByReferenceNames;

    /**
     * The set of reference names to process if the user chose to filter.
     */
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();

    /**
     * The name to embed in the wiggle file.
     */
    private String label;

    /**
     * Use to switch from the default "counts" archive to another count archive within the
     * same basename.
     */
    private String alternativeCountArchiveExtension;

    /**
     * Used to set the span for the wiggle.
     */
    private int resolution;

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

        final String includeReferenceNameCommas = jsapResult.getString("include-reference-names");
        if (includeReferenceNameCommas != null) {
            includeReferenceNames = new ObjectOpenHashSet<String>();
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameCommas.split("[,]")));
            System.out.println("Will write wiggles for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }

        // default output file has extension ".wig" for compatibility with the
        // http://www.broadinstitute.org/igv/
        if (outputFile == null) {
            outputFile = inputBasename + (includeReferenceNameCommas != null
                    ? ("-" + includeReferenceNameCommas) : "-all") + ".wig";
        }

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
            final WiggleWindow wiggleWindow = new WiggleWindow(writer, resolution, 0);

            for (int referenceIndex = 0; referenceIndex < reader.getNumberOfIndices(); referenceIndex++) {
                String referenceId = backwards.getId(referenceIndex).toString();
                boolean processThisSequence = true;
                // prepare reference ID for UCSC genome browser.
                if ("MT".equalsIgnoreCase(referenceId)) {
                    // patch chromosome name for UCSC genome browser:
                    referenceId = "M";
                }

                // ignore c22_H2, c5_H2, and other contigs but not things like chr1 (mm9)
                if (referenceId.startsWith("c") && !referenceId.startsWith("chr")) {
                    processThisSequence = false;
                }

                // ignore NT_*
                if (referenceId.startsWith("NT_")) {
                    processThisSequence = false;
                }

                if (filterByReferenceNames && !includeReferenceNames.contains(referenceId)) {
                    processThisSequence = false;
                }

                if (processThisSequence) {
                    // prepend the reference id with "chr" if it doesn't use that already
                    final String chromosome;
                    if (referenceId.startsWith("chr")) {
                        chromosome = referenceId;
                    } else {
                        chromosome = "chr" + referenceId;
                    }

                    long sumCount = 0;
                    int numCounts = 0;

                    CountsReader counts = reader.getCountReader(referenceIndex);
                    int lastLength = 0;
                    int lastPosition = 0;
                    while (counts.hasNextTransition()) {
                        counts.nextTransition();
                        lastPosition = counts.getPosition();
                        lastLength = counts.getLength();
                    }
                    final int maxWritePosition = (lastPosition + lastLength - 1);
                    wiggleWindow.reset();
                    wiggleWindow.setMaxDataSize(maxWritePosition);

                    writer.printf("variableStep chrom=%s span=%d\n", chromosome, resolution);
                    counts = reader.getCountReader(referenceIndex);

                    while (counts.hasNextTransition()) {
                        counts.nextTransition();
                        final int length = counts.getLength();

                        final int count = counts.getCount();
                        final int position = counts.getPosition();

                        wiggleWindow.addData(position, length, count);

                        sumCount += count;
                        numCounts++;
                    }
                    wiggleWindow.finish();
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
     * @throws JSAPException error parsing
     * @throws IOException error parsing or executing.
     */

    public static void main(final String[] args) throws JSAPException, IOException {
        new CountsArchiveToWiggleMode().configure(args).execute();
    }
}
