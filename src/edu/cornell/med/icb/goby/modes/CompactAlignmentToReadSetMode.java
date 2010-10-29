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
import edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.reads.ReadSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.io.File;
import java.io.IOException;

/**
 * Converts a compact alignment to a read set file.
 *
 * @author Fabien Campagne
 */
public class CompactAlignmentToReadSetMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "alignment-to-read-set";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Converts a compact alignment to a read set file.";

    /**
     * The output read set suffix.
     */
    private String suffix;

    private String[] basenames;
    private boolean matchingReads;
    private boolean nonMatchingReads;
    private boolean ambiguousReads;
    private boolean nonAmbiguousReads;
    private int k;
    private ReadSet preFilter;

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

        final String[] inputFiles = jsapResult.getStringArray("input");
        final ObjectSet<String> basenameSet = new ObjectOpenHashSet<String>();

        for (final String inputFile : inputFiles) {
            basenameSet.add(AlignmentReader.getBasename(inputFile));
        }

        basenames = basenameSet.toArray(new String[basenameSet.size()]);
        this.matchingReads = jsapResult.getBoolean("matching-reads");
        this.nonMatchingReads = jsapResult.getBoolean("non-matching-reads");


        this.ambiguousReads = jsapResult.getBoolean("ambiguous-reads");
        this.nonAmbiguousReads = jsapResult.getBoolean("non-ambiguous-reads");
        if (!ambiguousReads && !nonAmbiguousReads) {
            System.err.println("Combination of options will result in no reads selected since a read is either ambiguous or is not.");
            System.exit(1);
        }
        this.k = jsapResult.getInt("ambiguity-threshold", 2);
        suffix = jsapResult.getString("suffix");
        final File preFilterFile = jsapResult.getFile("pre-filter");
        preFilter = new ReadSet();
        if (preFilterFile == null) {
            preFilter = null;
        } else {
            preFilter.load(preFilterFile);
        }

        if (matchingReads) {
            System.out.println("Will output matching reads ");
        }
        if (nonMatchingReads) {
            System.out.println("Will output non matching reads ");
        }
        if (ambiguousReads) {
            System.out.println("Will output ambiguous reads ");
        }
        if (nonAmbiguousReads) {
            System.out.println("Will output non ambiguous reads ");
        }
        System.out.println("With ambiguity threshold k=" + k);

        return this;
    }

    /**
     * Execute this mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        for (final String basename : basenames) {
            if (suffix == null) {
                suffix = (matchingReads ? "matching" : nonMatchingReads ? "nonmatching" : "")
                        + (ambiguousReads ? "-ambiguous" : "-nonambigious")
                        + ("-k=" + Integer.toString(k));
            }
            if (preFilter != null) {
                System.out.println("Pre-filtering alignment is activated (--pre-filter option).");
            }
            System.out.println("Scanning " + basename);
            alignmentToReadSet(basename);
        }
    }

    private void alignmentToReadSet(final String basename) throws IOException {
        final AlignmentReader reader = new AlignmentReader(basename);
        reader.readHeader();
        final ReadSet outputSet = new ReadSet();
        final int numQueries = reader.getNumberOfQueries();
        final IntSet matchingIndices = new IntOpenHashSet();
        for (final Alignments.AlignmentEntry entry : reader) {
            final int queryIndex = entry.getQueryIndex();

            final int queryId = 1263;
            if (queryIndex == queryId) {
                System.out.println("found 1263");
            }
            matchingIndices.add(queryIndex);

        }

        reader.close();
        final AlignmentTooManyHitsReader tmhReader = new AlignmentTooManyHitsReader(basename);
        for (int queryIndex = 0; queryIndex < numQueries; ++queryIndex) {
            final int queryId = 1263;
            if (queryIndex == queryId) {
                System.out.println("found 1263");
            }
            if (matchingReads) {
                if (matchingIndices.contains(queryIndex)) {
                    if (passesTmhFilter(tmhReader, queryIndex)) {
                        if (preFilter == null || preFilter.contains(queryIndex)) {

                            outputSet.add(queryIndex, 1);
                        }
                        if (queryIndex == queryId) {
                            System.out.println("found 1263 in filter as well, matching read.");
                        }
                    }
                }
            }
            if (nonMatchingReads) {
                if (!matchingIndices.contains(queryIndex)) {
                    if (passesTmhFilter(tmhReader, queryIndex)) {
                        if (preFilter == null || preFilter.contains(queryIndex)) {

                            outputSet.add(queryIndex, 1);
                        }
                        if (queryIndex == queryId) {
                            System.out.println("found 1263 in filter as well, non matching read.");
                        }
                    }
                }
            }
        }

        outputSet.save(basename, suffix);
        System.out.printf("Wrote filter with %d elements", outputSet.size());

    }

    private boolean passesTmhFilter(final AlignmentTooManyHitsReader tmhReader, final int queryIndex) {
        if (ambiguousReads && nonAmbiguousReads) {
            return true;
        } else if (ambiguousReads) {
            return tmhReader.isQueryAmbiguous(queryIndex, k);
        }

        if (nonAmbiguousReads) {
            return !tmhReader.isQueryAmbiguous(queryIndex, k);
        } else {
            return false;
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
        new CompactAlignmentToReadSetMode().configure(args).execute();
    }
}
