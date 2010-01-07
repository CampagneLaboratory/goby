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
import edu.cornell.med.icb.goby.alignments.Merge;
import edu.cornell.med.icb.util.VersionUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Merge several compact alignment files to an alignment file.   Merge is used
 * when assembling results after searching by chromosome or by transcripts.  It
 * is especially useful if the reference does not fit into memory of a computer
 * used for alignment, for example, when searching against transcripts.
 *
 * Each input MAQ file must have been generated by searching the same set of
 * reads against different subset of reference sequences. Merging consists of
 * putting back the results as if the set of reads had been searched against
 * the combined reference.
 *
 * Several merging strategies are supported by this class. For instance, one
 * strategy assumes that each set of reads provided as input was searched
 * against a different chromosome (or contig). Another strategy assumes that
 * the reads were aligned to cDNA/transcript individual reference sequences.
 *
 * @author Kevin Dorff
 * @author Fabien Campagne
 */
public class MergeCompactAlignmentsMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    public static final String MODE_NAME = "merge-compact-alignments";
    public static final String MODE_DESCRIPTION = "Merge several compact alignment files to an alignment file.   Merge is used when assembling results after searching by chromosome or by transcripts.  It is especially useful if the reference does not fit into memory of a computer used for alignment, for example, when searching against transcripts.  Each input MAQ file must have been generated by searching the same set of reads against different subset of reference sequences. Merging consists of putting back the results as if the set of reads had been searched against the combined reference.  Several merging strategies are supported by this class. For instance, one strategy assumes that each set of reads provided as input was searched against a different chromosome (or contig). Another strategy assumes that the reads were aligned to cDNA/transcript individual reference sequences.";


    /**
     * Max number of duplicates allowed for the top quality score when there are dupes.
     */
    private static final int K_NUM_OF_BEST_QUAL_TO_KEEP = 2;

    /**
     * The input file.
     */
    private List<File> inputFiles;

    private String geneTranscriptMapFile;

    /**
     * The output file.
     */
    private String outputFile;

    /**
     * Max number of duplicates allowed for the top quality score when there are dupes.
     */
    private int k;

    /**
     * Map to override help / default values.
     */
    private static final Map<String, String> HELP_VALUES;

    static {
        HELP_VALUES = new HashMap<String, String>();
        HELP_VALUES.put("[K_VALUE]", Integer.toString(K_NUM_OF_BEST_QUAL_TO_KEEP));
    }

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * Create the MaqToolsDriver object, parse argument(s). This intentionaly doesn't check
     * for errors as we only care for "mode" and we KNOW there will be other command line
     * parameters. All other arguments are specifed by the specific mode being called
     * (including things like --help).
     *
     * @param args the arguments
     * @return this object for chaining
     * @throws java.io.IOException error parsing arguments
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing arguments
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {

        final JSAPResult jsapResult = parseJsapArguments(args, HELP_VALUES);

        final String[] inputFilesAr = jsapResult.getStringArray("input");

        inputFiles = new ArrayList<File>(inputFilesAr.length);
        for (final String file : inputFilesAr) {
            final String basename = AlignmentReader.getBasename(file);
            inputFiles.add(new File(basename));
        }
        outputFile = jsapResult.getString("output");
        k = jsapResult.getInt("k");
        geneTranscriptMapFile = jsapResult.getString("gene-transcript-map-file");
        System.out.println("Configured with k=" + k);
        return this;
    }

    /**
     * Run the merge mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        System.out.println("Version: "+ VersionUtils.getImplementationVersion(MergeCompactAlignmentsMode.class));
        final Merge merger = new Merge(geneTranscriptMapFile, k);
        merger.setSilent(false);
        merger.merge(inputFiles, outputFile);
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
        new MergeCompactAlignmentsMode().configure(args).execute();
    }
}
