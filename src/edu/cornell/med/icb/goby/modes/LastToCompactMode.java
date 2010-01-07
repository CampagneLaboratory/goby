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
import edu.cornell.med.icb.goby.alignments.AlignedSequence;
import edu.cornell.med.icb.goby.alignments.AlignmentStats;
import edu.cornell.med.icb.goby.alignments.AlignmentTooManyHitsWriter;
import edu.cornell.med.icb.goby.alignments.AlignmentWriter;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.LastParser;
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.iterators.TsvLineIterator;
import it.unimi.dsi.fastutil.ints.Int2FloatOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

/**
 * Convert a Last MAF file to the alignment compact format.
 *
 * Can also convert the associated Last COUNTS file (one line per query and count).
 *
 * Default behavior is to convert both maf and count files.  The configure() method
 * resets to the default behavior before parsing command-line options.
 *
 * The input file name can be the full name of the COUNTS file, the full name
 * of the MAF file, or the shared basename of these two files.
 *
 * If MAF and COUNTS files have different basenames, then this tool must be used
 * twice, first converting the maf file using the onlyMaf setting, and second
 * converting the counts file using onlyCounts.
 *
 * This tool assumes that the extensions are ".maf" and ".counts".
 *
 *
 * @author Fabien Campagne
 */
public class LastToCompactMode extends AbstractAlignmentToCompactMode {
    /**
     * The mode name and description
     */
    public static final String MODE_NAME = "last-to-compact";
    public static final String MODE_DESCRIPTION = "Convert a Last MAF file to the alignment compact format.  Also converts the associated Last COUNTS file (one line per query and count) to the alignment too-many-hits format.";

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(LastToCompactMode.class);

    /**
     * default behavior is to convert both maf and count files
     */
    protected boolean onlyMafFile = false;
    protected boolean onlyCountsFile = false;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * If onlyMafFile is set to true, then onlyCountsFile is set to false
     */
    public void setOnlyMafFile(final boolean onlyMafFile) {
        this.onlyMafFile = onlyMafFile;
        if (onlyMafFile) {
            this.onlyCountsFile = false;
        }
    }

    /**
     * If onlyCountsFile is set to true, then onlyMafFile is set to false
     */
    public void setOnlyCountsFile(final boolean onlyCountsFile) {
        this.onlyCountsFile = onlyCountsFile;
        if (onlyCountsFile) {
            this.onlyMafFile = false;
        }
    }



    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        // reset
        setOnlyMafFile(false);
        setOnlyCountsFile(false);

        // configure baseclass
        super.configure(args);

        final JSAPResult jsapResult = parseJsapArguments(args);

        onlyMafFile = jsapResult.getBoolean("only-maf");
        onlyCountsFile = jsapResult.getBoolean("only-counts");
        // default behavior is to convert both maf and count files
        if (onlyMafFile && onlyCountsFile) {
            System.err.println("Only one of --only-maf or --only-counts can be specified.");
            System.exit(1);
        }

        return this;
    }

    @Override
    protected int scan(final ReadSet readIndexFilter, final IndexedIdentifier targetIds, final AlignmentWriter writer,
                     final AlignmentTooManyHitsWriter tmhWriter) throws IOException {

        int numAligns = 0;

        // remove extension from inputFile
        if (FilenameUtils.isExtension(inputFile, new String[]{"maf", "counts"})) {
            inputFile = FilenameUtils.removeExtension(inputFile);
        }
        final String mafInputFile = inputFile + ".maf";
        final String countsInputFile = inputFile + ".counts";

        // convert MAF to compact alignment
        if (!onlyCountsFile) {

            // initialize minimum score & num hits maps
            final Int2FloatOpenHashMap queryIndexToMaxAlignmentScore = new Int2FloatOpenHashMap();
            final Int2IntOpenHashMap queryIndexToNumHitsAtMaxScore = new Int2IntOpenHashMap();

            final ProgressLogger progress = new ProgressLogger(LOG);
            final AlignmentStats stats = new AlignmentStats();
            final int[] readLengths = new int[numberOfReads];

            // first pass: collect minimum score to keep each queryEntry
            // second pass: write to compact alignment file for those entries with score above threshold
            for (final boolean writeAlignment : new boolean[]{false, true}) {
                //
                assert new File(mafInputFile).exists() : "Missing MAF file: "+mafInputFile;
                final LastParser parser = new LastParser(new FileReader(mafInputFile));

                //
                progress.start();
                numAligns = 0;

                while (parser.hasNext()) {

                    // parse maf alignment entry
                    parser.next();
                    final float score = parser.getScore();
                    final ObjectArrayList<AlignedSequence> alignedSequences = parser.getAlignedSequences();

                    // retrieve alignment properties
                    final AlignedSequence reference = alignedSequences.get(0);
                    final AlignedSequence query = alignedSequences.get(1);

                    final int queryIndex = Integer.parseInt(query.sequenceIdentifier.toString());
                    int targetIndex = -1;
                    try {
                        targetIndex = Integer.parseInt(reference.sequenceIdentifier.toString());
                    } catch (NumberFormatException e) {
                        if (targetIds != null) {
                            // BUG?! next line was missing "targetIndex = "
                            targetIndex = targetIds.get(reference.sequenceIdentifier);
                        } else {
                            System.out.println("Cannot convert reference identifier to index. " + reference.sequenceIdentifier);
                            System.exit(1);
                        }
                    }
                    final boolean reverseStrand = !(query.strand == reference.strand);
                    final int depth = query.sequenceLength;
                    final int targetPosition = reference.alignedStart;
                    // save length
                    readLengths[queryIndex] = depth;

                    // we have a multiplicity filter. Use it to determine multiplicity.
                    int multiplicity = 1;
                    if (readIndexFilter != null) {
                        /* Multiplicity of a read is the number of times the (exact) sequence
                         * of the read is identically repeated across a sample file.  The filter
                         * removes duplicates to avoid repeating the same alignments.  Once
                         * aligned, these are recorded multiplicity times.
                         */
                        multiplicity = readIndexFilter.getMultiplicity(queryIndex);
                    }

                    evaluateStatistics(reference, query, stats);

                    final Alignments.AlignmentEntry.Builder currentEntry = Alignments.AlignmentEntry.newBuilder();
                    currentEntry.setNumberOfIndels(stats.numberOfIndels);
                    currentEntry.setNumberOfMismatches(stats.numberOfMismatches);
                    currentEntry.setMatchingReverseStrand(reverseStrand);
                    currentEntry.setMultiplicity(multiplicity);
                    currentEntry.setPosition(targetPosition);
                    currentEntry.setQueryAlignedLength(query.alignedLength);
                    currentEntry.setQueryIndex(queryIndex);
                    currentEntry.setScore(score);
                    currentEntry.setTargetAlignedLength(reference.alignedLength);
                    currentEntry.setTargetIndex(targetIndex);

                    final Alignments.AlignmentEntry alignmentEntry = currentEntry.build();

                    if (qualityFilter.keepEntry(depth, alignmentEntry)) {
                        final float currentMax = queryIndexToMaxAlignmentScore.get(queryIndex);
                        final int currentNumHits = queryIndexToNumHitsAtMaxScore.get(queryIndex);
                        // on the first pass, writeAlignment=false
                        if (writeAlignment == false) {
                            // save the maximum score per read
                            //   and reset the counter to reflect numHits at this new value
                            if (score > currentMax) {
                                queryIndexToMaxAlignmentScore.put(queryIndex, score);
                                queryIndexToNumHitsAtMaxScore.put(queryIndex, 1);
                            }
                            // if query score equals the current max, add 1 to the counter
                            if (score == currentMax) {
                                queryIndexToNumHitsAtMaxScore.put(queryIndex, currentNumHits + 1);
                            }
                        } else {
                            // on the second pass, writeAlignment=true
                            // write the maximum scoring entry (or entries) per read
                            if (score == currentMax) {
                                // only write non-ambiguous entries i.e. currentNumHits <= mParameter
                                if (currentNumHits <= mParameter) {
                                    writer.appendEntry(alignmentEntry);
                                    numAligns += multiplicity;
                                }
                                // TMH writer adds the alignment entry only if hits > thresh
                            }
                        }
                    }

                    progress.lightUpdate();
                }
                parser.close();
                if (writeAlignment) {
                    if (readIndexFilter != null) {
                        writer.putStatistic("keep-filter-filename", readIndexFilterFile.getName());
                    }
                    writer.putStatistic("number-of-entries-written", numAligns);
                    writer.printStats(System.out);
                    /*
                      TODO discuss whether this is needed ...
                      ... it is a lot of additional memory if files are large or indices have not been compacted
                     */
                    // writer.setQueryLengths(readLengths);
                }
                progress.stop();
            }
        }


        // convert COUNTS to compact alignment
        if (!onlyMafFile) {
            assert new File(countsInputFile).exists() : "Missing COUNTS file: "+countsInputFile;

            System.out.println("Recording ambiguity-threshold="+mParameter);
            System.out.println("Will import length of match.");

            for (final Map<String, String> line : new TsvLineIterator(countsInputFile)) {
                final int queryNameToIndex= Integer.parseInt(line.get("query-name"));
                final int depth = Integer.parseInt(line.get("depth"));
                final int count = Integer.parseInt(line.get("number-of-matches"));
                // TMH writer adds the alignment entry only if hits > thresh
                tmhWriter.append(queryNameToIndex, count, depth);
            }
        }

        return numAligns;
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
        new LastToCompactMode().configure(args).execute();
    }

}
