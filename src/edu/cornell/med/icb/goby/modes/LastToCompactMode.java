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
import edu.cornell.med.icb.goby.reads.ReadSet;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.iterators.TsvLineIterator;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Convert a Last MAF file to the alignment compact format.
 * <p/>
 * Can also convert the associated Last COUNTS file (one line per query and count).
 * <p/>
 * Default behavior is to convert both maf and count files.  The configure() method
 * resets to the default behavior before parsing command-line options.
 * <p/>
 * The input file name can be the full name of the COUNTS file, the full name
 * of the MAF file, or the shared basename of these two files.
 * <p/>
 * If MAF and COUNTS files have different basenames, then this tool must be used
 * twice, first converting the maf file using the onlyMaf setting, and second
 * converting the counts file using onlyCounts.
 * <p/>
 * This tool assumes that the extensions are ".maf" and ".counts".
 *
 * @author Fabien Campagne
 */
public class LastToCompactMode extends AbstractAlignmentToCompactMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "last-to-compact";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Convert a Last MAF file to the alignment "
            + "compact format.  Also converts the associated Last COUNTS file (one line per "
            + "query and count) to the alignment too-many-hits format.";

    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(LastToCompactMode.class);

    /**
     * Default behavior is to convert both maf and count files.
     */
    protected boolean onlyMafFile;
    protected boolean onlyCountsFile;
    private boolean flipStrand;
    private boolean substituteCharacter;
    private Substitution[] substitutions;
    private boolean hasPair;
    private boolean firstInPair;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    /**
     * If {@link #onlyMafFile} is set to true, then {@link #onlyCountsFile} is set to false.
     */
    public void setOnlyMafFile(final boolean onlyMafFile) {
        this.onlyMafFile = onlyMafFile;
        if (onlyMafFile) {
            this.onlyCountsFile = false;
        }
    }

    /**
     * If {@link #onlyCountsFile} is set to true, then {@link #onlyMafFile} is set to false.
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
        flipStrand = jsapResult.getBoolean("flip-strand");
        final String substitutionString = jsapResult.getString("substitutions");
        if (substitutionString != null) {
            final String[] tokens = substitutionString.split(",");
            substitutions = new Substitution[tokens.length];
            int i = 0;
            for (final String token : tokens) {
                substitutions[i] = new Substitution();
                final String[] parts = token.split("/");
                if (parts[0].length() != 1 || parts[1].length() != 1) {
                    System.err.println("Substitutions must have length 1. Error parsing " + token);
                    System.exit(1);
                }
                substitutions[i].from = parts[0].charAt(0);
                substitutions[i].to = parts[1].charAt(0);
                i++;
            }
            substituteCharacter = true;
        } else {
            substitutions = new Substitution[0];
            substituteCharacter = false;
        }
        return this;
    }

    @Override
    protected int scan(final ReadSet readIndexFilter, final IndexedIdentifier targetIds,
                       final AlignmentWriter writer,
                       final AlignmentTooManyHitsWriter tmhWriter) throws IOException {

        int currentQueryIndex = -1;
        int currentFragmentIndex = -1;
        final List<Alignments.AlignmentEntry.Builder> sameQueryIndexAlignmentEntries =
                new ArrayList<Alignments.AlignmentEntry.Builder>();

        int numAligns = 0;

        // remove extension from inputFile
        if (FilenameUtils.isExtension(inputFile, new String[]{"maf", "counts"})) {
            inputFile = FilenameUtils.removeExtension(inputFile);
        }
        final String mafInputFile = inputFile + ".maf";
        final String countsInputFile = inputFile + ".counts";

        // convert MAF to compact alignment
        if (!onlyCountsFile) {

            // initialize minimum score & num hits maps. The following arrays are indexed on the fragmentIndex:
            final Int2FloatOpenHashMap queryIndexToMaxAlignmentScore[] = new Int2FloatOpenHashMap[]{new Int2FloatOpenHashMap(), new Int2FloatOpenHashMap()};
            final Int2IntOpenHashMap queryIndexToNumHitsAtMaxScore[] = new Int2IntOpenHashMap[]{new Int2IntOpenHashMap(), new Int2IntOpenHashMap()};


            final AlignmentStats stats = new AlignmentStats();
            //      final int[] readLengths = createReadLengthArray();
            IntArrayList targetLengths = new IntArrayList();
            // first pass: collect minimum score to keep each queryEntry
            // second pass: write to compact alignment file for those entries with score above threshold
            for (final boolean writeAlignment : new boolean[]{false, true}) {
                assert new File(mafInputFile).exists() : "Missing MAF file: " + mafInputFile;
                final LastParser parser = new LastParser(new FileReader(mafInputFile));

                //
                final ProgressLogger progress = new ProgressLogger(LOG);
                progress.start();
                numAligns = 0;
                int removedByQualityFilter = 0;
                int notBestScore = 0;
                while (parser.hasNext()) {
                    // parse maf alignment entry
                    parser.next();
                    final float score = parser.getScore();
                    final ObjectArrayList<AlignedSequence> alignedSequences = parser.getAlignedSequences();

                    // retrieve alignment properties
                    final AlignedSequence reference = alignedSequences.get(0);
                    final AlignedSequence query = alignedSequences.get(1);
                    if (flipStrand) {
                        // override the query strand with forceStrand if requested on the command line.
                        query.strand = query.strand == '+' ? '-' : query.strand == '-' ? '+' : '?';
                        flip(reference.alignment);
                        flip(query.alignment);
                    }
                    if (substituteCharacter) {
                        for (final Substitution sub : substitutions) {
                            query.alignment.replace(sub.from, sub.to);
                        }
                    }
                    String s = query.sequenceIdentifier.toString();
                    int fragmentIndex = 0;
                    if (s.endsWith("/1")) {
                        s = s.substring(0, s.length() - 2);
                        fragmentIndex = 0;
                        hasPair=true;
                        firstInPair=true;
                    }
                    if (s.endsWith("/2")) {
                        s = s.substring(0, s.length() - 2);
                        fragmentIndex = 1;
                        hasPair=true;
                        firstInPair=false;
                    }
                    final int queryIndex = Integer.parseInt(s);
                    if (currentQueryIndex == -1) {
                        currentQueryIndex = queryIndex;
                    }
                    if (currentFragmentIndex == -1) {
                        currentFragmentIndex = fragmentIndex;
                    }
                    largestQueryIndex = Math.max(queryIndex, largestQueryIndex);
                    int targetIndex = -1;
                    targetIndex = getTargetIndex(targetIds, reference.sequenceIdentifier, thirdPartyInput);
                    final boolean reverseStrand = !(query.strand == reference.strand);
                    final int depth = query.sequenceLength;
                    final int targetPosition = reference.alignedStart;

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
                    currentEntry.setFragmentIndex(fragmentIndex);
                    currentEntry.setScore(score);
                    currentEntry.setTargetAlignedLength(reference.alignedLength);

                    storeForPairsIndex(writeAlignment, queryIndex, fragmentIndex, targetIndex, targetPosition, reverseStrand);

                    if (writeAlignment && hasPair) {
                        final Alignments.RelatedAlignmentEntry.Builder builderForValue = Alignments.RelatedAlignmentEntry.newBuilder();
                        builderForValue.setFragmentIndex(flipFragmentIndex(fragmentIndex));
                        builderForValue.setTargetIndex(getPairTargetIndex(queryIndex, fragmentIndex));
                        builderForValue.setPosition(getPairPosition(queryIndex, fragmentIndex));
                        currentEntry.setPairAlignmentLink(builderForValue.build());
                        int pairFlags = EntryFlagHelper.paired();
                        if (firstInPair) pairFlags|= EntryFlagHelper.firstInPair();
                        if (getPairStrand(queryIndex, fragmentIndex) && firstInPair) pairFlags|= EntryFlagHelper.mateReverseStrand();
                        if (getPairStrand(queryIndex, fragmentIndex) && !firstInPair) pairFlags|= EntryFlagHelper.readReverseStrand();
                        if (getPairStrand(queryIndex, flipFragmentIndex(fragmentIndex)) && firstInPair) pairFlags|= EntryFlagHelper.readReverseStrand();
                        if (getPairStrand(queryIndex, flipFragmentIndex(fragmentIndex)) && !firstInPair) pairFlags|= EntryFlagHelper.mateReverseStrand();
                        currentEntry.setPairFlags(pairFlags);
                    }

                    if (targetLengths.size() <= targetIndex) {
                        targetLengths.size(targetIndex + 1);
                    }
                    targetLengths.set(targetIndex, reference.sequenceLength);
                    currentEntry.setTargetIndex(targetIndex);
                    final int queryLength = query.sequenceLength;
                    currentEntry.setQueryLength(queryLength);
                    final int readStartPosition = query.alignedStart;

                    parseSequenceVariations(currentEntry, reference, query, readStartPosition, queryLength, reverseStrand);

                    if (qualityFilter.keepEntry(depth, currentEntry)) {
                        final float currentMax = queryIndexToMaxAlignmentScore[fragmentIndex].get(queryIndex);
                        final int currentNumHits = queryIndexToNumHitsAtMaxScore[fragmentIndex].get(queryIndex);
                        // on the first pass, writeAlignment=false
                        if (!writeAlignment) {
                            // save the maximum score per read
                            //   and reset the counter to reflect numHits at this new value
                            if (score > currentMax) {
                                queryIndexToMaxAlignmentScore[fragmentIndex].put(queryIndex, score);
                                queryIndexToNumHitsAtMaxScore[fragmentIndex].put(queryIndex, 1);
                            }
                            // if query score equals the current max, add 1 to the counter
                            if (score == currentMax) {
                                queryIndexToNumHitsAtMaxScore[fragmentIndex].put(queryIndex, currentNumHits + 1);
                            }
                        } else {
                            // on the second pass, writeAlignment=true
                            // write the maximum scoring entry (or entries) per read
                            if (score == currentMax) {
                                // only write non-ambiguous entries i.e. currentNumHits <= mParameter
                                if (currentNumHits <= mParameter) {

                                    if (currentEntry.getQueryIndex() == currentQueryIndex ) {
                                        sameQueryIndexAlignmentEntries.add(currentEntry);

                                    } else {
                                        writeEntries(writer, sameQueryIndexAlignmentEntries);
                                        sameQueryIndexAlignmentEntries.add(currentEntry);
                                        currentQueryIndex = currentEntry.getQueryIndex();
                                        currentFragmentIndex = currentEntry.getFragmentIndex();
                                        numAligns += multiplicity;
                                    }
                                }
                                // TMH writer adds the alignment entry only if hits > thresh
                            } else {
                                notBestScore++;
                                //        System.out.println("Excluding entry "+alignmentEntry);
                            }
                        }
                    } else {
                        removedByQualityFilter++;
                    }

                    progress.lightUpdate();

                }
                parser.close();
                if (writeAlignment) {
                    // Write the remaining entries (last query index);
                    numAligns += writeEntries(writer, sameQueryIndexAlignmentEntries);
                    writer.setTargetLengths(targetLengths.toIntArray(new int[targetLengths.size()]));
                    if (readIndexFilter != null) {
                        writer.putStatistic("keep-filter-filename", readIndexFilterFile.getName());
                    }
                    writer.putStatistic("number-of-entries-written", numAligns);

                    writer.setNumQueries(numberOfReads);
                    writer.printStats(System.out);
                    System.out.printf("Removed by quality filter: %d%n", removedByQualityFilter);
                    System.out.printf("Not best score: %d%n", notBestScore);
                }
                progress.stop();
                currentQueryIndex = -1;
                currentFragmentIndex = -1;
            }
        }


        // convert COUNTS to compact alignment
        if (!onlyMafFile) {
            assert new File(countsInputFile).exists() : "Missing COUNTS file: " + countsInputFile;

            System.out.println("Recording ambiguity-threshold=" + mParameter);
            System.out.println("Will import length of match.");

            for (final Map<String, String> line : new TsvLineIterator(countsInputFile)) {
                final int queryNameToIndex = Integer.parseInt(line.get("query-name"));
                final int depth = Integer.parseInt(line.get("depth"));
                final int count = Integer.parseInt(line.get("number-of-matches"));
                // TMH writer adds the alignment entry only if hits > thresh
                tmhWriter.append(queryNameToIndex, count, depth);
            }
        }

        return numAligns;
    }

    // swap the fragment index to locate the other read in the pair:
    private int flipFragmentIndex(int fragmentIndex) {
        return fragmentIndex ^ 1;
    }

    private int getPairTargetIndex(int queryIndex, int fragmentIndex) {
        // swap the fragment index to locate the other read in the pair:
        fragmentIndex = flipFragmentIndex(fragmentIndex);
        return targetMap[fragmentIndex].get(queryIndex);
    }

    private int getPairPosition(int queryIndex, int fragmentIndex) {
        // swap the fragment index to locate the other read in the pair:
        fragmentIndex = flipFragmentIndex(fragmentIndex);
        return positionMap[fragmentIndex].get(queryIndex);
    }
    private boolean getPairStrand(int queryIndex, int fragmentIndex) {
        // swap the fragment index to locate the other read in the pair:
        fragmentIndex = flipFragmentIndex(fragmentIndex);
        return strandMap[fragmentIndex].get(queryIndex);
    }

    /**
     * The following arrays positionMap and targetMap are indexed on the fragment index.
     * positionMap stores the position of each queryIndex, fragmentIndex
     * targetMap stores the target of each queryIndex, fragmentIndex
     */
    Int2IntAVLTreeMap positionMap[] = new Int2IntAVLTreeMap[]{new Int2IntAVLTreeMap(), new Int2IntAVLTreeMap()};
    Int2IntAVLTreeMap targetMap[] = new Int2IntAVLTreeMap[]{new Int2IntAVLTreeMap(), new Int2IntAVLTreeMap()};
    Int2BooleanAVLTreeMap strandMap[] = new Int2BooleanAVLTreeMap[]{new Int2BooleanAVLTreeMap(), new Int2BooleanAVLTreeMap()};

    private void storeForPairsIndex(boolean writeAlignment, int queryIndex, int fragmentIndex, int targetIndex, int position, boolean matchesReverseStrand) {
        if (!writeAlignment) {
            positionMap[fragmentIndex].put(queryIndex, position);
            targetMap[fragmentIndex].put(queryIndex, targetIndex);
            strandMap[fragmentIndex].put(queryIndex, matchesReverseStrand);
        }
    }

    private int writeEntries(final AlignmentWriter writer,
                             final List<Alignments.AlignmentEntry.Builder> alignmentEntries) throws IOException {

        int numAligns = 0;
        final int size = alignmentEntries.size();
        if (size > 0) {
            for (final Alignments.AlignmentEntry.Builder alignmentEntry : alignmentEntries) {
                alignmentEntry.setQueryIndexOccurrences(size);
                alignmentEntry.setAmbiguity(1);
                writer.appendEntry(alignmentEntry.build());
                numAligns += alignmentEntry.getMultiplicity();
            }
            alignmentEntries.clear();
        }
        return numAligns;
    }

    private void flip(final MutableString alignment) {

        for (int i = 0; i < alignment.length(); i++) {
            char base = alignment.charAt(i);
            switch (base) {
                case 'A':
                    base = 'T';
                    break;

                case 'C':
                    base = 'G';
                    break;
                case 'T':
                    base = 'A';
                    break;
                case 'G':
                    base = 'C';
                    break;
            }
            alignment.setCharAt(i, base);
        }
    }

    private void parseSequenceVariations(final Alignments.AlignmentEntry.Builder currentEntry,
                                         final AlignedSequence reference,
                                         final AlignedSequence query,
                                         final int readStartPosition,
                                         final int queryLength,
                                         final boolean reverseStrand) {
        final int alignmentLength = reference.alignment.length();
        final MutableString referenceSequence = reference.alignment;
        final MutableString querySequence = query.alignment;

        extractSequenceVariations(currentEntry, alignmentLength, referenceSequence, querySequence, readStartPosition,
                queryLength, reverseStrand, query.getQualityScores().toByteArray());
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException error parsing
     * @throws java.io.IOException                    error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new LastToCompactMode().configure(args).execute();
    }

    private class Substitution {
        // substitute this character
        char from;
        // by this one:
        char to;
    }
}
