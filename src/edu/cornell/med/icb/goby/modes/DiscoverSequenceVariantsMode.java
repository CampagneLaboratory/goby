/*
 * Copyright (C) 2010 Institute for Computational Biomedicine,
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

import java.io.*;
import java.util.Map;
import java.util.Collections;
import java.util.Comparator;

import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.FisherExactRCalculator;
import edu.cornell.med.icb.goby.algorithmic.algorithm.SequenceVariationPool;
import edu.cornell.med.icb.goby.R.GobyRengine;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.lang.MutableString;
import org.rosuda.JRI.Rengine;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 * This mode discovers sequence variants within groups of samples or between groups of samples.
 * Within mode discovery is useful to identify sequence variations that cannot be explained by
 * sequencing error in a given sample. Between mode discovery identify those variants that are found
 * more often in one group versus another.
 *
 * @author Fabien Campagne
 *         Date: Aug 30, 2010
 *         Time: 12:04:59 PM
 */
public class DiscoverSequenceVariantsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "discover-sequence-variants";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Discover sequence variants within and across groups of samples. Identify variations significantly enriched" +
                    "in one group or the other. This mode will either (i) identify sequence variants within a group of sample\n" +
                    "  or (ii) identify variants whose frequency is significantly enriched in one of two groups. \n" +
                    "  This mode requires sorted/indexed alignments as input. (Since Goby 1.8) ";

    private static final Logger LOG = Logger.getLogger(DiscoverSequenceVariantsMode.class);
    private String[] inputFilenames;
    private String outputFile;
    private int[] readerIndexToGroupIndex;
    private int numberOfGroups;
    private MutableString currentReferenceId;
    private int thresholdDistinctReadIndices = 3;
    private int minimumVariationSupport = 10;
    private PrintWriter outWriter;
    private boolean fisherRInstalled;
    private int currentReferenceIndex;
    private String[] groups;
    /**
     * The maximum value of read index.
     */
    private int numberOfReadIndices = Integer.MIN_VALUE;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;
    }

    private final DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();
    private final DifferentialExpressionAnalysis deAnalyzer = new DifferentialExpressionAnalysis();


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
        inputFilenames = jsapResult.getStringArray("input");

        outputFile = jsapResult.getString("output");
        outWriter = "-".equals(outputFile) ? new PrintWriter(System.out) : new PrintWriter(outputFile);

        final String groupsDefinition = jsapResult.getString("groups");
        deAnalyzer.parseGroupsDefinition(groupsDefinition, deCalculator, inputFilenames);
        final String compare = jsapResult.getString("compare");

        deAnalyzer.parseCompare(compare);
        Map<String, String> sampleToGroupMap = deCalculator.getSampleToGroupMap();
        readerIndexToGroupIndex = new int[inputFilenames.length];

        IndexedIdentifier groupIds = new IndexedIdentifier();
        for (String group : sampleToGroupMap.values()) {
            groupIds.registerIdentifier(new MutableString(group));
        }
        minimumVariationSupport = jsapResult.getInt("minimum-variation-support");
        thresholdDistinctReadIndices = jsapResult.getInt("threshold-distinct-read-indices");
        CompactAlignmentToAnnotationCountsMode.parseEval(jsapResult, deAnalyzer);

        numberOfGroups = deAnalyzer.getGroups().length;
        groups = deAnalyzer.getGroups();
        variationPool = new SequenceVariationPool(numberOfGroups);
        for (String sample : sampleToGroupMap.keySet()) {
            final String group = sampleToGroupMap.get(sample);
            System.out.printf("sample: %s group %s%n", sample, group);
            for (int readerIndex = 0; readerIndex < inputFilenames.length; readerIndex++) {
                if (AlignmentReader.getBasename(inputFilenames[readerIndex]).endsWith(sample)) {
                    readerIndexToGroupIndex[readerIndex] = groupIds.get(new MutableString(group));

                }
            }
        }
        File statFile = jsapResult.getFile("variation-stats");
        if (statFile != null) {
            loadStatFile(statFile);
        } else {
            if (deAnalyzer.eval("within-groups")) {
                System.err.println("To evaluate statistics within-groups you must provide a --variation-stats argument.");
                System.exit(1);
            }

        }

        refCount = new int[numberOfGroups];
        variantsCount = new int[numberOfGroups];
        distinctReadIndexCount = new int[numberOfGroups];
        if (deAnalyzer.eval("within-groups") || deAnalyzer.eval("between-groups")) {
            //activate R only if we need it:
            final Rengine rEngine = GobyRengine.getInstance().getRengine();
            fisherRInstalled = rEngine != null && rEngine.isAlive();
        }

        return this;
    }

    class ReadIndexStats {
        public String basename;
        /**
         * The index of the alignment reader that is reading over this basename, will be populated when we know.
         */
        public int readerIndex;
        /**
         * Indexed by readIndex
         */
        public int[] countVariationBases;
        /**
         * Indexed by readIndex
         */
        public int[] countReferenceBases;
    }

    ObjectArrayList<ReadIndexStats> readIndexStats;

    private void loadStatFile(File statFile) {
        try {
            TSVReader reader = new TSVReader(new FileReader(statFile), '\t');
            reader.setCommentPrefix("basename");
            readIndexStats = new ObjectArrayList<ReadIndexStats>();
            ReadIndexStats stat = new ReadIndexStats();
            String lastBasename = null;
            IntArrayList countVariationBases = new IntArrayList();
            IntArrayList countReferenceBases = new IntArrayList();

            String basename = null;
            while (reader.hasNext()) {
                if (reader.isCommentLine() || reader.isEmptyLine()) {
                    // Do nothing, this is a comment or empty line
                    reader.skip();
                } else {
                    reader.next();
                    //if (reader.isCommentLine()) continue;
                    basename = reader.getString();

                    if (lastBasename != null && !lastBasename.equals(basename)) {
                        //we are now processing a new basename. Save the previous stat and start a new one.
                        stat.basename = lastBasename;
                        stat.countVariationBases = countVariationBases.toIntArray();
                        stat.countReferenceBases = countReferenceBases.toIntArray();
                        readIndexStats.add(stat);

                        stat = new ReadIndexStats();
                        countVariationBases.clear();
                        countReferenceBases.clear();
                        lastBasename = basename;

                    }
                    int readIndex = reader.getInt();
                    countVariationBases.add(reader.getInt());

                    assert readIndex == countVariationBases.size();
                    reader.getFloat(); // ignore
                    reader.getFloat(); // ignore
                    countReferenceBases.add(reader.getInt());
                    lastBasename = basename;
                }
            }
            stat.basename = basename;
            stat.countVariationBases = countVariationBases.toIntArray();
            stat.countReferenceBases = countReferenceBases.toIntArray();
            readIndexStats.add(stat);

        } catch (FileNotFoundException e) {
            throw new InternalError("This should never happen.");
        }

        catch (IOException e) {
            System.err.println("Cannot parse stats file. Details may be provided below." +
                    " The file should have been produced with --mode sequence-variation-stats");
            e.printStackTrace(System.err);
            System.exit(1);
        }

        for (ReadIndexStats stat : readIndexStats) {
            numberOfReadIndices = Math.max(numberOfReadIndices, stat.countReferenceBases.length);
        }
    }

    /**
     * Perform the concatenation.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {
        final String outputFilename = outputFile;

        final String[] basenames = AlignmentReader.getBasenames(inputFilenames);
        final boolean allSorted = ConcatenateAlignmentMode.isAllSorted(basenames);
        if (!allSorted) {
            System.out.println("Each input alignment must be sorted. Aborting.");
            System.exit(10);
        }
        if (readIndexStats != null) {
            // associate reader index to basename in the stats, then sort by readerIndex:
            int readerIndex = 0;
            for (String basename : basenames) {
                boolean found = false;
                for (ReadIndexStats stat : readIndexStats) {
                    if (FilenameUtils.getBaseName(basename).equals(stat.basename)) {
                        stat.readerIndex = readerIndex;
                        found = true;
                    }
                }
                if (!found) {
                    System.err.printf("Cannot find basename %s in stat file.", basename);
                }
                readerIndex++;
            }
            Collections.sort(readIndexStats, new Comparator<ReadIndexStats>() {
                public int compare(ReadIndexStats readIndexStats, ReadIndexStats readIndexStatsFirst) {
                    return readIndexStats.readerIndex - readIndexStatsFirst.readerIndex;
                }
            });
            // TODO incorrect, but useful for now. We assume the count of reference bases is uniform over the read indices.
            // TODO remove after mode sequence-variation-stats has been updated to estimate these quantities per read index.
            for (ReadIndexStats stat : readIndexStats) {
                for (int readIndex = 0; readIndex < numberOfReadIndices; readIndex++) {
                    stat.countReferenceBases[readIndex] /= numberOfReadIndices;
                }
            }

        }
        ConcatSortedAlignmentReader sortedReaders = new ConcatSortedAlignmentReader(basenames);
        int lastPosition = -1;
        int lastTargetIndex = -1;
        ObjectArrayList<Alignments.AlignmentEntry> entriesAtPosition = new ObjectArrayList<Alignments.AlignmentEntry>();
        IntArrayList readerIndices = new IntArrayList();

        DoubleIndexedIdentifier targetIds = new DoubleIndexedIdentifier(sortedReaders.getTargetIdentifiers());
        //sortedReaders.skipTo(0,209849405);
        outWriter.printf("referenceId\tposition\t");
        if (deAnalyzer.eval("between-groups")) {
            outWriter.printf("Odds-ratio[%s/%s]\tFisher-Exact-P-value[%s/%s]", groups[0], groups[1],
                    groups[0], groups[1]);
        }


        for (String group : groups) {
            outWriter.printf("\trefCount[%s]\tvarCount[%s]\tdistinct-read-index-count[%s]\taverage-variant-quality-scores[%s]",
                    group, group, group, group);
            
            if (deAnalyzer.eval("within-groups")) {
                outWriter.printf("\twithin-group-p-value[%s]", group);
            }

        }

        outWriter.printf("\tobserved variations at position ([frequency:from/to,]+)%n");


        for (Alignments.AlignmentEntry entry : sortedReaders) {
            if (currentReferenceId == null) {
                currentReferenceId = targetIds.getId(entry.getTargetIndex());
                currentReferenceIndex = entry.getTargetIndex();
            }
            final int position = entry.getPosition();
            final int targetIndex = entry.getTargetIndex();

            if (targetIndex != lastTargetIndex || position != lastPosition) {

                if (entriesAtPosition.size() != 0) {
                    pushVariations(entriesAtPosition, readerIndices);
                }
                processVariations(lastPosition);
                if (targetIndex != lastTargetIndex) {
                    // new reference: clear completely the previous sequence variations.
                    variationPool.reset();
                    currentReferenceId = targetIds.getId(entry.getTargetIndex());
                }
                // clear the variations at the previous position:
                variationPool.removePosition(lastPosition);

                // remove all the entries and readIndices that we have already scanned:
                entriesAtPosition.clear();
                readerIndices.clear();

                lastPosition = position;
                lastTargetIndex = targetIndex;

            }
            entriesAtPosition.add(entry);
            readerIndices.add(sortedReaders.getReaderIndex());

        }
        outWriter.flush();

    }

    int[] refCount;
    int[] variantsCount;
    int[] distinctReadIndexCount;

    private void processVariations(int position) {
        //     System.out.printf("Processing position %d %d %n", lastTargetIndex, position);
        ObjectArrayList<SequenceVariationPool.Variation> list = variationPool.getAtPosition(position);
        int sumVariantCounts = 0;
        for (int i = 0; i < numberOfGroups; i++) {
            refCount[i] = 0;
            variantsCount[i] = 0;
            distinctReadIndexCount[i] = 0;
        }
        if (list != null) {

            int sumDistintReadIndices = 0;
            for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
                refCount[groupIndex] += variationPool.getReferenceAlleleCount(position, groupIndex);

            }
            for (SequenceVariationPool.Variation var : list) {
                variantsCount[var.groupIndex] += 1;
                sumVariantCounts += 1;
            }
            for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {
                distinctReadIndexCount[groupIndex] = variationPool.getNumDistinctReadIndices(position, groupIndex);

                sumDistintReadIndices += distinctReadIndexCount[groupIndex];
            }
            if (sumDistintReadIndices >= thresholdDistinctReadIndices && sumVariantCounts > minimumVariationSupport) {
                int groupIndexA = 0;
                int groupIndexB = 1;


                outWriter.printf("%s\t%d\t", currentReferenceId, position);
                if (deAnalyzer.eval("between-groups")) {
                    final double denominator = (double) refCount[groupIndexA] * (double) variantsCount[groupIndexB];
                    double oddsRatio = denominator == 0 ? 1 :
                            ((double) refCount[groupIndexB] *
                                    (double) variantsCount[groupIndexA]) / denominator;
                    double fisherP = Double.NaN;

                    boolean ok = checkCounts();
                    if (ok) {
                        fisherP = fisherRInstalled ? FisherExactRCalculator.getFisherPValue(
                                refCount[groupIndexB], variantsCount[groupIndexB],
                                refCount[groupIndexA], variantsCount[groupIndexA]) : Double.NaN;
                    } else {
                        System.err.printf("An exception was caught evaluating the Fisher Exact test P-value. Details are provided below%n" +
                                "referenceId=%s referenceIndex=%d position=%d %n" +
                                "refCount[1]=%d variantsCount[1]=%d%n" +
                                "refCount[0]=%d, variantsCount[0]=%d",
                                currentReferenceId, currentReferenceIndex,
                                position,
                                refCount[groupIndexB], variantsCount[groupIndexB],
                                refCount[groupIndexA], variantsCount[groupIndexA]
                        );

                    }
                    outWriter.printf("%f\t%f\t", oddsRatio, fisherP);
                }

                for (int groupIndex = 0; groupIndex < numberOfGroups; groupIndex++) {


                    outWriter.printf("%d\t%d\t%d\t%g\t",
                            refCount[groupIndex],
                            variantsCount[groupIndex],
                            distinctReadIndexCount[groupIndex],
                            variationPool.getAverageVariantQualityScore(position, groupIndex)
                    );

                    if (deAnalyzer.eval("within-groups")) {
                        outWriter.printf("%g\t", estimateWithinGroupDiscoveryPalue(position, groupIndex, list, variantsCount));
                    }
                }

                summarizeVariations(outWriter, list);

                outWriter.println();
            }


        }

    }

    private double estimateWithinGroupDiscoveryPalue(int position, int groupIndex,
                                                     ObjectArrayList<SequenceVariationPool.Variation> list,
                                                     int[] variantsCount) {
        double pValue = 1;
        if (!deAnalyzer.eval("within-groups")) return Double.NaN;
        ObjectArrayList<SequenceVariationPool.ReadIndexInfo> readIndicesAtPosition = variationPool.getReadIndicesAt(position, groupIndex);
        int observedReferenceCount = 0;
        double expectedVariationRate = 0;
        int expectedVariationCount = 0;
        int expectedReferenceCount = 0;
        int observedVariationCount = 0;
        observedReferenceCount = variationPool.getReferenceAlleleCount(position, groupIndex);
        observedVariationCount = variantsCount[groupIndex];
        if (observedReferenceCount + observedVariationCount == 0) {
            // cannot call a variant if we do not observe this position in this group at all.
            return 1;
        }
        long sum = 0;
        if (readIndicesAtPosition == null) return 1;
        for (SequenceVariationPool.ReadIndexInfo readIndexInfo : readIndicesAtPosition) {
            final DiscoverSequenceVariantsMode.ReadIndexStats stats = readIndexStats.get(readIndexInfo.alignmnentReaderIndex);
            final int readIndex = readIndexInfo.readIndex;
            if (readIndex < 0 || readIndex >= stats.countVariationBases.length || readIndex >= stats.countReferenceBases.length) {
                // TODO this test is a kludge: readIndex should always be within the bounds of countVariationBases or countReferenceBases
                // the fact that we need this test indicates that there is a bug in the calculation of readIndex, probably when read insertions or deletions are present.
                continue;
            }
            int variationBases = stats.countVariationBases[readIndex];
            int referenceBases = stats.countReferenceBases[readIndex];
            //  System.out.printf("readIndex=%d variationBases=%d referenceBases=%d %n",readIndexInfo.readIndex, variationBases, referenceBases);
            expectedVariationRate += variationBases;
            sum += variationBases + referenceBases;
            //   refBaseCountAnyPosition += referenceBases;
        }
        expectedVariationRate /= sum;
        //   System.out.printf("Expected variation rate: %f%n", expectedVariationRate);
        expectedVariationCount = (int) Math.round(expectedVariationRate * (double) (observedVariationCount + observedReferenceCount));
        expectedReferenceCount = (int) Math.round((1 - expectedVariationRate) * (double) (observedVariationCount + observedReferenceCount));
        if (LOG.isTraceEnabled()) {
            LOG.trace(String.format("contingency position=%d: %n" +
                    "    [  exp    obs ]%n" +
                    "ref [ %d       %d ]%n" +
                    "var [ %d       %d ] %n",
                    position,
                    expectedVariationCount, observedVariationCount,
                    expectedReferenceCount, observedReferenceCount));
        }
        pValue = fisherRInstalled ? FisherExactRCalculator.getFisherOneTailedLesserPValue(
                expectedVariationCount, observedVariationCount,
                expectedReferenceCount, observedReferenceCount
        ) : Double.NaN;
        //  System.out.printf("position=%d P-Value=%f%n", position, pValue);
        return pValue;


    }

    private boolean checkCounts() {
        boolean ok = true;
        // detect if any count is negative (that's a bug)
        for (int count : refCount) {

            if (count < 0) ok = false;
        }
        for (int count : variantsCount) {
            if (count < 0) ok = false;
        }
        return ok;
    }

    private void summarizeVariations(PrintWriter outWriter, ObjectArrayList<SequenceVariationPool.Variation> list) {
        final Object2IntMap<MutableString> tally = new Object2IntArrayMap<MutableString>();
        tally.defaultReturnValue(0);
        for (SequenceVariationPool.Variation var : list) {
            MutableString variation = new MutableString(var.from);
            variation.append('/');
            variation.append(var.to);
            int count = tally.getInt(variation);
            tally.put(variation, count + 1);

        }
        // sort variations by decreasing tally:
        ObjectSet<MutableString> varSet = (tally.keySet());
        ObjectList<MutableString> varList = new ObjectArrayList<MutableString>();
        varList.addAll(varSet);
        Collections.sort(varList, new Comparator<MutableString>() {
            public int compare(MutableString variationA, MutableString variationB) {
                // sort by decreasing tally:
                return tally.getInt(variationB) - tally.getInt(variationA);
            }
        });

        for (MutableString variation : varList) {
            outWriter.printf("%d:%s,", tally.getInt(variation), variation);
        }
    }

    SequenceVariationPool variationPool;

    private void pushVariations(ObjectArrayList<Alignments.AlignmentEntry> entriesAtPosition,
                                IntArrayList readerIndices) {

        int alignmentReaderIndex = 0;
        // organize all variations found in these entries by their  position on the reference sequence:
        for (Alignments.AlignmentEntry entry : entriesAtPosition) {
            final int groupIndex = readerIndexToGroupIndex[readerIndices.get(alignmentReaderIndex)];

            variationPool.store(entry, groupIndex, alignmentReaderIndex);

        }
        alignmentReaderIndex += 1;
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
        new DiscoverSequenceVariantsMode().configure(args).execute();
    }

}
