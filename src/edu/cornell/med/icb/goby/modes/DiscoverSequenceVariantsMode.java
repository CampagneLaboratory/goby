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
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceCache;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;

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
            "Discover sequence variants within and across groups of samples. This mode will either (i) identify sequence variants within a group of sample\n" +
                    "  or (ii) identify variants whose frequency is significantly enriched in one of two groups. \n" +
                    "  This mode requires sorted/indexed alignments as input. (Since Goby 1.8) ";

    private static final Logger LOG = Logger.getLogger(DiscoverSequenceVariantsMode.class);
    private String[] inputFilenames;
    private String outputFile;
    private int[] readerIndexToGroupIndex;
    private int numberOfGroups;
    private CharSequence currentReferenceId;
    private int thresholdDistinctReadIndices = 3;
    private int minimumVariationSupport = 10;
    private PrintWriter outWriter;

    private String[] groups;
    /**
     * The maximum value of read index, indexed by readerIndex;
     */
    private int numberOfReadIndices[];
    private DifferentialExpressionCalculator diffExpCalculator;
    private String[] samples;
    private boolean outputVCF;

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
        boolean parallel = jsapResult.getBoolean("parallel", false);
        deAnalyzer.setRunInParallel(parallel);
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

        for (String sample : sampleToGroupMap.keySet()) {
            final String group = sampleToGroupMap.get(sample);
            System.out.printf("sample: %s group %s%n", sample, group);
            for (int readerIndex = 0; readerIndex < inputFilenames.length; readerIndex++) {
                if (AlignmentReaderImpl.getBasename(inputFilenames[readerIndex]).endsWith(sample)) {
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
        final String formatSring = jsapResult.getString("format");

        OutputFormat format = OutputFormat.valueOf(formatSring.toUpperCase());

        SequenceVariationOutputFormat formatter = null;
        switch (format) {
            case VARIANT_DISCOVERY:
                formatter = new BetweenGroupSequenceVariationOutputFormat();
                break;
            case COMPARE_GROUPS:
                formatter = new CompareGroupsVCFOutputFormat();
                break;
            case ALLELE_FREQUENCIES:
                formatter = new AlleleFrequencyOutputFormat();
                break;
            case GENOTYPES:
                formatter = new GenotypesOutputFormat();
                break;
            case METHYLATION:
                formatter = new MethylationRateVCFOutputFormat();
                // methylated bases match the reference. Do not filter on minimum variation support.
                this.minimumVariationSupport=-1;
                this.thresholdDistinctReadIndices=1;
                System.out.println("Methylation format ignores thresholdDistinctReadIndices and minimumVariationSupport.");
                break;
            default:
                ObjectArrayList<OutputFormat> values = ObjectArrayList.wrap(OutputFormat.values());
                System.err.printf("The format argument is not recognized. Allowed values include %s",
                        values.toString());
                System.exit(1);
        }

        final String genome = jsapResult.getString("genome");
        RandomAccessSequenceCache cache = null;
        if (genome != null) {
            try {
                System.err.println("Loading genome cache " + genome);
                cache = new RandomAccessSequenceCache();
                cache.load(genome);
            } catch (ClassNotFoundException e) {
                System.err.println("Could not load genome cache");
                e.printStackTrace();
                System.exit(1);
            }
        }

        outputVCF = jsapResult.getBoolean("vcf");
        int startFlapSize = jsapResult.getInt("start-flap-size", 100);

        sortedPositionIterator = new DiscoverVariantIterateSortedAlignments(formatter);

        sortedPositionIterator.setGenome(cache);
        sortedPositionIterator.setStartFlapLength(startFlapSize);
        sortedPositionIterator.parseIncludeReferenceArgument(jsapResult);
        sortedPositionIterator.setMinimumVariationSupport(minimumVariationSupport);
        sortedPositionIterator.setThresholdDistinctReadIndices(thresholdDistinctReadIndices);
        return this;
    }


    enum OutputFormat {
        VARIANT_DISCOVERY,
        ALLELE_FREQUENCIES,
        GENOTYPES,
        COMPARE_GROUPS,
        METHYLATION
    }

    public boolean outputVCF() {
        return outputVCF;
    }

    DiscoverVariantIterateSortedAlignments sortedPositionIterator;

    public DifferentialExpressionAnalysis getDiffExpAnalyzer() {
        return deAnalyzer;
    }

    public DifferentialExpressionCalculator getDiffExpCalculator() {
        return diffExpCalculator;
    }

    public String[] getGroups() {
        return groups;
    }

    public ObjectArrayList<ReadIndexStats> getReadIndexStats() {
        return readIndexStats;
    }

    public String[] getSamples() {
        // build a samples array with the correct order:
        int numberOfSamples = readIndexStats.size();
        samples = new String[numberOfSamples];
        if (readIndexStats.size()==0) {
            System.err.println("Cannot find any basename in stats file. Aborting.");
            System.exit(1);
        }
        for (ReadIndexStats stat : readIndexStats) {
            samples[stat.readerIndex] = stat.basename;
        }
        return samples;
    }

    public int[] getReaderIndexToGroupIndex() {
        return readerIndexToGroupIndex;
    }

    ObjectArrayList<ReadIndexStats> readIndexStats;

    private void loadStatFile(File statFile) {
        try {
            TSVReader reader = new TSVReader(new FileReader(statFile), '\t');
            reader.setCommentPrefix("basename");
            readIndexStats = new ObjectArrayList<ReadIndexStats>();
            ReadIndexStats stat = new ReadIndexStats();
            String lastBasename = null;
            LongArrayList countVariationBases = new LongArrayList();
            LongArrayList countReferenceBases = new LongArrayList();

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
                        stat.countVariationBases = countVariationBases.toLongArray();
                        stat.countReferenceBases = countReferenceBases.toLongArray();
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
                    reader.getLong(); // ignore
                    countReferenceBases.add(reader.getLong());
                    lastBasename = basename;
                }
            }
            stat.basename = basename;
            stat.countVariationBases = countVariationBases.toLongArray();
            stat.countReferenceBases = countReferenceBases.toLongArray();
            readIndexStats.add(stat);

        } catch (FileNotFoundException e) {
            System.err.printf("Error. The -v file argument cannot be found (%s)%n", statFile);
            System.exit(1);
        }

        catch (IOException e) {
            System.err.println("Cannot parse stats file. Details may be provided below." +
                    " The file should have been produced with --mode sequence-variation-stats");
            e.printStackTrace(System.err);
            System.exit(1);
        }
        if (readIndexStats.size() < inputFilenames.length) {
            System.err.printf("The stats file seems incomplete. Expected to find statistics for %d samples, but found only %d %n",
                    inputFilenames.length, readIndexStats.size());
            System.exit(1);
        }
    }

    /**
     * Perform the concatenation.
     *
     * @throws java.io.IOException
     */
    @Override
    public void execute() throws IOException {

        final String[] basenames = AlignmentReaderImpl.getBasenames(inputFilenames);
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
                    System.err.flush();
                }
                readerIndex++;
            }
            Collections.sort(readIndexStats, new Comparator<ReadIndexStats>() {
                public int compare(ReadIndexStats readIndexStats, ReadIndexStats readIndexStatsFirst) {
                    return readIndexStats.readerIndex - readIndexStatsFirst.readerIndex;
                }
            });

            // Determine the maximum read length for each input sample (fill numberOfReadIndices)

            final int sampleToGroupAssociationNumber = this.deCalculator.getSampleToGroupMap().keySet().size();
            // some samples provided on the input may not be associated with groups. Use the number of input filenames
            // to allocate the size of ths array.
            numberOfReadIndices = new int[inputFilenames.length];

            ObjectSet<ReadIndexStats> toRemove = new ObjectArraySet<ReadIndexStats>();

            for (ReadIndexStats stat : readIndexStats) {
                if (stat.readerIndex == -1) {
                    // this sample was not loaded, remove it from consideration
                    toRemove.add(stat);
                    continue;
                }
                numberOfReadIndices[stat.readerIndex] = Math.max(numberOfReadIndices[stat.readerIndex], stat.countReferenceBases.length);
            }

            readIndexStats.removeAll(toRemove);
        }
        sortedPositionIterator.allocateStorage(basenames.length, numberOfGroups);
        sortedPositionIterator.initialize(this, outWriter);
        // install a reader factory that filters out ambiguous reads:
        sortedPositionIterator.setAlignmentReaderFactory(new NonAmbiguousAlignmentReaderFactory());
        sortedPositionIterator.iterate(basenames);
        sortedPositionIterator.finish();
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
        System.exit(0);
    }

}
