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
import edu.cornell.med.icb.goby.Release1_9_7_2;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.alignments.*;
import edu.cornell.med.icb.goby.alignments.processors.*;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceCache;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

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
            "This mode is discovering sequence variants, calling genotypes, estimating allele frequencies, " +
                    "or methylation rates (RRBS or methyl-seq datasets)." +
                    " This mode requires sorted/indexed alignments as input. (Since Goby 1.8) ";

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
    private boolean groupsAreDefined;
    private ObjectArrayList<GenotypeFilter> genotypeFilters;
    private AlignmentProcessorFactory realignmentFactory = new DefaultAlignmentProcessorFactory();
    /**
     * A genome used for testing.
     */
    private RandomAccessSequenceInterface testGenome;
    /**
     * When this flag is true, the mode overrides reference alleles with that of the genome when there is a disagreement
     * between the genome and the allele.
     */
    private boolean overrideReferenceWithGenome = true;
    private FormatConfigurator formatConfigurator = new DummyFormatConfigurator();
    private ArrayList<GroupComparison> groupComparisonsList = new ArrayList<GroupComparison>();
    private int maxThresholdPerSite;


    public void setDisableAtLeastQuarterFilter(boolean disableAtLeastQuarterFilter) {
        this.disableAtLeastQuarterFilter = disableAtLeastQuarterFilter;
    }

    private boolean disableAtLeastQuarterFilter = false;

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

        String groupsDefinition = jsapResult.getString("groups");
        String groupsDefinitionFile = jsapResult.getString("group-file");
        if (groupsDefinition != null && groupsDefinitionFile != null) {
            System.err.println("--groups and --groups-file are mutually exclusive. Please provide only once such parameter.");
            System.exit(1);
        }
        if (groupsDefinitionFile != null) {
            groupsDefinition = parseGroupFile(groupsDefinitionFile, inputFilenames);
        }
        String compare = jsapResult.getString("compare");
        /*if (compare == null) {
            // make default groups and group definitions.  Each sample becomes its own group, compare group1 and group2.

            compare = "group1/group2";
            MutableString buffer = new MutableString();
            int groupIndex = 1;
            for (String inputFilename : inputFilenames) {
                buffer.append("group").append(String.valueOf(groupIndex++));
                buffer.append("=");
                buffer.append(AlignmentReaderImpl.getBasename(inputFilename));
                buffer.append('/');
            }
            groupsDefinition = buffer.substring(0, buffer.length() - 1).toString();
            System.out.println(groupsDefinition);
        } else {
            groupsAreDefined = true;
        }     */
        deAnalyzer.parseGroupsDefinition(groupsDefinition, deCalculator, inputFilenames);
        groupComparisonsList = deAnalyzer.parseCompare(compare);


        boolean parallel = jsapResult.getBoolean("parallel", false);
        deAnalyzer.setRunInParallel(parallel);
        Map<String, String> sampleToGroupMap = deCalculator.getSampleToGroupMap();
        readerIndexToGroupIndex = new int[inputFilenames.length];

        groups = deAnalyzer.getGroups();
        numberOfGroups = groups.length;
        IndexedIdentifier groupIds = new IndexedIdentifier();
        for (String group : groups) {
            groupIds.registerIdentifier(new MutableString(group));
        }
        maxThresholdPerSite = jsapResult.getInt("max-coverage-per-site");
        minimumVariationSupport = jsapResult.getInt("minimum-variation-support");
        thresholdDistinctReadIndices = jsapResult.getInt("threshold-distinct-read-indices");
        CompactAlignmentToAnnotationCountsMode.parseEval(jsapResult, deAnalyzer);

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

        realignmentFactory = configureProcessor(jsapResult);
        final String formatString = jsapResult.getString("format");

        final OutputFormat format = OutputFormat.valueOf(formatString.toUpperCase());

        SequenceVariationOutputFormat formatter = null;
        switch (format) {
            case VARIANT_DISCOVERY:
            case BETWEEN_GROUPS:
                stopWhenDefaultGroupOptions();
                formatter = new BetweenGroupSequenceVariationOutputFormat();
                break;
            case COMPARE_GROUPS:
                stopWhenDefaultGroupOptions();
                formatter = new CompareGroupsVCFOutputFormat();
                break;
            case ALLELE_FREQUENCIES:
                stopWhenDefaultGroupOptions();
                formatter = new AlleleFrequencyOutputFormat();
                break;
            case GENOTYPES:
                formatter = new GenotypesOutputFormat();
                break;
            case METHYLATION:
                formatter = new MethylationRateVCFOutputFormat();
                // methylated bases match the reference. Do not filter on minimum variation support.
                int tmp = minimumVariationSupport;
                this.minimumVariationSupport = -1;
                this.thresholdDistinctReadIndices = 1;
                // need at least so many methylation/non-methylation event to record site in output
                // the value configure put in minimumVariationSupport as minimum coverage for the site.
                ((MethylationRateVCFOutputFormat) formatter).setMinimumEventThreshold(tmp);
                System.out.println("Methylation format ignores thresholdDistinctReadIndices. Additionally, the minimum coverage needed for a site to be reported can be changed with --minimum-variation-support.");
                break;
            case INDEL_COUNTS:
                formatter = new IndelCountOutputFormat();
                break;
            default:
                ObjectArrayList<OutputFormat> values = ObjectArrayList.wrap(OutputFormat.values());
                System.err.printf("The format argument is not recognized. Allowed values include %s",
                        values.toString());
                System.exit(1);
        }

        // set base filters according to output format:
        genotypeFilters = new ObjectArrayList<GenotypeFilter>();
        switch (format) {

            case METHYLATION:
                // no filters at all for methylation. It seems that in some dataset quality scores are much lower for
                //bases that are not methylation and therefore converted. We don't want to filter these bases and therefore
                // do not install the quality score filter.
                if (Release1_9_7_2.callIndels) {
                    genotypeFilters.add(new RemoveIndelArtifactsFilter());
                }
                break;
            case COMPARE_GROUPS:
            case ALLELE_FREQUENCIES:
            case BETWEEN_GROUPS:
            case VARIANT_DISCOVERY:

                genotypeFilters.add(new QualityScoreFilter());
                genotypeFilters.add(new LeftOverFilter());
                if (Release1_9_7_2.callIndels) {
                    genotypeFilters.add(new RemoveIndelArtifactsFilter());
                }
                break;
            case GENOTYPES:

                genotypeFilters.add(new QualityScoreFilter());
                genotypeFilters.add(new LeftOverFilter());

                if (Release1_9_7_2.callIndels) {
                    genotypeFilters.add(new RemoveIndelArtifactsFilter());
                }
                if (!disableAtLeastQuarterFilter) {
                    genotypeFilters.add(new AtLeastAQuarterFilter());
                }
                break;
            case INDEL_COUNTS:
                genotypeFilters.add(new QualityScoreFilter());
                genotypeFilters.add(new LeftOverFilter());
                genotypeFilters.add(new RemoveIndelArtifactsFilter());
                if (!disableAtLeastQuarterFilter) {
                    genotypeFilters.add(new AtLeastAQuarterFilter());
                }
                break;
            default:
                throw new InternalError("Filters must be configured for new output format.");
        }
        System.out.println("Filtering reads that have these criteria:");
        for (GenotypeFilter filter : genotypeFilters) {
            System.out.println(filter.describe());
        }

        RandomAccessSequenceInterface genome = configureGenome(testGenome, jsapResult);

        int startFlapSize = jsapResult.getInt("start-flap-size", 100);
        formatConfigurator.configureFormatter(formatter);
        sortedPositionIterator = new DiscoverVariantIterateSortedAlignments(formatter);

        sortedPositionIterator.setGenome(genome);
        sortedPositionIterator.setStartFlapLength(startFlapSize);
        sortedPositionIterator.parseIncludeReferenceArgument(jsapResult);
        sortedPositionIterator.setMinimumVariationSupport(minimumVariationSupport);
        sortedPositionIterator.setThresholdDistinctReadIndices(thresholdDistinctReadIndices);
        return this;
    }


    /**
     * Parse the group-definition file, in the format sample-id=group-id (Java properties file)
     *
     * @param groupsDefinitionFile file with sample=group mapping information.
     * @param inputFilenames
     * @return group definition in the format expected by the --groups argument.
     */
    private String parseGroupFile(String groupsDefinitionFile, String[] inputFilenames) {
        Properties groupProps = new Properties();
        FileReader fileReader = null;
        try {
            fileReader = new FileReader(groupsDefinitionFile);
            groupProps.load(fileReader);
            Object2ObjectMap<String, MutableString> groupDefs = new Object2ObjectArrayMap<String, MutableString>();

            for (final Object key : groupProps.keySet()) {
                final String groupId = (String) groupProps.get(key);
                MutableString groupDef = groupDefs.get(groupId);
                if (groupDef == null) {
                    groupDef = new MutableString();
                    groupDefs.put(groupId, groupDef);
                }
                final String sampleId = (String) key;
                if (groupDef.indexOf(sampleId) == -1) {
                    if (isInputFilename(inputFilenames, sampleId)) {
                        groupDef.append(sampleId);
                        groupDef.append(',');
                    }
                }
            }
            MutableString result = new MutableString();
            for (final String groupId : groupDefs.keySet()) {
                result.append(groupId);
                result.append('=');
                result.append(groupDefs.get(groupId));
                // remove trailing coma:
                result.setLength(result.length() - 1);
                result.append('/');
            }
            result.setLength(result.length() - 1);
            System.out.println("generated group text: " + result);
            return result.toString();
        } catch (IOException e) {
            System.err.println("Cannot open or parse groups-file parameter " + groupsDefinitionFile);
            System.exit(1);
        } finally {

            IOUtils.closeQuietly(fileReader);
        }
        return null;
    }

    /**
     * Determine if an sampleId is provided on the command line.
     *
     * @param inputFilenames
     * @param sampleId
     * @return
     */
    private boolean isInputFilename(String[] inputFilenames, String sampleId) {
        for (final String input : inputFilenames) {
            String commandLineBasename = FilenameUtils.getBaseName(input);
            if (commandLineBasename.equals(FilenameUtils.getBaseName(sampleId))) {
                return true;
            }
        }
        return false;
    }

    public static RandomAccessSequenceInterface configureGenome(JSAPResult jsapResult) throws IOException {

        return configureGenome(null, jsapResult);
    }

    public static RandomAccessSequenceInterface configureGenome(RandomAccessSequenceInterface testGenome,
                                                                JSAPResult jsapResult) throws IOException {
        if (testGenome != null) {
            return testGenome;
        }
        String startOffsetArgument = jsapResult.getString("start-position");
        String endOffsetArgument = jsapResult.getString("end-position");
        String minIndex = getReferenceId(startOffsetArgument, "min");
        String maxIndex = getReferenceId(endOffsetArgument, "max");
        final String genome = jsapResult.getString("genome");
        RandomAccessSequenceCache cache = null;
        if (genome != null) {
            try {
                System.err.println("Loading genome cache " + genome);
                cache = new RandomAccessSequenceCache();
                cache.load(genome, minIndex, maxIndex);
                System.err.println("Done loading genome. ");
            } catch (ClassNotFoundException e) {
                System.err.println("Could not load genome cache");
                e.printStackTrace();
                System.exit(1);
            }
        }
        return cache;
    }

    private static String getReferenceId(String offsetArgument, String defaultValue) {
        if (offsetArgument == null) {
            return defaultValue;
        } else {
            String chr = offsetArgument.split(",")[0];
            return chr;
        }
    }

    public static AlignmentProcessorFactory configureProcessor(JSAPResult jsapResult) {
        final AlignmentProcessorNames processorName = AlignmentProcessorNames.valueOf(jsapResult.getString("processor").toUpperCase());
        AlignmentProcessorFactory realignmentFactory;
        switch (processorName) {
            case REALIGN_NEAR_INDELS:
                realignmentFactory = new AlignmentProcessorFactory() {
                    public AlignmentProcessorInterface create(final ConcatSortedAlignmentReader sortedReaders) {
                        return new LocalSortProcessor(new RealignmentProcessor(sortedReaders));
                    }
                };

                break;
            case NONE:
            default:
                System.err.println("Using processor NONE");
                realignmentFactory = new DefaultAlignmentProcessorFactory();

        }
        return realignmentFactory;
    }

    private void stopWhenDefaultGroupOptions() {
        if (groupComparisonsList.size() == 0) {
            System.err.println("Format group_comparison requires that arguments --group and --compare be defined.");
            System.exit(1);
        }
    }


    public void setTestGenome(final RandomAccessSequenceTestSupport testGenome) {
        this.testGenome = testGenome;
        overrideReferenceWithGenome = false;
    }

    /**
     * Install a format configurator, responsible for configuring the output format.
     *
     * @param configurator
     */
    public void setFormatConfigurator(FormatConfigurator configurator) {
        this.formatConfigurator = configurator;
    }

    public ArrayList<GroupComparison> getGroupComparisons() {
        return groupComparisonsList;
    }


    enum OutputFormat {
        VARIANT_DISCOVERY,
        ALLELE_FREQUENCIES,
        GENOTYPES,
        COMPARE_GROUPS,
        METHYLATION,
        BETWEEN_GROUPS,
        INDEL_COUNTS
    }

    enum AlignmentProcessorNames {
        NONE,
        REALIGN_NEAR_INDELS
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
        if (readIndexStats != null) {
            // build a samples array with the correct order:
            int numberOfSamples = readIndexStats.size();
            samples = new String[numberOfSamples];
            if (readIndexStats.size() == 0) {
                System.err.println("Cannot find any basename in stats file. Aborting.");
                System.exit(1);
            }
            for (ReadIndexStats stat : readIndexStats) {
                samples[stat.readerIndex] = stat.basename;
            }
            return samples;
        } else {
            // since we don't need to map basename order to readerIndex to just create a samples array from basenames
            // listed on the command line:
            samples = AlignmentReaderImpl.getBasenames(inputFilenames);
            // also remove the path to the file to keep only filenames:
            for (int i = 0; i < samples.length; i++) {
                samples[i] = FilenameUtils.getBaseName(samples[i]);
            }
            return samples;
        }
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
        } catch (IOException e) {
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
        sortedPositionIterator.initialize(this, outWriter, genotypeFilters);
        // install a reader factory that filters out ambiguous reads:
        sortedPositionIterator.setAlignmentReaderFactory(new NonAmbiguousAlignmentReaderFactory());
        sortedPositionIterator.setAlignmentProcessorFactory(realignmentFactory);
        sortedPositionIterator.setOverrideReferenceWithGenome(overrideReferenceWithGenome);
        sortedPositionIterator.setMaxThreshold(maxThresholdPerSite);
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
