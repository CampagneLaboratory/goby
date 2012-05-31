/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.algorithmic.data.xml.AnnotationLength;
import edu.cornell.med.icb.goby.algorithmic.data.xml.InfoOutput;
import edu.cornell.med.icb.goby.algorithmic.data.xml.SampleTotalCount;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionAnalysis;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionCalculator;
import edu.cornell.med.icb.goby.stats.DifferentialExpressionResults;
import edu.cornell.med.icb.goby.stats.NormalizationMethod;
import edu.cornell.med.icb.io.TSVReader;
import edu.rit.pj.ParallelTeam;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.io.*;
import java.util.ArrayList;

/**
 * Estimate statistics for a table produced with alignment-to-counts. Varied statistics can be produced, including
 * RPKMs, statistics of differential expression, etc.
 *
 * @author Fabien Campagne
 */
public class StatsMode extends AbstractGobyMode {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(StatsMode.class);

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "stats";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Estimate statistics for a table produced with alignment-to-counts";

    private String inputFilename;

    private boolean writeAnnotationCounts = true;
    private boolean omitNonInformativeColumns;
    private String outputFilename;
    private ParallelTeam team;
    private boolean parallel;
    private boolean doComparison;
    private final DifferentialExpressionCalculator deCalculator = new DifferentialExpressionCalculator();
    private final DifferentialExpressionAnalysis deAnalyzer = new DifferentialExpressionAnalysis();
    /**
     * The set of normalization methods to use for the comparison.
     */

    private ObjectArraySet<NormalizationMethod> normalizationMethods;
    private String[] sampleIds;
    private int numElements;
    private String infoFilename;
    private ArrayList<GroupComparison> groupComparisonList;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getShortModeName() {
        return null;
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

        parallel = jsapResult.getBoolean("parallel", false);
        inputFilename = jsapResult.getString("input");

        outputFilename = jsapResult.getString("output");
        final String groupsDefinition = jsapResult.getString("groups");
        sampleIds = getSampleIds(inputFilename);
        deAnalyzer.parseGroupsDefinition(groupsDefinition, deCalculator, sampleIds);
        final String compare = jsapResult.getString("compare");
        doComparison = compare != null;
        if (doComparison) {
            groupComparisonList=deAnalyzer.parseCompare(compare);
        }
        deAnalyzer.setRunInParallel(parallel);

        normalizationMethods = deAnalyzer.parseNormalization(jsapResult);
        parseEval(jsapResult, deAnalyzer);
        infoFilename = jsapResult.getString("info");
        return this;
    }

    private int countLines(final String inputFilename) throws FileNotFoundException {
        int lineNum = 0;
        final LineIterator lineIt = new LineIterator(new FastBufferedReader(new FileReader(inputFilename)));
        while (lineIt.hasNext()) {
            lineIt.next();
            lineNum++;
        }

        return lineNum;
    }

    private String[] getSampleIds(final String inputFilename) throws FileNotFoundException {
        final LineIterator lineIt = new LineIterator(new FastBufferedReader(new FileReader(inputFilename)));
        if (!lineIt.hasNext()) {
            System.err.println("Input file must have at least one line, tab delimited, in the format: element-id element-type element-length (sample-id)+");
        }
        final MutableString firstLine = lineIt.next();

        final String[] tokens = firstLine.toString().split("\t");
        final String[] result = new String[tokens.length - 2];

        System.arraycopy(tokens, 2, result, 0, tokens.length - 2);
        return result;
    }

    public static void parseEval(final JSAPResult jsapResult, final DifferentialExpressionAnalysis deAnalyzer) {
        final String evalString = jsapResult.getString("eval");

        final String[] evalArray = evalString.split(",");
        final ObjectSet<String> evalSet = new ObjectOpenHashSet<String>();
        for (final String evalName : evalArray) {
            evalSet.add(evalName.trim().toLowerCase().intern());
        }
        deAnalyzer.setEvalNames(evalSet);
    }


    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        numElements = countLines(inputFilename) - 1;
        loadInput(inputFilename, numElements, sampleIds);

        if (infoFilename != null) {
            loadInfo();

        } else {
            System.err.println("Warning: the info option must be provided to obtain accurate statistics, including RPKMs and statistics derived from RPKMs. Info files can be generated with the alignment-to-counts mode. ");
        }

        // estimate statistics:
        final DifferentialExpressionResults results =
                deAnalyzer.evaluateDifferentialExpressionStatistics(deCalculator, doComparison, normalizationMethods);

        // write results:
        PrintWriter statsWriter = null;
        try {
            statsWriter = new PrintWriter(outputFilename);
            results.write(statsWriter, '\t', deCalculator);
        } finally {
            IOUtils.closeQuietly(statsWriter);
        }

    }

    private void loadInfo() {
        try {
            final JAXBContext jc = JAXBContext.newInstance(InfoOutput.class);

            final Unmarshaller m = jc.createUnmarshaller();
            InfoOutput info = (InfoOutput) m.unmarshal(new File(infoFilename));
            for (AnnotationLength ae : info.lengths) {
                int index = deCalculator.getElementIndex(ae.id);
                if (index != -1) {

                    deCalculator.defineElementLength(index, ae.length);
                } else {

                    // OK since some elements will yield zero counts and not be in the input.
                }
            }
            for (SampleTotalCount tc : info.totalCounts) {
                deCalculator.setNumAlignedInSample(tc.sampleId, tc.totalCount);
            }
        } catch (JAXBException e) {
            System.err.printf("An error occurred loading the content of info file %s %n", infoFilename);
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Load the input file into deCalculator.
     *
     * @param inputFilename tab-delimited intut file produced by alignment-to-annotation-counts
     * @param sampleIds     sample identifiers.
     * @throws IOException if an error occurs reading the input file.
     */
    private void loadInput(final String inputFilename, final int numElements, final String[] sampleIds) throws IOException {
        final TSVReader reader = new TSVReader(new FileReader(inputFilename), '\t');
        // skip header line:
        if (reader.hasNext()) {
            reader.next();
        }

        deCalculator.reserve(numElements, sampleIds.length);
        int numLines = 0;
        // read counts for each next line:
        while (reader.hasNext()) {
            if (reader.isCommentLine()) {
                reader.skip();
            } else {
                reader.next();
                numLines++;
                final String elementId = reader.getString();
                final String elementType = reader.getString();
                if (!"element-id".equals(elementId)) {
                    deCalculator.defineElement(elementId, DifferentialExpressionCalculator.ElementType.valueOf(elementType));

                    for (final String sample : sampleIds) {
                        deCalculator.observe(sample, elementId, reader.getDouble());
                    }
                }
            }
        }

    }


    /**
     * Returns the log2 of the sum of the argument plus one.
     *
     * @param x argument.
     * @return log2p1(x)=Math.log1p(x)/Math.log(2)
     */
    private double log2p1(final double x) {
        return StrictMath.log1p(x) / LOG_2;
    }

    /**
     * The natural log of the number two.
     */
    private static final double LOG_2 = StrictMath.log(2);

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     * @throws java.io.IOException error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new StatsMode().configure(args).execute();
    }

}
