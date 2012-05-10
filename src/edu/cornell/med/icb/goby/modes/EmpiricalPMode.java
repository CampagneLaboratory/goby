/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.EstimatedDistribution;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.ObservationWriter;
import edu.cornell.med.icb.goby.stats.EmpiricalPValueEstimator;
import edu.cornell.med.icb.goby.stats.EstimatedTestDistributions;
import edu.cornell.med.icb.goby.stats.FormatFieldCounter;
import edu.cornell.med.icb.goby.util.dynoptions.DynamicOptionClient;
import edu.cornell.med.icb.goby.util.dynoptions.RegisterThis;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.LineIterator;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Random;

/**
 * @author Fabien Campagne
 *         Date: 2/27/12
 *         Time: 11:00 AM
 */
public class EmpiricalPMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "empirical-p";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Estimate null distribution or empirical p-values based on a recorded null distribution. The input file must be tab delimited with the following fields: " +
                    " - WITHIN_GROUP_PAIR|BETWEEN_GROUP_PAIR keyword" +
                    " - comparison_description " +
                    " - [element identifier]*" +
                    " - VALUES_A " +
                    " - [integer value to evaluate the statistic for element in sample A]+" +
                    " - VALUES_B " +
                    " - [integer value to evaluate the statistic for element in sample B]+" +
                    " - COVARIATES_A keyword" +
                    " - integer codes for covariates for sample A" +
                    " - COVARIATES_B keyword" +
                    " - integer codes for covariates for sample B";
    /**
     * Used to log debug and informational messages.
     */
    private static final Logger LOG = Logger.getLogger(EmpiricalPValueEstimator.class);

    private String inputFilename;
    private String outputFilename;
    private String[] dymamicOptions;
    private String statisticName;
    private int index = 0;
    // number of false negatives (assuming a simulation where the first 1000 elements tested are positives.
    private int numFN;
    private int numFP;
    private int numTN;
    private int numTP;
    private int numOther;
    private EstimatedTestDistributions testDensities;
    /**
     * When True, write an estimate of false discovery rate, rather than the empirical P (FWER).
     */
    private boolean fdr;
    /**
     * The number of test pairs to consider.
     */
    private int testN;


    public static DynamicOptionClient doc() {
        return doc;
    }


    public static final
    @RegisterThis
    DynamicOptionClient doc = new DynamicOptionClient(EmpiricalPMode.class,
            EmpiricalPValueEstimator.LOCAL_DYNAMIC_OPTIONS

    );
    private EstimatedDistribution density;
    private String densityFilename;
    private boolean useExistingDensity;
    private boolean forceEstimation;
    private PrintWriter outputWriter;
    private ObjectArrayList<String> previousElementIdList;
    private String label;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;

    }


    @Override
    public AbstractCommandLineMode configure(String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        label = jsapResult.getString("label");
        outputFilename = jsapResult.getString("output");
        statisticName = jsapResult.getString("statistic");
        forceEstimation = jsapResult.getBoolean("force-estimation");
        //fdr = jsapResult.getBoolean("fdr");

        testN = jsapResult.getInt("test-n",Integer.MAX_VALUE);

        {
            densityFilename = jsapResult.getString("density-filename");
            if (densityFilename != null) {

                if (new File(densityFilename).exists()) {
                    try {
                        // if we are given a serialized statistic, we override the stat name with that of the actual
                        // statistic used to estimate the density:
                        density = EstimatedDistribution.load(densityFilename);
                        //testDensities = EstimatedTestDistributions.load(EmpiricalPValueEstimator.suggestFilename(densityFilename));
                        statisticName = density.getStatAdaptor().statName();
                        /*  assert density.getStatAdaptor().statName().equals(testDensities.getStatAdaptor().statName()) :
                        "null and test distributions must be estimated with the same statistic";
                        */
                    } catch (ClassNotFoundException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        System.exit(1);
                    }
                    useExistingDensity = true;
                    // force loading of the pre-existing density and estimation of p-values:
                    doc.acceptsOption("EmpiricalPMode:serialized-estimator-filename=" + densityFilename);
                    if (forceEstimation) {
                        doc.acceptsOption("EmpiricalPMode:estimate-empirical-P=true");
                        doc.acceptsOption("EmpiricalPMode:estimate-intra-group-differences=false");
                    } else {
                        doc.acceptsOption("EmpiricalPMode:estimate-empirical-P=false");
                        doc.acceptsOption("EmpiricalPMode:estimate-intra-group-differences=true");
                    }
                }
            }


        }

        if (statisticName != null) {
            boolean result = doc.acceptsOption("EmpiricalPMode:statistic=" + statisticName);
            assert result : "EmpiricalPMode:statistic= definition must be accepted.";
        }

        return this;
    }

    private void setupOutput() throws FileNotFoundException {
        if (outputFilename == null) {
            String binningName = estimator.getNullDistribution().getBinningStrategy().getName();
            outputFilename = FilenameUtils.getBaseName(inputFilename) + "-" + statisticName + "-" + binningName + "-" + doc.getString("combinator") + ".tsv";
            System.out.println("Output will be written to " + outputFilename);
        }
        outputWriter = new PrintWriter(outputFilename);
    }

    private void constructDensityFilename(String densityFilename) {
        if (this.densityFilename == null) {
            // construct a density filename that reflects the arguments that control the density:
            String binningName = estimator.getNullDistribution().getBinningStrategy().getName();
            this.densityFilename = FilenameUtils.getBaseName(inputFilename) + "-" + statisticName + "-" + binningName + "-density.bin";
        } else {
            this.densityFilename = densityFilename;
        }
    }


    private EmpiricalPValueEstimator estimator;
    private FormatFieldCounter counter;

    @Override
    public void execute() throws IOException {
        estimator = new EmpiricalPValueEstimator();
        estimator.configure(1, doc);
        estimator.setEstimateFdr(fdr);
        scan();

        if (forceEstimation) {
            write(previousElementIdList);
        }
        outputWriter.flush();
        if (!useExistingDensity) {

            constructDensityFilename(densityFilename);
            System.err.println("Writing estimated null distribution to: " + densityFilename);
            EstimatedDistribution.store(estimator.getNullDistribution(), densityFilename);

            final String testFilename = EmpiricalPValueEstimator.suggestFilename(densityFilename);
            /*    System.err.println("Writing estimated test distributions to: " + testFilename);
            EstimatedTestDistributions.store(estimator.getTestDistributions(), testFilename);
            */

        }
        System.out.printf("Rate of true negatives: %3.2f%%%n", 100d * ((numTN + 0d) / (numOther + 0d)));
        System.out.printf("Rate of true positives: %3.2f%%%n", 100d * ((numTP + 0d) / (1000d)));

        System.out.printf("Rate of false negatives: %3.2f%%%n", 100d * ((numFN + 0d) / 1000d));
        System.out.printf("Rate of false positives: %3.2f%%%n", 100d * ((numFP + 0d) / (numOther + 0d)));
        System.out.printf("Overall false discovery rate: %3.2f%%%n", (1d - (0d + numTP) / (0d + numTP + numFP)) * 100d);

    }

    private void scan() throws FileNotFoundException {
        LineIterator iterator = new LineIterator(new FastBufferedReader(new FileReader(inputFilename)));
        int lineNumber = 0;
        ObjectArrayList<String> elementIds = new ObjectArrayList<String>();

        IntArrayList valuesA = new IntArrayList();

        IntArrayList valuesB = new IntArrayList();
        IntArrayList covariatesA = new IntArrayList();
        IntArrayList covariatesB = new IntArrayList();

        counter = new FormatFieldCounter(0, 2, 2, new String[]{"ALL"});
        setupOutput();
        // ignore the header line:
        iterator.next();
        ProgressLogger pg = new ProgressLogger(LOG);
        pg.displayFreeMemory = true;
        pg.itemsName = "pairs";
        pg.expectedUpdates = countLines(inputFilename) - 1;
        pg.start("Starting to scan pairs.");
        while (iterator.hasNext()) {
            String next = iterator.nextLine();
            String[] tokens = next.split("\t");
            boolean pastIds = false;
            boolean pastValues = false;

            String typeOfPairString = tokens[0];
            ObservationWriter.TypeOfPair typeOfPair = ObservationWriter.TypeOfPair.UNDEFINED;

            for (int i = 0; i < tokens.length; i++) {

                try {
                    typeOfPair = ObservationWriter.TypeOfPair.valueOf(typeOfPairString);
                } catch (IllegalArgumentException e) {
                    System.err.println("First token of every line should be WITHIN_GROUP_PAIR or BETWEEN_GROUP_PAIR. Found " + typeOfPairString + " on line " + lineNumber);
                    System.exit(1);
                }

                elementIds.clear();
                valuesA.clear();
                valuesB.clear();
                covariatesA.clear();
                covariatesB.clear();
                int j;
                String groupComparison = tokens[1];
                elementIds.add(groupComparison);
                for (j = 2; !"VALUES_A".equals(tokens[j]); j++) {
                    if (j == tokens.length) {
                        break;
                    }
                    elementIds.add(tokens[j]);
                }
                if (j == tokens.length) {
                    System.err.println("Every line must contain the VALUES keyword. Keyword not found on line " + lineNumber);
                    System.exit(1);
                }
                j++;
                for (; !"VALUES_B".equals(tokens[j]); j++) {
                    if (j == tokens.length) {
                        break;
                    }
                    valuesA.add(Integer.parseInt(tokens[j]));
                }
                j++;
                for (; !"COVARIATES_A".equals(tokens[j]); j++) {
                    if (j == tokens.length) {
                        break;
                    }
                    valuesB.add(Integer.parseInt(tokens[j]));
                }
                if (j == tokens.length) {
                    System.err.println("Every line must contain the COVARIATES_A keyword. Keyword not found on line " + lineNumber);
                    System.exit(1);
                }
                j++;
                for (; !"COVARIATES_B".equals(tokens[j]); j++) {
                    if (j == tokens.length) {
                        break;
                    }
                    covariatesA.add(Integer.parseInt(tokens[j]));
                }
                if (j == tokens.length) {
                    System.err.println("Every line must contain the COVARIATES_B keyword. Keyword not found on line " + lineNumber);
                    System.exit(1);
                }
                j++;
                for (; j < tokens.length; j++) {
                    covariatesB.add(Integer.parseInt(tokens[j]));
                }


            }
            lineNumber++;
            final String groupComparison = elementIds.get(0);
            process(typeOfPair, groupComparison, elementIds, valuesA, valuesB, covariatesA, covariatesB);
            pg.lightUpdate();
        }
        pg.done(lineNumber);
    }

    private int countLines(String inputFilename) throws FileNotFoundException {
        int lineCount = 0;
        LineIterator it = new LineIterator(new FileReader(inputFilename));
        while (it.hasNext()) {
            Object next = it.next();
            lineCount++;
        }
        it.close();
        return lineCount;
    }


    private void process(ObservationWriter.TypeOfPair typeOfPair, String groupComparison,
                         ObjectArrayList<String> elementIds, IntArrayList valuesA, IntArrayList valuesB,
                         IntArrayList covariatesA, IntArrayList covariatesB) {
        switch (typeOfPair) {
            case BETWEEN_GROUP_PAIR:
                if (!useExistingDensity && testDensities != null) {
                    observeTestDistributions(groupComparison, valuesA, valuesB, covariatesA, covariatesB);
                }
                estimateP(elementIds, valuesA, valuesB, covariatesA, covariatesB);

                break;
            case WITHIN_GROUP_PAIR:
                if (useExistingDensity && forceEstimation) {
                    estimateP(elementIds, valuesA, valuesB, covariatesA, covariatesB);
                } else {
                    observeNullDistribution(valuesA, valuesB, covariatesA, covariatesB);
                }
                break;
            default:
                System.out.println("typeofPair is unsupported at this time: " + typeOfPair);
                System.exit(1);
        }
    }

    private void observeNullDistribution(IntArrayList valuesA, IntArrayList valuesB, IntArrayList covariatesA, IntArrayList covariatesB) {

        estimator.estimateNullDensity(valuesA, valuesB, covariatesA, covariatesB);
    }

    private void observeTestDistributions(String groupComparison, IntArrayList valuesA, IntArrayList valuesB, IntArrayList covariatesA, IntArrayList covariatesB) {

        estimator.estimateTestDensity(groupComparison, valuesA, valuesB, covariatesA, covariatesB);
    }

    String previousElementId = "";

    ObjectArrayList<IntArrayList> valuesACollector = new ObjectArrayList<IntArrayList>();
    ObjectArrayList<IntArrayList> valuesBCollector = new ObjectArrayList<IntArrayList>();
    ObjectArrayList<IntArrayList> covariatesACollector = new ObjectArrayList<IntArrayList>();
    ObjectArrayList<IntArrayList> covariatesBCollector = new ObjectArrayList<IntArrayList>();
    boolean first = true;

    private void estimateP(ObjectArrayList<String> elementIds, IntArrayList valuesA, IntArrayList valuesB, IntArrayList covariatesA, IntArrayList covariatesB) {

        //  System.out.println(elementIds);
        if (first) {
            outputWriter.print("LABEL");
            previousElementId = elementIds.toString();
            previousElementIdList = elementIds.clone();
            first = false;
            outputWriter.print("\tignore");
            int index = 1;
            for (String id : elementIds) {

                outputWriter.print("\tid" + index++);
            }
            outputWriter.println("\tp-value");
        }
        if (!previousElementId.equals(elementIds.toString())) {

            write(previousElementIdList);
            valuesACollector.clear();
            valuesBCollector.clear();
            covariatesACollector.clear();
            covariatesBCollector.clear();
            previousElementId = elementIds.toString();
            previousElementIdList = elementIds.clone();
        }


        // continue to collect pair observation for the current element id
        valuesACollector.add(valuesA.clone());
        valuesBCollector.add(valuesB.clone());
        covariatesACollector.add(covariatesA.clone());
        covariatesBCollector.add(covariatesB.clone());


    }

    Random random = new Random();

    private void subsample(ObjectArrayList<IntArrayList> valuesA, ObjectArrayList<IntArrayList> covariatesA, int testN) {

        int size = valuesA.size();

        assert size == covariatesA.size() :        "value and covariate size must match!";
        boolean[] remove = new boolean[size];
        boolean done = false;
        int removedCount = 0;
        while (!done) {
            int index =Math.round( random.nextFloat() * (size - 1.0f));
            remove[index] = true;

            if (removedCount(remove) == size-testN) {
                done = true;
            }
        }
        int removed=0;
        for (int i = 0; i < size; i++) {
            int index=i-removed++;
            if (remove[i]) {
                valuesA.remove(index);
                covariatesA.remove(index);
            }
        }
  //      System.out.println("STOP");
    }

    private int removedCount(boolean[] remove) {
        int count=0;
        for (boolean removed: remove) {
            count+=removed?1:0;
        }
        return count;
    }

    private void write(ObjectArrayList<String> elementIds) {
        outputWriter.print(label);
        outputWriter.print("\tP-VALUE");
        boolean stop = false;
        for (final String elementId : elementIds) {

            /*  if ("ENSMUST00000034529_NM_027807".equals(elementId)) {
               System.out.println("STOP");
               stop = true;
           } else {
               stop = false;
           } */
            outputWriter.print('\t');
            outputWriter.print(elementId);
        }      /*
        if (stop) {
            System.out.println("STOP");
        }        */
        if (elementIds.get(5).equals("A104")) {
            System.out.println(elementIds);
        }
        final String groupComparison = elementIds.get(0);

        if (valuesACollector.size() >= testN) {
            subsample(valuesACollector, covariatesACollector, testN);
            subsample(valuesBCollector, covariatesBCollector, testN);
        }

        final double p = estimator.estimateEmpiricalP(groupComparison,
                valuesACollector.toArray(new IntArrayList[valuesACollector.size()]),
                valuesBCollector.toArray(new IntArrayList[valuesBCollector.size()]),
                covariatesACollector.toArray(new IntArrayList[covariatesACollector.size()]),
                covariatesBCollector.toArray(new IntArrayList[covariatesBCollector.size()]));

        final double alpha = 0.05;
        if (index < 1000) {
            if (p > alpha) {
                numFN++;
            } else {
                numTP++;
            }
        }
        if (index >= 1000) {
            if (p <= alpha) {
                System.out.println("STOP FP "+elementIds+" p="+p);
                numFP++;
            } else {
                numTN++;
            }
            numOther++;

        }

        index++;
        outputWriter.printf("\t%g%n", p);
    }

    protected EmpiricalPMode(final String jarFilename) {
        super(jarFilename);
    }

    protected EmpiricalPMode() {
        super();
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws JSAPException error parsing
     * @throws IOException   error parsing or executing.
     */

    public static void main(final String[] args) throws JSAPException, IOException {
        new EmpiricalPMode().configure(args).execute();
    }
}
