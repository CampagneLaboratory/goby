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
import com.sun.org.apache.xpath.internal.operations.Mod;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.DensityEstimator;
import edu.cornell.med.icb.goby.algorithmic.algorithm.dmr.ObservationWriter;
import edu.cornell.med.icb.goby.stats.EmpiricalPValueEstimator;
import edu.cornell.med.icb.goby.stats.FormatFieldCounter;
import edu.cornell.med.icb.goby.util.DynamicOptionClient;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.io.FastBufferedReader;
import org.apache.commons.io.LineIterator;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

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
                    " - [element identifier]*" +
                    " - VALUES_A " +
                    " - [integer value to evaluate the statistic for element in sample A]+" +
                    " - VALUES_B " +
                    " - [integer value to evaluate the statistic for element in sample B]+" +
                    " - COVARIATES_A keyword" +
                    " - integer codes for covariates for sample A" +
                    " - COVARIATES_B keyword" +
                    " - integer codes for covariates for sample B";

    private String inputFilename;
    private String outputFilename;
    private String[] dymamicOptions;
    private String statisticName;
    public static final DynamicOptionClient doc = new DynamicOptionClient(EmpiricalPMode.class,
            EmpiricalPValueEstimator.LOCAL_DYNAMIC_OPTIONS
    );
    private DensityEstimator density;
    private String densityFilename;
    private boolean useExistingDensity;
    private boolean forceEstimation;

    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getModeDescription() {
        return MODE_DESCRIPTION;

    }

    static ObjectArrayList<DynamicOptionClient> registeredDOClients = new ObjectArrayList<DynamicOptionClient>();

    @Override
    public AbstractCommandLineMode configure(String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        outputFilename = jsapResult.getString("output");
        statisticName = jsapResult.getString("statistic");
        forceEstimation = jsapResult.getBoolean("force-estimation");
        registeredDOClients.add(MethylationRegionsOutputFormat.doc);
        {
            densityFilename = jsapResult.getString("density-filename");
            if (densityFilename != null) {

                if (new File(densityFilename).exists()) {
                    useExistingDensity = true;
                    // force loading of the pre-existing density and estimation of p-values:
                    doc.acceptsOption("EmpiricalPMode:serialized-estimator-filename=" + densityFilename);
                    if (forceEstimation) {
                        doc.acceptsOption("EmpiricalPMode:estimate-empirical-P=true");
                    }
                }

            }
        }
        if (statisticName != null) {
            boolean result = doc.acceptsOption("EmpiricalPMode:statistic=" + statisticName);
            assert result : "EmpiricalPMode:statistic= definition must be accepted.";
        }
        // parse dynamic options:
        dymamicOptions = jsapResult.getStringArray("dynamic-options");
        registeredDOClients.add(EmpiricalPMode.doc);
        for (final String dymamicOption : dymamicOptions) {
            boolean parsed = false;
            for (final DynamicOptionClient doc : registeredDOClients) {

                if (doc.acceptsOption(dymamicOption)) {

                    parsed = true;
                    break;
                }
            }
            if (!parsed) {
                System.err.println("Error: none of the installed tools could parse dynamic option: " + dymamicOption);
                System.exit(1);
            }
        }
        return this;
    }


    EmpiricalPValueEstimator estimator;
    FormatFieldCounter counter;

    @Override
    public void execute() throws IOException {
        LineIterator iterator = new LineIterator(new FastBufferedReader(new FileReader(inputFilename)));
        int lineNumber = 0;
        ObjectArrayList<String> elementIds = new ObjectArrayList<String>();

        IntArrayList valuesA = new IntArrayList();

        IntArrayList valuesB = new IntArrayList();
        IntArrayList covariatesA = new IntArrayList();
        IntArrayList covariatesB = new IntArrayList();
        estimator = new EmpiricalPValueEstimator();
        estimator.configure(1, doc);
        counter = new FormatFieldCounter(0, 2, 2, new String[]{"ALL"});
        // ignore the header line:
        iterator.next();
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
                for (j = 1; !"VALUES_A".equals(tokens[j]); j++) {
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
            process(typeOfPair, elementIds, valuesA, valuesB, covariatesA, covariatesB);
        }
        if (!useExistingDensity) {
            System.err.println("Writting estimated statistic to: " + densityFilename);
            DensityEstimator.store(estimator.getEstimator(), densityFilename);
        }

    }


    private void process(ObservationWriter.TypeOfPair typeOfPair, ObjectArrayList<String> elementIds, IntArrayList valuesA, IntArrayList valuesB,
                         IntArrayList covariatesA, IntArrayList covariatesB) {
        switch (typeOfPair) {
            case BETWEEN_GROUP_PAIR:
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

    String previousElementId = "";
    ObjectArrayList<IntArrayList> valuesACollector = new ObjectArrayList<IntArrayList>();
    ObjectArrayList<IntArrayList> valuesBCollector = new ObjectArrayList<IntArrayList>();
    ObjectArrayList<IntArrayList> covariatesACollector = new ObjectArrayList<IntArrayList>();
    ObjectArrayList<IntArrayList> covariatesBCollector = new ObjectArrayList<IntArrayList>();
    boolean first = true;

    private void estimateP(ObjectArrayList<String> elementIds, IntArrayList valuesA, IntArrayList valuesB, IntArrayList covariatesA, IntArrayList covariatesB) {
        //  System.out.println(elementIds);
        if (first) {
            previousElementId = elementIds.toString();
            first = false;
            System.out.print("ignore");
            int index=1;
            for (String id: elementIds) {

                System.out.print("\tid"+index++);
            }
            System.out.println("\tp-value");
        }
        if (!previousElementId.equals(elementIds.toString())) {

            final double p = estimator.estimateEmpiricalPValue(
                    valuesACollector.toArray(new IntArrayList[valuesACollector.size()]),
                    valuesBCollector.toArray(new IntArrayList[valuesBCollector.size()]),
                    covariatesACollector.toArray(new IntArrayList[covariatesACollector.size()]),
                    covariatesBCollector.toArray(new IntArrayList[covariatesBCollector.size()]));


            System.out.print("P-VALUE");
            for (final String elementId : elementIds) {

                System.out.print('\t');
                System.out.print(elementId);
            }
            System.out.printf("\t%g%n", p);
            valuesACollector.clear();
            valuesBCollector.clear();
            covariatesACollector.clear();
            covariatesBCollector.clear();
            previousElementId = elementIds.toString();

        }


        // continue to collect pair observation for the current element id
        valuesACollector.add(valuesA.clone());
        valuesBCollector.add(valuesB.clone());
        covariatesACollector.add(covariatesA.clone());
        covariatesBCollector.add(covariatesB.clone());


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
