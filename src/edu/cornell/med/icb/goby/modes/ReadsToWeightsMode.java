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
import edu.cornell.med.icb.goby.reads.Reads;
import edu.cornell.med.icb.goby.reads.ReadsReader;
import edu.cornell.med.icb.goby.algorithmic.data.HeptamerInfo;
import edu.cornell.med.icb.goby.algorithmic.data.WeightsInfo;
import edu.cornell.med.icb.goby.algorithmic.algorithm.*;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.*;

/**
 * Create the read to weight map. This class scans a compact reads file to determine which heptamer
 * occurs at the beginning of each read. When a heptamer is found that exists in the heptamer to weight
 * data structure provided as input (see heptamer-frequencies mode to generate this data structure),
 * the read is associated to the heptamer weight in a map. This map is written as Java serialized
 * file for use with modes that estimate gene/transcript/exon/other counts.
 *
 * @author Fabien Campagne
 *         Date: May 17 2009
 *         Time: 11:15 AM
 */
public class ReadsToWeightsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "reads-to-weights";
    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION = "Create a data structure that maps reads to a weights.";

    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(ReadsToWeightsMode.class);


    private String inputFilenames[];

    private String mapFilename;
    private String heptamerInfoFilename;
    private String estimationMethod;


    public String getModeName() {
        return MODE_NAME;
    }

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
    public AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilenames = jsapResult.getStringArray("input");
        heptamerInfoFilename = jsapResult.getString("heptamer-info");
        mapFilename = jsapResult.getString("map");
        estimationMethod = jsapResult.getString("method");

        return this;
    }

    enum WeightCalculationMethod {
        HEPTAMERS,
        G_PROPORTION,
        C_PROPORTION,
        A_PROPORTION,
        T_PROPORTION,
        GC_PROPORTION,
        AT_PROPORTION,
        ATGC_CORRECTION
    }

    @Override
    public void execute() throws IOException {

        HeptamerInfo heptamers = null;
        try {
            heptamers = HeptamerInfo.load(heptamerInfoFilename);

        } catch (ClassNotFoundException e) {
            System.err.println("Cannot load heptamer information from file " + heptamerInfoFilename);
            System.exit(1);
        }

        final ProgressLogger progress = new ProgressLogger();
        progress.start();
        progress.displayFreeMemory = true;
        WeightCalculator calculator = null;
        WeightCalculationMethod method = WeightCalculationMethod.valueOf(estimationMethod.toUpperCase());

        switch (method) {
            case HEPTAMERS:
                calculator = new HeptamerWeight(heptamers);

                break;
            case G_PROPORTION: {
                BaseProportionWeight calc = new BaseProportionWeight(heptamers);
                calc.setBase('G');
                calculator = calc;
                break;
            }

            case C_PROPORTION: {
                BaseProportionWeight calc = new BaseProportionWeight(heptamers);
                calc.setBase('C');
                calculator = calc;
                break;
            }
            case A_PROPORTION: {
                BaseProportionWeight calc = new BaseProportionWeight(heptamers);
                calc.setBase('A');
                calculator = calc;
                break;
            }
            case T_PROPORTION: {
                BaseProportionWeight calc = new BaseProportionWeight(heptamers);
                calc.setBase('T');
                calculator = calc;
                break;
            }
            case GC_PROPORTION:
                calculator = new GCProportionWeight(heptamers);
                break;
            case AT_PROPORTION:
                calculator = new ATProportionWeight(heptamers);
                break;
            case ATGC_CORRECTION:
                calculator = new ATGCCorrectionWeight(heptamers);
                break;
        }

        for (String inputFilename : inputFilenames) {
            // for each reads file:
            LOG.info("Now scanning " + inputFilename);

            if (inputFilenames.length > 1) {
                // if we process more than one reads file, build the map filename dynamically for each input file.
                mapFilename = FilenameUtils.removeExtension(inputFilename) + "." + calculator.id() + "-weights";
            }
            ReadsReader reader = new ReadsReader(new FileInputStream(inputFilename));
            try {
                WeightsInfo weights = new WeightsInfo();
                final MutableString sequence = new MutableString();
                int numberOfReads = 0;
                for (final Reads.ReadEntry readEntry : reader) {
                    ReadsReader.decodeSequence(readEntry, sequence);
                    final int readIndex = readEntry.getReadIndex();
                    float weight = calculator.weight(sequence);

                    weights.setWeight(readIndex, weight);

                    progress.lightUpdate();
                    numberOfReads++;
                }
                weights.size(numberOfReads);
                progress.stop();
                weights.save(mapFilename);

            } finally {

                if (reader != null) {
                    try {
                        reader.close();

                    } catch (IOException e) { // NOPMD
                        // silently ignore
                    }
                }
            }
        }


    }


    public static void main(final String[] args) throws IOException, JSAPException {
        new ReadsToWeightsMode().configure(args).execute();
    }
}
