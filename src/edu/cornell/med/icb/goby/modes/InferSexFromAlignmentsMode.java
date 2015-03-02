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
import edu.cornell.med.icb.goby.util.Timer;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Converts a compact alignment to a compressed count archive.
 *
 * @author Fabien Campagne
 */
public class InferSexFromAlignmentsMode extends AbstractGobyMode {

    /**
     * The mode name.
     */
    private static final String MODE_NAME = "infer-sex";

    /**
     * The overridden short mode name.
     */
    private static final String SHORT_MODE_NAME = "is";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Infer sex from data for human samples.";

    /**
     * The output file.
     */
    private String outputFile;


    private String[] basenames;
    private String optionalOutputFile;


    private static final Logger LOG = Logger.getLogger(InferSexFromAlignmentsMode.class);
    private boolean verbose;


    @Override
    public String getModeName() {
        return MODE_NAME;
    }

    @Override
    public String getShortModeName() {
        return SHORT_MODE_NAME;
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
     * @throws java.io.IOException                    error parsing
     * @throws com.martiansoftware.jsap.JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        final String[] inputFiles = jsapResult.getStringArray("input");
        final ObjectSet<String> basenameSet = new ObjectOpenHashSet<String>();

        for (final String inputFile : inputFiles) {
            basenameSet.add(AlignmentReaderImpl.getBasename(inputFile));
        }

        basenames = basenameSet.toArray(new String[basenameSet.size()]);
        optionalOutputFile = jsapResult.getString("output");
        if (optionalOutputFile == null) {
            optionalOutputFile = "sex.tsv";
        }
        outputWriter = new PrintWriter(new FileWriter(optionalOutputFile));
        outputWriter.print("basename\thitsOnChrX\thitsOnChrY\tratioXOverY\tfractionOnX\tfractionOnY\tsex\n");
        return this;
    }

    PrintWriter outputWriter;

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        final ProgressLogger progress = new ProgressLogger(LOG);
        progress.priority = Level.INFO;
        progress.expectedUpdates = basenames.length;
        progress.itemsName = "basenames";
        progress.start();
        if (basenames.length > 1) {
            verbose = true;
        }
        for (final String basename : basenames) {
            if (optionalOutputFile == null) {
                outputFile = basename;
            }

            processFullGenomeAlignment(basename);
            progress.info = "last processed: " + basename;
            progress.lightUpdate();

        }
        progress.done();
        outputWriter.close();
    }

    private void processFullGenomeAlignment(final String basename) throws IOException {
        final AlignmentReaderFactory factory = new DefaultAlignmentReaderFactory();
        System.out.println("Processing " + basename);
        // don't use the factory here, we just need to read the header and AlignmentReaderImpl will be faster:
        AlignmentReader reader = new AlignmentReaderImpl(basename);

        reader.readHeader();

        final DoubleIndexedIdentifier referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
        reader.close();


        final AlignmentReader referenceReader = new AlignmentReaderImpl(basename);
        referenceReader.readHeader();
        referenceReader.readIndex();
        int chrXIndex = referenceIds.getIndex("chrX");
        int chrYIndex = referenceIds.getIndex("chrY");
        if (chrXIndex == -1) {
            chrXIndex = referenceIds.getIndex("X");
        }
        if (chrYIndex == -1) {
            chrYIndex = referenceIds.getIndex("Y");
        }


        long hitsOnChrX = 0;
        long hitsOnChrY = 0;

        final Timer timer = new Timer();
        timer.start();
        // jump to the first sex chromosome:
        int skipToIndex = Math.min(chrXIndex, chrYIndex);
        Alignments.AlignmentEntry alignmentEntry = null;

        while ((alignmentEntry = referenceReader.skipTo(skipToIndex, 0)) != null) {
            final int referenceIndex = alignmentEntry.getTargetIndex();
            if (referenceIndex == chrXIndex) {
                hitsOnChrX += 1;
            } else if (referenceIndex == chrYIndex) {
                hitsOnChrY += 1;
            } else {
                final int maxSexChrIndex = Math.max(chrXIndex, chrYIndex);
                if (referenceReader.isSorted() && referenceIndex < maxSexChrIndex) {
                    // jump to the second sex chromosome:
                    skipToIndex = maxSexChrIndex;
                }
            }
        }

        double ratioXOverY = ((double) hitsOnChrX) / (double) hitsOnChrY;
        reader.close();
        double fractionOnX=(double)hitsOnChrX/ (double)referenceReader.getNumberOfAlignedReads();
        double fractionOnY=(double)hitsOnChrY/ (double)referenceReader.getNumberOfAlignedReads();
        boolean inferFemale = ratioXOverY >= 100;
        outputWriter.printf("%s\t%d\t%d\t%g\t%g\t%g\t%s%n", FilenameUtils.getBaseName(basename), hitsOnChrX, hitsOnChrY,
                ratioXOverY, fractionOnX, fractionOnY,
                inferFemale ? "female" : "male");
        outputWriter.flush();
        timer.stop();
        System.out.println(timer);
    }

    /**
     * Main method.
     *
     * @param args command line args.
     * @throws com.martiansoftware.jsap.JSAPException error parsing
     * @throws java.io.IOException                    error parsing or executing.
     */
    public static void main(final String[] args) throws JSAPException, IOException {
        new InferSexFromAlignmentsMode().configure(args).execute();
    }


}

