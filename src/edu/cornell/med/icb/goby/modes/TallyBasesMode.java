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
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.counts.AnyTransitionCountsIterator;
import edu.cornell.med.icb.goby.counts.CountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountsReader;
import edu.cornell.med.icb.goby.counts.OffsetCountsReader;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceCache;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;
import java.util.zip.GZIPInputStream;

/**
 *
 */
public class TallyBasesMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    public static final String MODE_NAME = "tally-bases";

    /**
     * The mode description help text.
     */
    public static final String MODE_DESCRIPTION = "";

    /**
     * The input file.
     */
    private String[] inputFilenames;

    /**
     * The output file.
     */
    private String outputFilename;
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();
    private String[] basenames;
    private String alternativeCountArhive;
    private String genomeFilename;
    private String genomeCacheFilename;
    private final int windowSize = 10;
    private String offsetString;
    private double cutoff;

    /**
     * From 0 to 1.
     */
    private double sampleRate;

    @Override
    public String getModeName() {
        return MODE_NAME;
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
     * @throws IOException error parsing
     * @throws JSAPException error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilenames = jsapResult.getStringArray("input");
        basenames = AlignmentReader.getBasenames(inputFilenames);
        alternativeCountArhive = jsapResult.getString("alternative-count-archive");
        outputFilename = jsapResult.getString("output");
        genomeFilename = jsapResult.getString("genome");
        genomeCacheFilename = jsapResult.getString("genome-cache");
        offsetString = jsapResult.getString("offset");
        cutoff = jsapResult.getDouble("cutoff");
        sampleRate = jsapResult.getDouble("sample-rate", 100d) / 100d;
        if (sampleRate == 0d) {
            System.err.println("Sample rate must be larger than zero or no positions will be reported.");
            System.exit(1);
        }
        if ("auto".equals(genomeCacheFilename)) {
            String filename = FilenameUtils.getBaseName(genomeFilename);
            filename = FilenameUtils.removeExtension(filename);
            genomeCacheFilename = filename + "-cache";
        }


        final String includeReferenceNameCommas = jsapResult.getString("include-reference-names");
        if (includeReferenceNameCommas != null) {
            includeReferenceNames = new ObjectOpenHashSet<String>();
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameCommas.split("[,]")));
            System.out.println("Will tally bases for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
        }
        return this;
    }

    /**
     * Run the tally bases mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() throws IOException {
        if (basenames.length != 2) {
            System.err.println("Exactly two basenames are supported at this time.");
            System.exit(1);
        }
        final CountsArchiveReader[] archives = new CountsArchiveReader[basenames.length];
        int i = 0;
        for (final String basename : basenames) {
            archives[i++] = new CountsArchiveReader(basename, alternativeCountArhive);
        }

        final CountsArchiveReader archiveA = archives[0];
        final CountsArchiveReader archiveB = archives[1];
        // keep only common reference sequences between the two input count archives.
        final ObjectSet<String> identifiers = new ObjectOpenHashSet<String>();
        identifiers.addAll(archiveA.getIdentifiers());
        identifiers.retainAll(archiveB.getIdentifiers());
        // find the optimal offset A vs B:
        final int offset = offsetString.equals("auto") ? optimizeOffset(archiveA, archiveB, identifiers) : Integer.parseInt(offsetString);
        System.out.println("offset: " + offset);

        final RandomAccessSequenceCache cache = new RandomAccessSequenceCache();
        if (cache.canLoad(genomeCacheFilename)) {
            try {
                cache.load(genomeCacheFilename);
            } catch (ClassNotFoundException e) {
                System.err.println("Cannot load cache from disk. Consider deleting the cache and rebuilding.");
                e.printStackTrace();
                System.exit(1);
            }
        } else {
            Reader reader = null;
            try {
                if (genomeFilename.endsWith(".fa") || genomeFilename.endsWith(".fasta")) {
                    reader = new FileReader(genomeFilename);
                    cache.loadFasta(reader);
                } else if (genomeFilename.endsWith(".fa.gz") || genomeFilename.endsWith(".fasta.gz")) {
                    reader = new InputStreamReader(new GZIPInputStream(new FileInputStream(genomeFilename)));
                    cache.loadFasta(reader);
                } else {
                    System.err.println("The format of the input file is not supported at this time.");
                    System.exit(1);
                }
            } finally {
                IOUtils.closeQuietly(reader);
            }
        }

        System.out.println("Will use genome cache basename: " + genomeCacheFilename);
        cache.save(genomeCacheFilename);
        final Random random = new Random(new Date().getTime());

        final double delta = cutoff;
        final int countThreshold = 30;
        final PrintStream output = new PrintStream(outputFilename);
        writeHeader(output, windowSize);
        for (final String referenceSequenceId : identifiers) {
            if (isReferenceIncluded(referenceSequenceId)) {

                final int referenceIndex = cache.getReferenceIndex(referenceSequenceId);
                if (referenceIndex != -1) {
                    // sequence in cache.
                    System.out.println("Processing sequence " + referenceSequenceId);
                    final double sumA = getSumOfCounts(archiveA.getCountReader(referenceSequenceId));
                    final double sumB = getSumOfCounts(archiveB.getCountReader(referenceSequenceId));
                    final int referenceSize = cache.getSequenceSize(referenceIndex);
                    // process this sequence:
                    final AnyTransitionCountsIterator iterator = new AnyTransitionCountsIterator(
                            archiveA.getCountReader(referenceSequenceId),
                            new OffsetCountsReader(archiveB.getCountReader(referenceSequenceId), offset));

                    while (iterator.hasNextTransition()) {
                        iterator.nextTransition();
                        final int position = iterator.getPosition();
                        final int countA = iterator.getCount(0);
                        final int countB = iterator.getCount(1);

                        if (countA + countB >= countThreshold) {
                            final double foldChange = Math.log1p(countA) - Math.log1p(countB) - Math.log(sumA) + Math.log(sumB);
                            if (foldChange >= delta || foldChange <= -delta) {
                                if (random.nextDouble() < sampleRate) {
                                    tallyPosition(cache, referenceIndex, position, foldChange, windowSize,
                                            referenceSize,
                                            referenceSequenceId,
                                            output, countA, countB, sumA, sumB);
                                }
                            }
                        }
                    }
                    iterator.close();
                }
            }
            output.flush();
        }
        output.close();
    }

    private int optimizeOffset(final CountsArchiveReader archiveA, final CountsArchiveReader archiveB, final ObjectSet<String> identifiers) throws IOException {

        long minDifferenceSeen = Long.MAX_VALUE;
        int optimalOffset = -1;
        for (int offset = -10; offset <= 10; offset++) {
            System.out.println("Trying offset=" + offset);
            long differenceSeen = 0;
            for (final String referenceSequenceId : identifiers) {
                if (isReferenceIncluded(referenceSequenceId)) {


                    // process this sequence:
                    final AnyTransitionCountsIterator iterator = new AnyTransitionCountsIterator(
                            archiveA.getCountReader(referenceSequenceId),
                            new OffsetCountsReader(archiveB.getCountReader(referenceSequenceId), offset));
                    while (iterator.hasNextTransition()) {
                        iterator.nextTransition();

                        final int countA = iterator.getCount(0);
                        final int countB = iterator.getCount(1);
                        differenceSeen += Math.max(countA, countB) - Math.min(countB, countA);

                    }
                    iterator.close();
                }
            }
            System.out.printf("difference: %d offset=%d %n", differenceSeen, offset);
            if (differenceSeen < minDifferenceSeen) {
                minDifferenceSeen = differenceSeen;
                optimalOffset = offset;
            }
        }
        return optimalOffset;

    }

    private double getSumOfCounts(final CountsReader countReader) throws IOException {
        double sum = 0;
        while (countReader.hasNextTransition()) {
            countReader.nextTransition();
            sum += countReader.getCount();
        }
        countReader.close();
        return sum;
    }

    private boolean isReferenceIncluded(final String referenceSequenceId) {
        return this.includeReferenceNames.size() == 0 || this.includeReferenceNames.contains(referenceSequenceId);
    }

    private void writeHeader(final PrintStream output, final int windowSize) {
        final MutableString buffer = new MutableString();
        buffer.append("position");
        buffer.append('\t');
        buffer.append("referenceId");
        buffer.append('\t');
        for (int i = -windowSize; i < windowSize; i++) {
            buffer.append("position").append(i < 0 ? "" : '+').append(String.valueOf(i));
            buffer.append('\t');
        }
        buffer.append("foldChange");
        buffer.append('\t');
        buffer.append("countA");
        buffer.append('\t');
        buffer.append("countB");
        buffer.append('\t');
        buffer.append("sumA");
        buffer.append('\t');
        buffer.append("sumB");
        buffer.append('\n');
        output.print(buffer);
    }

    private void tallyPosition(final RandomAccessSequenceCache cache,
                               final int referenceIndex, final int position, final double foldChange, final int windowSize,
                               final int referenceSize, final String referenceId,
                               final PrintStream out, final int countA, final int countB, final double sumA, final double sumB) {
        int minPosition = position - windowSize;
        int maxPosition = position + windowSize;
        if (minPosition < 0) {
            minPosition = 0;
        }
        if (maxPosition >= referenceSize) {
            maxPosition = referenceSize - 2;
        }
        final MutableString buffer = new MutableString();
        buffer.append(position);
        buffer.append('\t');
        buffer.append(referenceId);
        buffer.append('\t');
        for (int i = position - windowSize; i < position + windowSize; i++) {
            final char base;
            if (i >= minPosition && i <= maxPosition) {
                base = cache.get(referenceIndex, i);
            } else {
                base = '-';
            }

            buffer.append(base);
            buffer.append('\t');

        }
        buffer.append(foldChange);
        buffer.append('\t');
        buffer.append(countA);
        buffer.append('\t');
        buffer.append(countB);
        buffer.append('\t');
        buffer.append(sumA);
        buffer.append('\t');
        buffer.append(sumB);
        buffer.append('\n');
        out.print(buffer);
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
        new TallyBasesMode().configure(args).execute();
    }
}
