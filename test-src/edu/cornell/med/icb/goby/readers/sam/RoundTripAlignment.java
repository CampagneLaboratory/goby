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

package edu.cornell.med.icb.goby.readers.sam;

/**
 * Class to help with testing round trip alignment comparisons.
 */

import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.modes.CompactToSAMMode;
import edu.cornell.med.icb.goby.modes.SAMToCompactMode;
import edu.cornell.med.icb.goby.modes.SortMode;
import edu.cornell.med.icb.goby.reads.DualRandomAccessSequenceCache;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceInterface;
import edu.cornell.med.icb.goby.reads.RandomAccessSequenceTestSupport;
import it.unimi.dsi.logging.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

/**
 * Configuration for round trip comparison.
 */
public class RoundTripAlignment {
    private static final Logger LOG = Logger.getLogger(RoundTripAlignment.class);

    // make true to quickly debug code (the read sequences will be made of Ns only)
    boolean FAST_GENOME = false;

    String inputGenomeFilename;
    String sourceBamFilename;
    String destGobyBasename;
    String destBamFilename;
    boolean convertBamToGoby = true;
    boolean convertGobyToBam = true;
    boolean keepQualityScores = true;
    boolean keepSoftClips = true;
    boolean stopAtOneFailure = false;
    boolean canonicalMdzForComparison = true;
    /* If the compact-alignment should be sorted after creation. */
    boolean sortGoby = false;
    /* If the sam/bam file should should be sorted after creation. */
    boolean sortSam = false;
    /* If read names should be written to the compact-alignment file */
    boolean preserveReadNames = true;
    /* Perform the record comparisons. */
    boolean compare = true;
    /* To dump differences in read names, positions. */
    boolean dumpReadNames = false;
    /* For spliced reads, this should be set to false. */
    boolean compareWithGoby = true;
    /* If a comparison is made at each source position change. If set to false, all alignments will
     * completely be loaded into memory before comparisons are made. */
    boolean compareAtEachNewPosition = true;
    /* If the source base can be and and thd destination base be a non-N. */
    boolean allowSourceNs = false;

    public void testRoundTripAny() throws IOException {

        // IMPORTANT!!
        // ** These two should always be set to true unless you are doing MANUAL testing and want to not do
        // ** one or both of the conversions.

        if (destGobyBasename != null && !new File(destGobyBasename + ".entries").exists()) {
            convertBamToGoby = true;
        }
        if (!new File(destBamFilename).exists()) {
            convertGobyToBam = true;
        }

        RandomAccessSequenceInterface genome = null;
        if (convertGobyToBam || convertBamToGoby) {
            assertNotNull("Could not locate random-access-genome in specified locations", inputGenomeFilename);
            genome = new DualRandomAccessSequenceCache();
            try {
                System.out.print("Loading random access genome...");
                if (FAST_GENOME) {
                    genome = new RandomAccessSequenceTestSupport(new String[]{"NNNNN"}) {
                        @Override
                        public int getReferenceIndex(final String referenceId) {
                            return 0;
                        }
                    };
                } else {
                    ((DualRandomAccessSequenceCache) genome).load(inputGenomeFilename);
                }
                System.out.println(" done");
            } catch (ClassNotFoundException e) {
                throw new IOException("Could not load genome", e);
            }
        }

        if (convertBamToGoby) {
            LOG.info("Converting bam to compact alignment");
            final SAMToCompactMode importer = new SAMToCompactMode();
            importer.setPreserveReadName(preserveReadNames);
            importer.setInputFile(sourceBamFilename);
            importer.setPreserveSoftClips(keepSoftClips);
            importer.setPreserveAllTags(true);
            importer.setOutputFile(destGobyBasename);
            importer.setGenome(genome);
            importer.setPreserveReadQualityScores(keepQualityScores);
            importer.execute();
            if (sortGoby) {
                LOG.info("Sorting Goby alignment");
                final SortMode sorter = new SortMode();
                sorter.setInput(destGobyBasename);
                sorter.setOutput(destGobyBasename + "-sorted");
                sorter.execute();

                destGobyBasename += "-sorted";
            }
        } else {
            // We didn't just create the goby file, it was created before
            if (sortGoby) {
                // But it was also sorted before, so change the basename
                destGobyBasename += "-sorted";
            }
        }
        final String sortedSamFilename = addBeforeExtension(destBamFilename, "-sorted");
        if (convertGobyToBam) {
            LOG.info("Converting compact alignment to bam");
            final CompactToSAMMode exporter = new CompactToSAMMode();
            exporter.setGenome(genome);
            exporter.setInputBasename(destGobyBasename);
            exporter.setOutput(destBamFilename);
            exporter.execute();

            if (sortSam) {
                final String outputFilename = addBeforeExtension(destBamFilename, "-sorted");
                final String[] samSortArgs = {
                        "I=" + destBamFilename,
                        "O=" + sortedSamFilename,
                        "SO=coordinate",
                        "VALIDATION_STRINGENCY=SILENT"};
                SortSamNoExit.main(samSortArgs);
                destBamFilename = sortedSamFilename;
            }
        } else {
            if (sortSam) {
                destBamFilename = sortedSamFilename;
            }
        }

        LOG.info("Comparing source bam and destination bam");
        final SAMFileReader sourceParser = new SAMFileReader(new FileInputStream(sourceBamFilename));
        final SAMFileReader destParser = new SAMFileReader(new FileInputStream(destBamFilename));
        // We need to set the validation to silent because an incomplete file (if the source isn't the entire file)
        // we can see errors that wouldn't exist in a real conversion.
        sourceParser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        destParser.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        final SAMRecordIterator sourceIterator = sourceParser.iterator();
        final SAMRecordIterator destIterator = destParser.iterator();
        AlignmentReaderImpl gobyReader = null;
        if (compareWithGoby && destGobyBasename != null) {
            gobyReader = new AlignmentReaderImpl(destGobyBasename);
        }
        final ProgressLogger progress = new ProgressLogger(LOG);
        progress.displayFreeMemory = true;
        final SamPerPositionComparison samComparison = new SamPerPositionComparison();
        progress.start();
        samComparison.setMappedQualitiesPreserved(keepQualityScores);
        samComparison.setSoftClipsPreserved(keepSoftClips);
        samComparison.setReadNamesPreserved(preserveReadNames);
        samComparison.setCheckMate(false);
        samComparison.setCanonicalMdzForComparison(canonicalMdzForComparison);
        samComparison.setCompareAtEachNewPosition(compareAtEachNewPosition);
        samComparison.setAllowSourceNs(allowSourceNs);
        int i = 0;
        while (sourceIterator.hasNext()) {
            final SAMRecord expected = sourceIterator.next();
            if (expected.getReadUnmappedFlag()) {
                // We don't store unmapped reads, so skip this source record
                continue;
            }

            SAMRecord actual;
            while (true) {
                assertTrue("Not enough records in destination SAM/BAM file", destIterator.hasNext());
                actual = destIterator.next();
                if (!actual.getReadUnmappedFlag()) {
                    break;
                }
            }

            Alignments.AlignmentEntry gobyActual = null;
            if (gobyReader != null) {
                assertTrue("Not enough records in goby compact-alignment file", gobyReader.hasNext());
                gobyActual = gobyReader.next();
            }


            if (compare) {
                final int compared = samComparison.compare(expected, actual, gobyActual);
                if (stopAtOneFailure && compared > 0) {
                    fail("Comparison failed");
                }
            }
            if (dumpReadNames) {
                boolean performDump = true;
                /*
                if (!expected.getReadName().equals(actual.getReadName())) {
                    performDump = true;
                }
                */
                if (expected.getAlignmentStart() != actual.getAlignmentStart()) {
                    performDump = true;
                }
                if (gobyActual != null) {
                    /*
                    if (!expected.getReadName().equals(gobyActual.getReadName())) {
                        performDump = true;
                    }
                    */
                    if (expected.getAlignmentStart() != gobyActual.getPosition() + 1) {
                        performDump = true;
                    }
                }
                if (performDump) {
                    System.out.println("--" + i + "----------------------------");
                    System.out.println("e=" + expected.getReadName() + " / " + expected.getAlignmentStart());
                    System.out.println("a=" + actual.getReadName() + " / " + actual.getAlignmentStart());
                    if (gobyActual != null) {
                        System.out.println("g=" + gobyActual.getReadName() + " / " + (gobyActual.getPosition() + 1));
                    }
                }
            }
            progress.lightUpdate();
            i++;
        }
        samComparison.finished();
        progress.stop();
        if (!stopAtOneFailure && samComparison.getComparisonFailureCount() > 0) {
            System.out.println("Number of records processed: " + samComparison.getReadNum());
            fail("Number of comparison failures: " + samComparison.getComparisonFailureCount());
        } else {
            System.out.println("Number of records processed: " + samComparison.getReadNum());
            System.out.println("Number of comparison failures: " + samComparison.getComparisonFailureCount());
        }
    }

    private String addBeforeExtension(final String filename, final String toAdd) {
        final String fullPath = FilenameUtils.getFullPath(filename);
        final String basename = FilenameUtils.getBaseName(filename);
        final String ext = FilenameUtils.getExtension(filename);
        return fullPath + basename + toAdd + ext;
    }

    public void createCompactFromSam() throws IOException {

        assertNotNull("Could not locate random-access-genome in specified locations", inputGenomeFilename);
        final RandomAccessSequenceInterface genome = new DualRandomAccessSequenceCache();
        try {
            System.out.print("Loading random access genome...");
            ((DualRandomAccessSequenceCache) genome).load(inputGenomeFilename);
            System.out.println(" done");
        } catch (ClassNotFoundException e) {
            throw new IOException("Could not load genome", e);
        }

        LOG.info("Converting bam to compact alignment");
        final SAMToCompactMode importer = new SAMToCompactMode();
        importer.setInputFile(sourceBamFilename);
        importer.setPreserveSoftClips(keepSoftClips);
        importer.setPreserveAllTags(true);
        importer.setOutputFile(destGobyBasename);
        importer.setGenome(genome);
        importer.setPreserveReadQualityScores(keepQualityScores);
        importer.execute();
    }

    public static byte[] byteArray(final int... bytes) {
        final byte[] result = new byte[bytes.length];
        for (int i = 0; i < bytes.length; i++) {
            result[i] = (byte) bytes[i];
        }
        return result;
    }
}
