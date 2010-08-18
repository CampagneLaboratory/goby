/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.reads.ReadsWriter;
import it.unimi.dsi.fastutil.ints.Int2BooleanOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2FloatOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import org.apache.commons.io.FileUtils;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class TestLastToCompact {
    @Test
    public void testLastToCompact1() throws IOException {

        // test LastToCompact convert+filtering

        // read fake data and convert+filter
        final LastToCompactMode processor = new LastToCompactMode();
        final int LAST_TO_COMPACT_M_PARAM = 2;
        processor.setAmbiguityThreshold(LAST_TO_COMPACT_M_PARAM);
        processor.setInputFile("test-results/alignments/last-to-compact/last-101.maf");
        processor.setOutputFile("test-results/alignments/last-to-compact/last-101.compact");
        processor.setTargetReferenceIdsFilename("test-results/alignments/last-to-compact/last-reference.compact-reads");
        processor.setOnlyMafFile(true);
        processor.setNumberOfReads(2857819);
        processor.setLargestQueryIndex(2857819);
        processor.setSmallestQueryIndex(0);
        processor.setPropagateQueryIds(false);
        processor.setPropagateTargetIds(true);
        processor.execute();


        // read compact alignment results
        final AlignmentReader reader = new AlignmentReader(processor.getOutputFile());
        reader.readHeader();
        assertEquals(2857819, reader.getNumberOfQueries());
        assertEquals(1, reader.getNumberOfTargets());


        // lookup tables
        final Int2IntOpenHashMap queryIndex2NumberOfHits = new Int2IntOpenHashMap();

        final Int2FloatOpenHashMap queryIndex2Score = new Int2FloatOpenHashMap();
        final Int2IntOpenHashMap queryIndex2Multiplicity = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2NumberOfIndels = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2NumberOfMismatches = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2Position = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2QueryAlignedLength = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2QueryPosition = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2TargetIndex = new Int2IntOpenHashMap();
        final Int2BooleanOpenHashMap queryIndex2MatchingReverseStrand = new Int2BooleanOpenHashMap();
        final Int2IntOpenHashMap queryIndex2TargetAlignedLength = new Int2IntOpenHashMap();

        // enter alignment data
        int qii;
        while (reader.hasNext()) {

            final Alignments.AlignmentEntry aln = reader.next();
            qii = aln.getQueryIndex();

            final int numHits = queryIndex2NumberOfHits.get(qii);

            queryIndex2NumberOfHits.put(qii, numHits + 1);

            queryIndex2Score.put(qii, aln.getScore());
            queryIndex2Multiplicity.put(qii, aln.getMultiplicity());
            queryIndex2NumberOfIndels.put(qii, aln.getNumberOfIndels());
            queryIndex2NumberOfMismatches.put(qii, aln.getNumberOfMismatches());
            queryIndex2Position.put(qii, aln.getPosition());
            queryIndex2QueryAlignedLength.put(qii, aln.getQueryAlignedLength());
            queryIndex2QueryPosition.put(qii, aln.getQueryPosition());
            queryIndex2TargetIndex.put(qii, aln.getTargetIndex());
            queryIndex2MatchingReverseStrand.put(qii, aln.getMatchingReverseStrand());
            queryIndex2TargetAlignedLength.put(qii, aln.getTargetAlignedLength());

        }

        //
        // validate alignment data using values from getMafInput() below
        //

        // there are a total of 5 entries with ID 2857818
        // 2 entries have score = 35 (the maximum score for this ID)
        // 3 entries are filtered b/c their scores are below 35
        qii = 2857818;
        assertEquals(queryIndex2NumberOfHits.get(qii), 2);
        assertEquals((int) queryIndex2Score.get(qii), 35);
        assertEquals(queryIndex2Multiplicity.get(qii), 1);
        assertEquals(queryIndex2NumberOfIndels.get(qii), 0);
        assertEquals(queryIndex2NumberOfMismatches.get(qii), 0);
        assertEquals(queryIndex2Position.get(qii), 1614); // last entry added
        assertEquals(queryIndex2QueryAlignedLength.get(qii), 35);
        assertEquals(queryIndex2QueryPosition.get(qii), 0);
        assertEquals(queryIndex2TargetIndex.get(qii), 0);
        assertEquals(queryIndex2MatchingReverseStrand.get(qii), false);
        assertEquals(queryIndex2TargetAlignedLength.get(qii), 35);

        // there are 5 entries with the score = 35 (the maximum score for this ID)
        // filtered due to ambiguity
        qii = 577287;
        assertEquals(queryIndex2NumberOfHits.get(qii), 0);
        assertEquals((int) queryIndex2Score.get(qii), 0);
        assertEquals(queryIndex2Multiplicity.get(qii), 0);
        assertEquals(queryIndex2NumberOfIndels.get(qii), 0);
        assertEquals(queryIndex2NumberOfMismatches.get(qii), 0);
        assertEquals(queryIndex2Position.get(qii), 0);
        assertEquals(queryIndex2QueryAlignedLength.get(qii), 0);
        assertEquals(queryIndex2QueryPosition.get(qii), 0);
        assertEquals(queryIndex2TargetIndex.get(qii), 0);
        assertEquals(queryIndex2MatchingReverseStrand.get(qii), false);
        assertEquals(queryIndex2TargetAlignedLength.get(qii), 0);

        //
        reader.close();
    }

    @Test
    public void testLastToCompact2() throws IOException {

        // test LastToCompact convert+filtering

        // read fake data and convert+filter
        final LastToCompactMode processor = new LastToCompactMode();
        final int LAST_TO_COMPACT_M_PARAM = 2;
        processor.setAmbiguityThreshold(LAST_TO_COMPACT_M_PARAM);
        processor.setInputFile("test-results/alignments/last-to-compact/last-102.maf");
        processor.setOutputFile("test-results/alignments/last-to-compact/last-102.compact");
        processor.setTargetReferenceIdsFilename("test-results/alignments/last-to-compact/last-reference.compact-reads");
        processor.setOnlyMafFile(true);
        processor.setNumberOfReads(3538282);
        processor.setLargestQueryIndex(3538282);
        processor.setSmallestQueryIndex(0);
        processor.setPropagateQueryIds(false);
        processor.setPropagateTargetIds(true);
        processor.execute();


        // read compact alignment results
        final AlignmentReader reader = new AlignmentReader(processor.getOutputFile());
        reader.readHeader();
        assertEquals(3538282, reader.getNumberOfQueries());
        assertEquals(1, reader.getNumberOfTargets());


        // lookup tables
        final Int2IntOpenHashMap queryIndex2NumberOfHits = new Int2IntOpenHashMap();

        final Int2FloatOpenHashMap queryIndex2Score = new Int2FloatOpenHashMap();
        final Int2IntOpenHashMap queryIndex2Multiplicity = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2NumberOfIndels = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2NumberOfMismatches = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2Position = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2QueryAlignedLength = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2QueryPosition = new Int2IntOpenHashMap();
        final Int2IntOpenHashMap queryIndex2TargetIndex = new Int2IntOpenHashMap();
        final Int2BooleanOpenHashMap queryIndex2MatchingReverseStrand = new Int2BooleanOpenHashMap();
        final Int2IntOpenHashMap queryIndex2TargetAlignedLength = new Int2IntOpenHashMap();

        // enter alignment data
        int qii;
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry aln = reader.next();
            qii = aln.getQueryIndex();

            final int numHits = queryIndex2NumberOfHits.get(qii);

            queryIndex2NumberOfHits.put(qii, numHits + 1);

            queryIndex2Score.put(qii, aln.getScore());
            queryIndex2Multiplicity.put(qii, aln.getMultiplicity());
            queryIndex2NumberOfIndels.put(qii, aln.getNumberOfIndels());
            queryIndex2NumberOfMismatches.put(qii, aln.getNumberOfMismatches());
            queryIndex2Position.put(qii, aln.getPosition());
            queryIndex2QueryAlignedLength.put(qii, aln.getQueryAlignedLength());
            queryIndex2QueryPosition.put(qii, aln.getQueryPosition());
            queryIndex2TargetIndex.put(qii, aln.getTargetIndex());
            queryIndex2MatchingReverseStrand.put(qii, aln.getMatchingReverseStrand());
            queryIndex2TargetAlignedLength.put(qii, aln.getTargetAlignedLength());

        }

        //
        // validate alignment data using values from getMafInput() below
        //

        // there are a total of 5 entries with ID 2857818
        // 2 entries have score = 35 (the maximum score for this ID)
        // 3 entries are filtered b/c their scores are below 35
        qii = 2857818;
        assertEquals(queryIndex2NumberOfHits.get(qii), 2);
        assertEquals((int) queryIndex2Score.get(qii), 35);
        assertEquals(queryIndex2Multiplicity.get(qii), 1);
        assertEquals(queryIndex2NumberOfIndels.get(qii), 0);
        assertEquals(queryIndex2NumberOfMismatches.get(qii), 0);
        assertEquals(queryIndex2Position.get(qii), 1614);  // last entry added
        assertEquals(queryIndex2QueryAlignedLength.get(qii), 35);
        assertEquals(queryIndex2QueryPosition.get(qii), 0);
        assertEquals(queryIndex2TargetIndex.get(qii), 0);
        assertEquals(queryIndex2MatchingReverseStrand.get(qii), false);
        assertEquals(queryIndex2TargetAlignedLength.get(qii), 35);

        // there are 5 entries with the score = 35 (the maximum score for this ID)
        // filtered due to ambiguity
        qii = 577287;
        assertEquals(queryIndex2NumberOfHits.get(qii), 0);
        assertEquals((int) queryIndex2Score.get(qii), 0);
        assertEquals(queryIndex2Multiplicity.get(qii), 0);
        assertEquals(queryIndex2NumberOfIndels.get(qii), 0);
        assertEquals(queryIndex2NumberOfMismatches.get(qii), 0);
        assertEquals(queryIndex2Position.get(qii), 0);
        assertEquals(queryIndex2QueryAlignedLength.get(qii), 0);
        assertEquals(queryIndex2QueryPosition.get(qii), 0);
        assertEquals(queryIndex2TargetIndex.get(qii), 0);
        assertEquals(queryIndex2MatchingReverseStrand.get(qii), false);
        assertEquals(queryIndex2TargetAlignedLength.get(qii), 0);


        // check more values ...

        qii = 1431626;
        assertEquals(queryIndex2NumberOfHits.get(qii), 1);
        assertEquals((int) queryIndex2Score.get(qii), 35);
        assertEquals(queryIndex2Multiplicity.get(qii), 1);
        assertEquals(queryIndex2NumberOfIndels.get(qii), 0);
        assertEquals(queryIndex2NumberOfMismatches.get(qii), 0);
        assertEquals(queryIndex2Position.get(qii), 1916);
        assertEquals(queryIndex2QueryAlignedLength.get(qii), 35);
        assertEquals(queryIndex2QueryPosition.get(qii), 0);
        assertEquals(queryIndex2TargetIndex.get(qii), 0);
        assertEquals(queryIndex2MatchingReverseStrand.get(qii), false);
        assertEquals(queryIndex2TargetAlignedLength.get(qii), 35);

        qii = 1135141;
        assertEquals(queryIndex2NumberOfHits.get(qii), 1);
        assertEquals((int) queryIndex2Score.get(qii), 35);
        assertEquals(queryIndex2Multiplicity.get(qii), 1);
        assertEquals(queryIndex2NumberOfIndels.get(qii), 0);
        assertEquals(queryIndex2NumberOfMismatches.get(qii), 0);
        assertEquals(queryIndex2Position.get(qii), 1995);
        assertEquals(queryIndex2QueryAlignedLength.get(qii), 35);
        assertEquals(queryIndex2QueryPosition.get(qii), 0);
        assertEquals(queryIndex2TargetIndex.get(qii), 0);
        assertEquals(queryIndex2MatchingReverseStrand.get(qii), false);
        assertEquals(queryIndex2TargetAlignedLength.get(qii), 35);

        /*
        qii = 2527313;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),      1       );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),     35       );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),      1       );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 2142468;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),      1       );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),     35       );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),      1       );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 2492089;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 3392286;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii =   24734;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii =  301127;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 2771094;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 1750522;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 1434206;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii = 3538281;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );

        qii =  751638;
        assertEquals(      queryIndex2NumberOfHits          .get ( qii  ),              );
        assertEquals((int) queryIndex2Score                 .get ( qii  ),              );
        assertEquals(      queryIndex2Multiplicity          .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfIndels        .get ( qii  ),              );
        assertEquals(      queryIndex2NumberOfMismatches    .get ( qii  ),              );
        assertEquals(      queryIndex2Position              .get ( qii  ),              );
        assertEquals(      queryIndex2QueryAlignedLength    .get ( qii  ),              );
        assertEquals(      queryIndex2QueryPosition         .get ( qii  ),              );
        assertEquals(      queryIndex2TargetIndex           .get ( qii  ),              );
        assertEquals(      queryIndex2MatchingReverseStrand .get ( qii  ),  false       );
        assertEquals(      queryIndex2TargetAlignedLength   .get ( qii  ),              );
        */

        //
        reader.close();
    }


    @Before
    public void setUp() throws IOException {
        final File dir = new File("test-results/alignments/last-to-compact/");
        dir.mkdirs();

        final FileWriter writer1 = new FileWriter("test-results/alignments/last-to-compact/last-101.maf");
        writer1.write(getMafInput1());
        writer1.close();

        final FileWriter writer2 = new FileWriter("test-results/alignments/last-to-compact/last-102.maf");
        writer2.write(getMafInput2());
        writer2.close();

        final FileWriter writer3 = new FileWriter("test-results/alignments/last-to-compact/last-103-variations.maf");
        writer3.write(getMafInput3Variations());
        writer3.close();

        final ReadsWriter referenceWriter = new ReadsWriter(FileUtils.openOutputStream(
                new File("test-results/alignments/last-to-compact/last-reference.compact-reads")));
        referenceWriter.setIdentifier("0");
        referenceWriter.appendEntry();
        referenceWriter.close();
    }

    private String getMafInput1() {
        return "# LAST version 18959\n" +
                "#\n" +
                "# a=2 b=1 c=100000 e=31 d=18 x=22 y=10\n" +
                "# u=0 s=2 m=10 l=1 k=1 i=134217728 w=1000 t=-1 g=1 j=3 z=0\n" +
                "# /scratchLocal/hadoop/hadoop-hadoop/mapred/local/taskTracker/jobcache/job_200905080202_0031/attempt_200905080202_0031_m_000001_0/work/tmp/Homo_sapiens.NCBI3\n" +
                "6.52.dna.chromosome.1.db\n" +
                "#\n" +
                "#    A  C  G  T\n" +
                "# A  1 -1 -1 -1\n" +
                "# C -1  1 -1 -1\n" +
                "# G -1 -1  1 -1\n" +
                "# T -1 -1 -1  1\n" +
                "#\n" +
                "# Coordinates are 0-based.  For - strand matches, coordinates\n" +
                "# in the reverse complement of the 2nd sequence are used.\n" +
                "#\n" +
                "# name start alnSize strand seqSize alignment\n" +
                "#\n" +
                "a score=35\n" +
                "s 0       1514 35 + 247249719 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=35\n" +
                "s 0       1614 35 + 247249719 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1714 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1814 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1914 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1598 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1698 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1798 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1898 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1998 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "";
    }

    private String getMafInput2() {
        return "# LAST version 18959\n" +
                "#\n" +
                "# a=2 b=1 c=100000 e=31 d=18 x=22 y=10\n" +
                "# u=0 s=2 m=10 l=1 k=1 i=134217728 w=1000 t=-1 g=1 j=3 z=0\n" +
                "# /scratchLocal/hadoop/hadoop-hadoop/mapred/local/taskTracker/jobcache/job_200905080202_0031/attempt_200905080202_0031_m_000001_0/work/tmp/Homo_sapiens.NCBI3\n" +
                "6.52.dna.chromosome.1.db\n" +
                "#\n" +
                "#    A  C  G  T\n" +
                "# A  1 -1 -1 -1\n" +
                "# C -1  1 -1 -1\n" +
                "# G -1 -1  1 -1\n" +
                "# T -1 -1 -1  1\n" +
                "#\n" +
                "# Coordinates are 0-based.  For - strand matches, coordinates\n" +
                "# in the reverse complement of the 2nd sequence are used.\n" +
                "#\n" +
                "# name start alnSize strand seqSize alignment\n" +
                "#\n" +
                "a score=35\n" +
                "s 0       1514 35 + 247249719 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=35\n" +
                "s 0       1614 35 + 247249719 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1714 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1814 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1914 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1598 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1698 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1798 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1898 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1998 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0       1916 35 + 247249719 TGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAG\n" +
                "s 1431626    0 35 +        35 TGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAG\n" +
                "\n" +
                "a score=35\n" +
                "s 0       1995 35 + 247249719 GTTGTCTGCATGTAACTTAATACCACAACCAGGCA\n" +
                "s 1135141    0 35 +        35 GTTGTCTGCATGTAACTTAATACCACAACCAGGCA\n" +
                "\n" +
                "a score=35\n" +
                "s 0       2389 35 + 247249719 CACAGCCAAGAGGGCTGAAGAAATGGTAGAACGGA\n" +
                "s 2527313    0 35 +        35 CACAGCCAAGAGGGCTGAAGAAATGGTAGAACGGA\n" +
                "\n" +
                "a score=35\n" +
                "s 0       2539 35 + 247249719 GCCGACTTGGATCACACTCTTGTGAGTGTCCCCAG\n" +
                "s 2142468    0 35 +        35 GCCGACTTGGATCACACTCTTGTGAGTGTCCCCAG\n" +
                "\n" +
                "a score=35\n" +
                "s 0       2544 35 + 247249719 CTTGGATCACACTCTTGTGAGTGTCCCCAGTGTTG\n" +
                "s 2492089    0 35 +        35 CTTGGATCACACTCTTGTGAGTGTCCCCAGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0       2579 35 + 247249719 CAGAGGTGAGAGGAGAGTAGACAGTGAGTGGGAGT\n" +
                "s 3392286    0 35 +        35 CAGAGGTGAGAGGAGAGTAGACAGTGAGTGGGAGT\n" +
                "\n" +
                "a score=35\n" +
                "s 0     2598 35 + 247249719 GACAGTGAGTGGGAGTGGCGTCGCCCCTAGGGCTC\n" +
                "s 24734    0 35 +        35 GACAGTGAGTGGGAGTGGCGTCGCCCCTAGGGCTC\n" +
                "\n" +
                "a score=35\n" +
                "s 0      3108 35 + 247249719 TTCTGCCAGCATAGTGCTCCTGGACCAGTGATACA\n" +
                "s 301127    0 35 +        35 TTCTGCCAGCATAGTGCTCCTGGACCAGTGATACA\n" +
                "\n" +
                "a score=35\n" +
                "s 0       3242 35 + 247249719 CTCCTGCCTTTTCCTTTCCCTAGAGCCTCCACCAC\n" +
                "s 2771094    0 35 +        35 CTCCTGCCTTTTCCTTTCCCTAGAGCCTCCACCAC\n" +
                "\n" +
                "a score=35\n" +
                "s 0       3858 35 + 247249719 CCCTAACCTGCCCCACAGCCTTGCCTGGATTTCTA\n" +
                "s 1750522    0 35 +        35 CCCTAACCTGCCCCACAGCCTTGCCTGGATTTCTA\n" +
                "\n" +
                "a score=35\n" +
                "s 0       4269 35 + 247249719 CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTC\n" +
                "s 1434206    0 35 +        35 CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTC\n" +
                "\n" +
                "a score=35\n" +
                "s 0       4274 35 + 247249719 CAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGG\n" +
                "s 3538281    0 35 +        35 CAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      4295 35 + 247249719 GCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTG\n" +
                "s 751638    0 35 +        35 GCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTG\n" +
                "";
    }


    /**
     * Test the parsing of sequence variations from Last MAF file.
     *
     * @throws IOException
     */
    @Test
    public void testLastToCompactVariations() throws IOException {

        // test LastToCompact convert+filtering

        // read fake data and convert+filter
        final LastToCompactMode processor = new LastToCompactMode();
        final int LAST_TO_COMPACT_M_PARAM = 2;
        processor.setAmbiguityThreshold(LAST_TO_COMPACT_M_PARAM);
        processor.setInputFile("test-results/alignments/last-to-compact/last-103-variations.maf");
        processor.setOutputFile("test-results/alignments/last-to-compact/last-103-variations");
        processor.setTargetReferenceIdsFilename("test-results/alignments/last-to-compact/last-reference.compact-reads");
        processor.setOnlyMafFile(true);
        processor.setNumberOfReads(2857822);
        processor.setLargestQueryIndex(2857822);
        processor.setSmallestQueryIndex(0);
        processor.setPropagateQueryIds(false);
        processor.setPropagateTargetIds(true);
        processor.setQualityFilterParameters("threshold=1"); // allow everything to pass through, important to detect all variations in this test.
        processor.execute();


        // read compact alignment results
        final AlignmentReader reader = new AlignmentReader(processor.getOutputFile());
        reader.readHeader();
        assertEquals(2857822, reader.getNumberOfQueries());
        assertEquals(1, reader.getNumberOfTargets());
        int entryIndex = 0;
        while (reader.hasNext()) {
            final Alignments.AlignmentEntry alignmentEntry = reader.next();

            System.out.println(entryIndex + " entry : " + alignmentEntry);
            switch (entryIndex) {
                case 0:
                    assertEquals(2, alignmentEntry.getSequenceVariationsCount());
                    final Alignments.SequenceVariation var1 = alignmentEntry.getSequenceVariations(0);
                    assertEquals("C", var1.getTo());
                    assertEquals("G", var1.getFrom());
                    assertEquals(1, var1.getPosition());

                    final Alignments.SequenceVariation var2 = alignmentEntry.getSequenceVariations(1);
                    assertEquals("C", var2.getTo());
                    assertEquals("A", var2.getFrom());
                    assertEquals(11, var2.getPosition());
                    break;

                case 1:
                    assertEquals(1, alignmentEntry.getSequenceVariationsCount());
                    final Alignments.SequenceVariation var1_0 = alignmentEntry.getSequenceVariations(0);
                    assertEquals("---", var1_0.getTo());
                    assertEquals("TTT", var1_0.getFrom());
                    assertEquals(3, var1_0.getPosition());
                    break;

                case 2:
                    assertEquals(2, alignmentEntry.getSequenceVariationsCount());
                    final Alignments.SequenceVariation var2_1 = alignmentEntry.getSequenceVariations(0);
                    assertEquals("A", var2_1.getTo());
                    assertEquals("G", var2_1.getFrom());
                    assertEquals(2, var2_1.getPosition());

                    final Alignments.SequenceVariation var2_2 = alignmentEntry.getSequenceVariations(1);
                    assertEquals("A", var2_2.getTo());
                    assertEquals("-", var2_2.getFrom());
                    assertEquals(29, var2_2.getPosition());
                    break;
                default:
                    break;
            }
            entryIndex++;
        }
        assertEquals(4, entryIndex);
    }

    private String getMafInput3Variations() {
        return "# LAST version 18959\n" +
                "#\n" +
                "# a=2 b=1 c=100000 e=31 d=18 x=22 y=10\n" +
                "# u=0 s=2 m=10 l=1 k=1 i=134217728 w=1000 t=-1 g=1 j=3 z=0\n" +
                "# /scratchLocal/hadoop/hadoop-hadoop/mapred/local/taskTracker/jobcache/job_200905080202_0031/attempt_200905080202_0031_m_000001_0/work/tmp/Homo_sapiens.NCBI3\n" +
                "6.52.dna.chromosome.1.db\n" +
                "#\n" +
                "#    A  C  G  T\n" +
                "# A  1 -1 -1 -1\n" +
                "# C -1  1 -1 -1\n" +
                "# G -1 -1  1 -1\n" +
                "# T -1 -1 -1  1\n" +
                "#\n" +
                "# Coordinates are 0-based.  For - strand matches, coordinates\n" +
                "# in the reverse complement of the 2nd sequence are used.\n" +
                "#\n" +
                "# name start alnSize strand seqSize alignment\n" +
                "#\n" +
                "a score=35\n" +
                "s 0       1514 35 + 247249719 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 CATTTTTGCCCGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=36\n" +
                "s 0       1614 35 + 247249719 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857819    0 35 +        35 GA---TTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1714 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCC-CTGGAGA\n" +
                "s 2857820    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCACTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1814 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857821    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=34\n" +
                "s 0       1914 35 + 247249719 GGTTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "s 2857818    0 35 +        35 GATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGA\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1598 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1698 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1798 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1898 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "\n" +
                "a score=35\n" +
                "s 0      1998 35 + 247249719 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "s 577287    0 35 +        35 TTTCCACTGATGATTTTGCTGCATGGCCGGTGTTG\n" +
                "";
    }
}
