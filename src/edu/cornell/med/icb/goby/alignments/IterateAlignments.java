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

package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.lang.MutableString;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;

import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.goby.algorithmic.algorithm.ComputeCount;
import edu.cornell.med.icb.goby.algorithmic.algorithm.ComputeStartCount;
import edu.cornell.med.icb.goby.counts.CountsArchiveWriter;
import com.martiansoftware.jsap.JSAPResult;

/**
 * A helper class to iterate through a set of alignments and process only a subset of references in each alignment.
 *
 * @author Fabien Campagne
 *         Date: Mar 10, 2010
 *         Time: 11:12:23 AM
 */
public abstract class IterateAlignments {
    boolean filterByReferenceNames;
    private ObjectSet<String> includeReferenceNames = new ObjectOpenHashSet<String>();
    private IntSet referencesToProcess;
    private DoubleIndexedIdentifier referenceIds;


    /**
     * Iterate through a set of alignments. Iterations are performed through these steps:
     * <UL>
     * <LI>Iterate will call processHeader on each input alignment, giving the opportunity
     * to the client to read the header opf the alignment and extract information from it.
     * </LI>
     * </UL>
     *
     * @param basenames
     * @throws FileNotFoundException
     */
    public void iterate(String[] basenames) throws IOException {


        for (String basename : basenames) {
            final AlignmentReader reader = new AlignmentReader(basename);
            reader.readHeader();
            final int numberOfReferences = reader.getNumberOfTargets();

            referenceIds = new DoubleIndexedIdentifier(reader.getTargetIdentifiers());
            reader.close();
            System.out.println(String.format("Alignment contains %d reference sequences", numberOfReferences));
            processNumberOfReferences(basename, numberOfReferences);
            //  CountsWriter writers[] = new CountsWriter[numberOfReferences];
            referencesToProcess = new IntOpenHashSet();

            // setup referencesToProcess data structure according to the command line (filterByReferenceNames and includeReferenceNames)
            for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {


                final MutableString referenceId = referenceIds.getId(referenceIndex);
                assert referenceId!=null: "reference id cannot be null for reference index="+referenceIndex;
                final String referenceName = referenceId.toString();
                if (filterByReferenceNames) {
                    if (includeReferenceNames.contains(referenceName)) {
                        // subset of reference names selected by the command line:
                        referencesToProcess.add(referenceIndex);
                    }
                } else {
                    // process each sequence:
                    referencesToProcess.add(referenceIndex);
                }
            }
            // Give the client the ability to prepare data structures for each reference that will be processed.
            for (int referenceIndex = 0; referenceIndex < numberOfReferences; referenceIndex++) {
                if (referencesToProcess.contains(referenceIndex)) {
                    prepareDataStructuresForReference(referenceIndex);
                }
            }

            final AlignmentReader alignmentReader = new AlignmentReader(basename);
            alignmentReader.readHeader();

            // read the alignment:
            System.out.println("Loading the alignment " + basename);
            for (final Alignments.AlignmentEntry alignmentEntry : alignmentReader) {
                final int referenceIndex = alignmentEntry.getTargetIndex();
                if (referencesToProcess.contains(referenceIndex)) {
                    processAlignmentEntry(alignmentReader, alignmentEntry);
                }
            }
            reader.close();


        }
    }

    public abstract void processAlignmentEntry(AlignmentReader alignmentReader, Alignments.AlignmentEntry alignmentEntry);

    public void prepareDataStructuresForReference(int referenceIndex) {

    }

    /**
     * Will be called to let the client know how many references will be processed in a given alignment.
     *
     * @param basename
     * @param numberOfReferences
     * @throws IOException
     */
    public void processNumberOfReferences(String basename, int numberOfReferences) throws IOException {
    }

    /**
     * Parse the string of reference sequences to include the iteration.
     *
     * @param jsapResult The jsapResult available to the mode.
     */
    public void parseIncludeReferenceArgument(final JSAPResult jsapResult) {

        final String includeReferenceNameComas = jsapResult.getString("include-reference-names");
        if (includeReferenceNameComas != null) {
            includeReferenceNames = new ObjectOpenHashSet<String>();
            includeReferenceNames.addAll(Arrays.asList(includeReferenceNameComas.split("[,]")));
            System.out.println("Will write counts for the following sequences:");
            for (final String name : includeReferenceNames) {
                System.out.println(name);
            }
            filterByReferenceNames = true;
        }

    }

    /**
     * Return the reference sequence id given its index.
     * @param targetIndex
     * @return
     */
    protected CharSequence getReferenceId(int targetIndex) {
        return referenceIds.getId(targetIndex);
    }
}
