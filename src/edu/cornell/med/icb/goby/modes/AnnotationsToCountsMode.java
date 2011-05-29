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
import edu.cornell.med.icb.goby.counts.CountsArchiveWriter;
import edu.cornell.med.icb.goby.counts.CountsWriter;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * Write annotations corresponding to consensus peaks found in each sequence of count archives.
 */
public class AnnotationsToCountsMode extends AbstractGobyMode {
    /**
     * The mode name.
     */
    private static final String MODE_NAME = "annotations-to-counts";

    /**
     * The mode description help text.
     */
    private static final String MODE_DESCRIPTION =
            "Convert annotations to counts. Annotations can be used to represent segment of the reference sequences. " +
                    "This is can used to represent boundaries of exons on the chromosome sequences or the boundary of any " +
                    "other chromosomal region. This mode converts annotation files to the count format, setting count" +
                    "within the regions defined by the annotations to 1 and zero outside of the annotated regions.";
    /**
     * filename of the annotation file.
     */
    private String inputFilename;
    /**
     * Basename of the output counts archive.
     */
    private String countBasename;
    private int flankingSize = 0;
    private boolean verbose;
    /**
     * The total number of transitions written over the all  targets.
     */
    private int numTransitions;


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
     * @throws java.io.IOException error parsing
     * @throws com.martiansoftware.jsap.JSAPException
     *                             error parsing
     */
    @Override
    public AbstractCommandLineMode configure(final String[] args) throws IOException, JSAPException {
        final JSAPResult jsapResult = parseJsapArguments(args);

        inputFilename = jsapResult.getString("input");
        countBasename = jsapResult.getString("output");
        verbose = jsapResult.getBoolean("verbose");
        flankingSize = jsapResult.getInt("flanking-size");
        return this;
    }

    class AnnotationSegment {
        String chromosome;
        int start;
        int end;

        @Override
        public String toString() {
            return String.format("%s %d-%d ", chromosome, start, end);
        }

        public boolean canCombine(AnnotationSegment last) {
            boolean result = (last.end >= start && chromosome.equals(last.chromosome));
            // System.out.println("canCombine: "+result);
            //   if (result==true) {
            //      System.out.println("last: "+last +" this: "+this);
            //  }
            return result;
        }

        /**
         * Extend this annotation with the boundaries of last.
         *
         * @param last
         */
        public void extendWith(AnnotationSegment last) {
            if (verbose) {
                System.out.printf("joining two annotation segments that overlap first: %s second: %s%n", last, this);
            }
            //     System.out.printf("before: other: %s this: %s",last,this);
            start = Math.min(start, last.start);
            end = Math.max(end, last.end);
            //    System.out.printf("now: other: %s this: %s",last,this);
        }
    }

    /**
     * Run the mode.
     *
     * @throws java.io.IOException error reading / writing
     */
    @Override
    public void execute() {
        CountsArchiveWriter writer = null;
        try {
            writer = new CountsArchiveWriter(countBasename);
        } catch (IOException e) {
            System.err.println("Cannot write counts to basename " + countBasename);
            e.printStackTrace();
        }
        try {
            TSVReader reader = new TSVReader(new FileReader(inputFilename));


            String previousChromosome = "k21o-k2w";
            CountsWriter countWriter = null;
            int lastPosition = 0;
            int refIndex = 0;
            reader.setCommentPrefix("#");
            ObjectArrayList<AnnotationSegment> buffer = new ObjectArrayList<AnnotationSegment>(10);

            while (reader.hasNext()) {
                if (reader.isCommentLine()) {
                    reader.skip();
                } else {
                    reader.next();

                    String chromosome = reader.getString();
                    int startPosition = reader.getInt();
                    int endPosition = reader.getInt();
                    AnnotationSegment a = new AnnotationSegment();
                    a.chromosome = chromosome;
                    a.start = startPosition - flankingSize;
                    a.end = endPosition + flankingSize;
                    if (!buffer.isEmpty()) {
                        AnnotationSegment last = buffer.get(buffer.size() - 1);
                        if (a.canCombine(last)) {
                            last.extendWith(a);
                        } else {
                            buffer.add(a);
                        }
                    } else {
                        buffer.add(a);
                    }
                }
            }

            for (int index = 0; index < buffer.size(); index++) {
                AnnotationSegment a = buffer.get(index);
                //   System.out.printf("size: %d %s %d-%d %n",buffer.size(), a.chromosome, a.start, a.end);
                if (!a.chromosome.equals(previousChromosome)) {

                    if (countWriter != null) {
                        numTransitions += countWriter.getNumberOfTransitions();
                        writer.returnWriter(countWriter);
                        System.out.println("finished writting " + previousChromosome);
                    }
                    countWriter = writer.newCountWriter(refIndex++, a.chromosome);
                    lastPosition = 0;
                    previousChromosome = a.chromosome;
                }

                if (countWriter != null) {
                    if (a.end - a.start <= 0) {
                        System.out.printf("error: annotation has zero or negative length: %d %d %d%n",
                                refIndex, a.start, a.end);
                    } else {
                        //      System.out.println("appending 0 for "+(startPosition - lastPosition));
                        //     System.out.println("appending 1 for "+(endPosition - startPosition));

                        final int length = a.start - lastPosition;
                        if (length > 0) {
                           // System.out.printf("appending 0 for length=%d%n", length);

                            countWriter.appendCount(0, length);
                         //   System.out.printf("appending 1 for length=%d%n", a.end - a.start);
                            countWriter.appendCount(1, a.end - a.start);
                            lastPosition = a.end;
                        }
                    }
                }
            }

            if (countWriter != null) {
                countWriter.appendCount(0,1);
                countWriter.close();
                numTransitions += countWriter.getNumberOfTransitions();
                writer.returnWriter(countWriter);
            }

            writer.close();
            System.out.println("Overall number of transitions written: " + numTransitions);
        }

        catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        catch (IOException e) {
            e.printStackTrace();
        }

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
        new AnnotationsToCountsMode().configure(args).execute();
    }
}