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

package edu.cornell.med.icb.goby.reads;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;
import edu.cornell.med.icb.goby.modes.CompactToFastaMode;
import edu.cornell.med.icb.parsers.FastaParser;
import it.unimi.dsi.fastutil.chars.Char2IntArrayMap;
import it.unimi.dsi.fastutil.chars.Char2IntMap;
import it.unimi.dsi.fastutil.chars.Char2ObjectArrayMap;
import it.unimi.dsi.fastutil.chars.Char2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2CharArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Date;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

/**
 * Convert a sequence to color space. Bases that are not recognized result in code 7
 * (i.e., digram 'A?' will be encoded as 7).
 *
 * @author Fabien Campagne
 *         Date: May 18, 2009
 *         Time: 5:32:52 PM
 */
public final class ColorSpaceConverter {
    private static final Char2ObjectMap<Char2IntMap> CONVERSION_MAP;
    private static final Int2ObjectMap<Int2CharArrayMap> DECODING_MAP;
    public static final int UNKNOWN = 7;
    public static final char UNKNOWN_BASE = 'N';
     static {
        CONVERSION_MAP = new Char2ObjectArrayMap<Char2IntMap>();
        CONVERSION_MAP.put('A', new Char2IntArrayMap());
        CONVERSION_MAP.put('C', new Char2IntArrayMap());
        CONVERSION_MAP.put('T', new Char2IntArrayMap());
        CONVERSION_MAP.put('G', new Char2IntArrayMap());
        CONVERSION_MAP.put('N', new Char2IntArrayMap());
        for (final Char2IntMap map : CONVERSION_MAP.values()) {
            map.defaultReturnValue(UNKNOWN);
        }
        DECODING_MAP = new Int2ObjectArrayMap<Int2CharArrayMap>();
        DECODING_MAP.put(0, new Int2CharArrayMap());
        DECODING_MAP.put(1, new Int2CharArrayMap());
        DECODING_MAP.put(2, new Int2CharArrayMap());
        DECODING_MAP.put(3, new Int2CharArrayMap());
        DECODING_MAP.put(4, new Int2CharArrayMap());
        DECODING_MAP.put(5, new Int2CharArrayMap());
        DECODING_MAP.put(6, new Int2CharArrayMap());

        for (final Int2CharArrayMap map : DECODING_MAP.values()) {
            map.defaultReturnValue(UNKNOWN_BASE);
        }
        push("AA", 0);
        push("CC", 0);
        push("GG", 0);
        push("TT", 0);
        push("AC", 1);
        push("CA", 1);
        push("GT", 1);
        push("TG", 1);
        push("AG", 2);
        push("CT", 2);
        push("GA", 2);
        push("TC", 2);
        push("AT", 3);
        push("CG", 3);
        push("GC", 3);
        push("TA", 3);
        push("AN", 4);
        push("CN", 4);
        push("GN", 4);
        push("TN", 4);
        push("NA", 5);
        push("NC", 5);
        push("NG", 5);
        push("NT", 5);
        push("NN", 6);



    }
    /**
     * Private constructor for utility class.
     */
    private ColorSpaceConverter() {
        super();
    }

    private static void push(final String s, final int colorCode) {
        final char firstBase = s.charAt(0);
        final char secondBase = s.charAt(1);
        CONVERSION_MAP.get(firstBase).put(secondBase, colorCode);
        DECODING_MAP.get(colorCode).put(firstBase, secondBase);
    }

    /**
     * Return the base corresponding to the previous base and color code transition.
     *
     * @param previousBase Base in sequence space at the previous position.
     * @param colorCode    Transition in color space.
     * @return
     */
    public static char decodeColor(final char previousBase, final char colorCode) {
        return DECODING_MAP.get(colorCode).get(previousBase);
    }

    private static int getColorCode(final char firstBase, final char secondBase) {
        final Char2IntMap intMap = CONVERSION_MAP.get(firstBase);
        if (intMap == null) {
            return UNKNOWN;
        }
        return intMap.get(secondBase);
    }



    /**
     * Converts a sequence into the equivalent sequence in color space.
     *
     * @param input  The sequence to be converted
     * @param output where the converted sequence should be placed
     */
    public static void convert(final CharSequence input, final MutableString output) {
        convert(input, output, false);
    }

    /**
     * Converts a sequence into the equivalent sequence in color space.
     *
     * @param input  The sequence to be converted
     * @param output where the converted sequence should be placed
     * @param anyN   if true, converts color space codes larger or equal to 4 to 'N' characters.
     */
    public static void convert(final CharSequence input, final MutableString output,
                               final boolean anyN) {
        assert output != null : "The output location must not be null";

        output.setLength(0);
        if (input != null) {
            int position = 0;
            final int length = input.length() - 1;    // -1 since we enumerate digrams
            if (input.length() > 0) {

                output.setLength(position + 1);
                output.setCharAt(position++, input.charAt(0));
            }
            for (int index = 0; index < length; ++index) {
                final char code = Character.forDigit(getColorCode(input.charAt(index), input.charAt(index + 1)), 10);
                output.setLength(position + 1);
                output.setCharAt(position++, (code >= '4') ? (anyN ? 'N' : code) : code);
            }
            output.setLength(position);
        }
    }


    public static String getColorSpaceSubstitutionMatrix() {
        return "# Substitution matrix to align in color space. Perfect \n" +
                "# matches score 2, mimatches between valid base bigrams (0-3) score -1.\n" +
                "\n" +
                "      0     1     2     3\n" +
                "0     2   -1    -1    -1\n" +
                "1    -1    2    -1    -1\n" +
                "2    -1   -1     2    -1\n" +
                "3    -1   -1    -1     2\n";
    }

    public static void main(final String[] args) throws JSAPException, IOException {
        final JSAP jsap = new JSAP();

        final FlaggedOption sequenceOption = new FlaggedOption("input");
        sequenceOption.setRequired(true);
        sequenceOption.setLongFlag("input");
        sequenceOption.setShortFlag('i');
        sequenceOption.setStringParser(
                FileStringParser.getParser().setMustBeFile(true).setMustExist(true));
        sequenceOption.setHelp("The input file (in Fasta format) to convert");
        jsap.registerParameter(sequenceOption);

        final FlaggedOption outputOption = new FlaggedOption("output");
        outputOption.setRequired(false);
        outputOption.setLongFlag("output");
        outputOption.setShortFlag('o');
        outputOption.setStringParser(FileStringParser.getParser().setMustBeFile(true));
        outputOption.setHelp("The output file to write to (default = stdout)");
        jsap.registerParameter(outputOption);


        final FlaggedOption titleOption = new FlaggedOption("title");
        titleOption.setRequired(false);
        titleOption.setLongFlag("title");
        titleOption.setShortFlag('t');
        titleOption.setHelp("Title for this conversion");
        jsap.registerParameter(titleOption);

        final Switch verboseOption = new Switch("verbose");
        verboseOption.setLongFlag("verbose");
        verboseOption.setShortFlag('v');
        verboseOption.setHelp("Verbose output");
        jsap.registerParameter(verboseOption);

        final Switch helpOption = new Switch("help");
        helpOption.setLongFlag("help");
        helpOption.setShortFlag('h');
        helpOption.setHelp("Print this message");
        jsap.registerParameter(helpOption);

        jsap.setUsage("Usage: " + ColorSpaceConverter.class.getName() + " " + jsap.getUsage());

        final JSAPResult result = jsap.parse(args);

        if (result.getBoolean("help")) {
            System.out.println(jsap.getHelp());
            System.exit(0);
        }

        if (!result.success()) {
            final Iterator<String> errors = result.getErrorMessageIterator();
            while (errors.hasNext()) {
                System.err.println(errors.next());
            }
            System.err.println(jsap.getUsage());
            System.exit(1);
        }

        final boolean verbose = result.getBoolean("verbose");

        final File sequenceFile = result.getFile("input");
        if (verbose) {
            System.out.println("Reading sequence from: " + sequenceFile);
        }


        // extract the title to use for the output header
        final String title;
        if (result.contains("title")) {
            title = result.getString("title");
        } else {
            title = sequenceFile.getName();
        }

        Reader inputReader = null;
        PrintWriter outputWriter = null;

        try {
            if ("gz".equals(FilenameUtils.getExtension(sequenceFile.getName()))) {
                inputReader = new InputStreamReader(
                        new GZIPInputStream(FileUtils.openInputStream(sequenceFile)));
            } else {
                inputReader = new FileReader(sequenceFile);
            }
            final FastaParser fastaParser = new FastaParser(inputReader);

            final File outputFile = result.getFile("output");
            final OutputStream outputStream;
            if (outputFile != null) {
                outputStream = FileUtils.openOutputStream(outputFile);
                if (verbose) {
                    System.out.println("Writing sequence : " + outputFile);
                }
            } else {
                outputStream = System.out;
            }
            outputWriter = new PrintWriter(outputStream);

            // write the header portion of the output
            outputWriter.print("# ");
            outputWriter.print(new Date());
            outputWriter.print(' ');
            outputWriter.print(ColorSpaceConverter.class.getName());
            for (final String arg : args) {
                outputWriter.print(' ');
                outputWriter.print(arg);
            }
            outputWriter.println();
            outputWriter.print("# Cwd: ");
            outputWriter.println(new File(".").getCanonicalPath());
            outputWriter.print("# Title: ");
            outputWriter.println(title);

            // now parse the input sequence
            long sequenceCount = 0;
            final MutableString descriptionLine = new MutableString();
            final MutableString sequence = new MutableString();
            final MutableString colorSpaceSequence = new MutableString();

            while (fastaParser.hasNext()) {
                fastaParser.next(descriptionLine, sequence);
                outputWriter.print('>');
                outputWriter.println(descriptionLine);


                convert(sequence, colorSpaceSequence);
                CompactToFastaMode.writeSequence(outputWriter, colorSpaceSequence);

                sequenceCount++;
                if (verbose && sequenceCount % 10000 == 0) {
                    System.out.println("Converted " + sequenceCount + " entries");
                }
            }

            if (verbose) {
                System.out.println("Conversion complete!");
            }
        } finally {
            IOUtils.closeQuietly(inputReader);
            IOUtils.closeQuietly(outputWriter);
        }
    }
}
