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

import com.martiansoftware.jsap.Flagged;
import com.martiansoftware.jsap.IDMap;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Option;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.Switch;
import org.apache.commons.lang.BooleanUtils;
import org.apache.commons.lang.ClassUtils;
import org.apache.commons.lang.StringEscapeUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.WordUtils;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Map;

/**
 * Base abstract class for Goby modes.
 *
 * @author Kevin Dorff
 */
public abstract class AbstractCommandLineMode {
    private String jarFilename;

    protected AbstractCommandLineMode(final String jarFilename) {
        setJarFilename(jarFilename);
    }

    public String getJarFilename() {
        return jarFilename;
    }

    public void setJarFilename(final String jarFilename) {
        this.jarFilename = jarFilename;
    }

    /**
     * Returns the mode name defined by subclasses.
     *
     * @return The name of the mode
     */
    public abstract String getModeName();

    /**
     * Returns the mode description defined by subclasses.
     *
     * @return A description of the mode
     */
    public abstract String getModeDescription();

    /**
     * Returns true when the id represents a help option.
     *
     * @param id The id to check
     * @return true if the id is a known help option.
     */
    private boolean isJSAPHelpId(final String id) {
        return "help".equals(id) || "wikihelp".equals(id) || "htmlhelp".equals(id);
    }

    /**
     * Print usage for the mode.
     *
     * @param jsap the configured JSAP object
     */
    public void printUsage(final JSAP jsap) {
        final String modeName = getModeName();
        final String modeDescription = getModeDescription();
        final Parameter modeId = jsap.getByID("mode");
        if (modeId != null) {
            jsap.unregisterParameter(modeId);
        }
        System.err.println("Usage: ");
        System.err.println("java -jar " + jarFilename
                + " (-m|--mode) " + modeName + " " + jsap.getUsage());
        System.err.println();
        System.err.println("Description: ");
        System.err.println(WordUtils.wrap(modeDescription, JSAP.DEFAULT_SCREENWIDTH - 1));
        System.err.println();
        System.err.println("Options: ");
        System.err.println(jsap.getHelp(JSAP.DEFAULT_SCREENWIDTH - 1));
    }

    /**
     * Print the usage for a mode in the format appropriate for a
     * <a href="http://www.wikimedia.org/">Wikimedia</a> table.
     *
     * @param jsap The parameters for the mode
     */
    public void printUsageAsWikiTable(final JSAP jsap) {
        final PrintStream stream = System.out;

        // Table definition/header
        stream.println("{| class=\"wikitable\" style=\"margin: 1em auto 1em auto; background:#efefef;\" valign=\"top\" align=\"center\" border=\"1\" cellpadding=\"5\" width=\"80%\"");
        stream.println("|- valign=\"bottom\" align=\"left\" style=\"background:#ffdead;\"");
        stream.println("!width=\"15%\"| Flag");
        stream.println("!width=\"10%\" | Arguments");
        stream.println("!width=\"5%\"| Required");
        stream.println("! Description");
        stream.println();

        // iterate through help entries for the mode
        final IDMap idMap = jsap.getIDMap();
        final Iterator iterator = idMap.idIterator();
        while (iterator.hasNext()) {
            final String id = (String) iterator.next();
            if (!"mode".equals(id) && !isJSAPHelpId(id)) { // skip the mode and help ids
                final Parameter parameter = jsap.getByID(id);
                stream.println("|- valign=\"top\" align=\"left\"");
                stream.print("| <nowiki>");
                if (parameter instanceof Flagged) {
                    final Flagged flagged = (Flagged) parameter;
                    final Character characterFlag = flagged.getShortFlagCharacter();
                    final String longFlag = flagged.getLongFlag();
                    if (characterFlag != null && StringUtils.isNotBlank(longFlag)) {
                        stream.print("(-" + characterFlag + "|--" + longFlag + ")");
                    } else if (characterFlag != null) {
                        stream.print("-" + characterFlag);
                    } else if (StringUtils.isNotBlank(longFlag)) {
                        stream.print("--" + longFlag);
                    }
                } else {
                    stream.print("n/a");
                }
                stream.println("</nowiki>");
                stream.print("| ");
                if (parameter instanceof Switch) {
                    stream.println("n/a");
                } else {
                    stream.println(id);
                }

                final boolean required;
                if (parameter instanceof Option) {
                    final Option option = (Option) parameter;
                    required = option.required();
                } else {
                    required = !(parameter instanceof Switch);
                }
                stream.println("| " + BooleanUtils.toStringYesNo(required));
                stream.println("| " + parameter.getHelp());
                if (parameter.getDefault() != null) {
                    stream.print("Default value: ");
                    for (final String defaultValue : parameter.getDefault()) {
                        stream.print(" ");
                        stream.print(defaultValue);
                    }
                }
                stream.println();
            }
        }

        // table close
        stream.println("|}");
    }

    /**
     * Print the usage for a mode in the format appropriate for a table in HTML.
     *
     * @param jsap The parameters for the mode
     */
    public void printUsageAsHtmlTable(final JSAP jsap) {
        final PrintStream writer = System.out;

        // Table definition/header
        writer.println("<table style=\"background: #efefef;\" border=\"1\" cellpadding=\"5\">");
        writer.println("<tbody>");
        writer.println("<tr style=\"background: #ffdead;\" align=\"left\" valign=\"bottom\">");
        writer.println("<th style=\"width: 15%\">Flag</th>");
        writer.println("<th style=\"width: 10%\">Arguments</th>");
        writer.println("<th style=\"width: 5%\">Required</th>");
        writer.println("<th>Description</th>");
        writer.println("</tr>");

        // iterate through help entries for the mode
        final IDMap idMap = jsap.getIDMap();
        final Iterator iterator = idMap.idIterator();
        while (iterator.hasNext()) {
            final String id = (String) iterator.next();
            if (!"mode".equals(id) && !isJSAPHelpId(id)) { // skip the mode and help ids
                final Parameter parameter = jsap.getByID(id);
                writer.println("<tr align=\"left\" valign=\"top\">");
                writer.print("<td><code>");
                if (parameter instanceof Flagged) {
                    final Flagged flagged = (Flagged) parameter;
                    final Character characterFlag = flagged.getShortFlagCharacter();
                    final String longFlag = flagged.getLongFlag();
                    if (characterFlag != null && StringUtils.isNotBlank(longFlag)) {
                        writer.print("(-" + characterFlag + "|--" + longFlag + ")");
                    } else if (characterFlag != null) {
                        writer.print("-" + characterFlag);
                    } else if (StringUtils.isNotBlank(longFlag)) {
                        writer.print("--" + longFlag);
                    }
                } else {
                    writer.print("n/a");
                }
                writer.println("</code></td>");
                writer.print("<td>");
                if (parameter instanceof Switch) {
                    writer.print("n/a");
                } else {
                    writer.print(id);
                }
                writer.println("</td>");
                final boolean required;
                if (parameter instanceof Option) {
                    final Option option = (Option) parameter;
                    required = option.required();
                } else {
                    required = !(parameter instanceof Switch);
                }
                writer.println("<td>" + BooleanUtils.toStringYesNo(required) + "</td>");
                writer.print("<td>");

                // convert any strange characters to html codes
                final String htmlHelpString = StringEscapeUtils.escapeHtml(parameter.getHelp());
                // and also convert "-" to the hyphen character code
                writer.print(StringUtils.replace(htmlHelpString, "-", "&#45;"));
                 if (parameter.getDefault() != null) {
                    writer.print(" Default value: ");
                    for (final String defaultValue : parameter.getDefault()) {
                        writer.print(" ");
                        writer.print(StringUtils.replace(defaultValue, "-", "&#45;"));
                    }
                }
                writer.println("</td>");
                writer.println("</tr>");
            }
        }

        // table close
        writer.println("</tbody>");
        writer.println("</table>");
    }

    /**
     * Load JSAP from resource with no modifications to help / defaults.
     *
     * @return the configured JSAP object
     * @throws IOException   error reading or configuring
     * @throws JSAPException error reading or configuring
     */
    public JSAP loadJsapFromResource() throws IOException, JSAPException {
        return loadJsapFromResource(null);
    }

    /**
     * Load JSAP XML configuration for the current class. With potential
     * modifications to help / defaults stored in the helpValues map.
     *
     * @param helpValues a map of values to replace in the jsap help with other values
     * @return the configured JSAP object
     * @throws IOException   error reading or configuring
     * @throws JSAPException error reading or configuring
     */
    @SuppressWarnings("unchecked")
    public JSAP loadJsapFromResource(final Map<String, String> helpValues) throws IOException, JSAPException {
        final Class thisClass = this.getClass();
        final String className = ClassUtils.getShortCanonicalName(thisClass);
        JSAP jsap;
        try {
            jsap = new JSAP(thisClass.getResource(className + ".jsap"));
        } catch (NullPointerException e) {  // NOPMD
            // No JSAP for "this", no additional args for the specific driver
            jsap = new JSAP();
        }

        // add options from driver to mode specific JSAP:
        final JSAP jsapSource = new JSAP(GenericToolsDriver.class.getResource(
                ClassUtils.getShortCanonicalName(GenericToolsDriver.class) + ".jsap"));
        final Iterator<String> idsIt = jsapSource.getIDMap().idIterator();
        while (idsIt.hasNext()) {
            final String id = idsIt.next();
            // remove parameters if they already exist in the destination.
            // This makes sure we use the source parameters, such as mode, help,
            // which contain variables to be replaced.
            if (jsap.getByID(id) != null) {
                jsap.unregisterParameter(jsap.getByID(id));
            }
            jsap.registerParameter(jsapSource.getByID(id));

        }
        if (helpValues == null) {
            return jsap;
        }
        final IDMap idMap = jsap.getIDMap();
        final Iterator<String> idIterator = idMap.idIterator();
        while (idIterator.hasNext()) {
            final String id = idIterator.next();
            final Parameter param = jsap.getByID(id);
            final String help = param.getHelp();
            final String[] defaults = param.getDefault();
            for (final Map.Entry<String, String> entry : helpValues.entrySet()) {
                // replace values in help
                if (help.contains(entry.getKey())) {
                    param.setHelp(StringUtils.replace(help, entry.getKey(), entry.getValue()));
                }
                // replace values in defaults
                if (defaults != null) {
                    for (int i = 0; i < defaults.length; i++) {
                        if (defaults[i].contains(entry.getKey())) {
                            defaults[i] = StringUtils.replace(
                                    defaults[i], entry.getKey(), entry.getValue());
                        }
                    }
                }
            }
        }

        return jsap;
    }

    /**
     * Parse the JSAP to JSAPResult. Handle if "help" was requested.
     *
     * @param jsap the parser to use when parsing arguments
     * @param args the arguments to parse
     * @return the JSAPResult
     */
    public JSAPResult parseJsap(final JSAP jsap, final String[] args) {
        final JSAPResult jsapResult = jsap.parse(args);
        if (getModeName() != null) {
            if (jsap.getByID("help") != null && jsapResult.getBoolean("help")) {
                printUsage(jsap);
                System.exit(1);
            } else if (jsap.getByID("htmlhelp") != null && jsapResult.getBoolean("htmlhelp")) {
                printUsageAsHtmlTable(jsap);
                System.exit(2);
            } else if (jsap.getByID("wikihelp") != null && jsapResult.getBoolean("wikihelp")) {
                printUsageAsWikiTable(jsap);
                System.exit(3);
            }
        }
        return jsapResult;
    }

    /**
     * Configure the mode via command line arguments.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException   error configuring
     * @throws JSAPException error configuring
     */
    public abstract AbstractCommandLineMode configure(final String[] args)
            throws IOException, JSAPException;

    /**
     * Execute the mode.
     *
     * @throws IOException io error
     */
    public abstract void execute() throws IOException;

    /**
     * Parse the JSAP arguments defined for the mode. Different arguments can be defined for
     * each mode.
     *
     * @param args       command line arguments.
     * @param helpValues A map of option names to help strings
     * @return Parsed arguments.
     * @throws IOException   if there was a problem configuring the parser
     * @throws JSAPException if there was a problem parsing the args
     * @see #loadJsapFromResource()
     * @see #loadJsapFromResource(java.util.Map)
     */
    protected JSAPResult parseJsapArguments(final String[] args, final Map<String, String> helpValues)
            throws IOException, JSAPException {
        final JSAP jsap = loadJsapFromResource(helpValues);
        final JSAPResult jsapResult = parseJsap(jsap, args);
        abortOnError(jsap, jsapResult);
        return jsapResult;
    }

    /**
     * Parse the JSAP arguments defined for the mode. Different arguments can be defined for
     * each mode.
     *
     * @param args command line arguments.
     * @return Parsed arguments.
     * @throws IOException   if there was a problem parsing the help text
     * @throws JSAPException if there was a problem parsing the args
     * @see #loadJsapFromResource()
     * @see #loadJsapFromResource(java.util.Map)
     */
    protected JSAPResult parseJsapArguments(final String[] args) throws IOException, JSAPException {
        final JSAP jsap = loadJsapFromResource();
        final JSAPResult jsapResult = parseJsap(jsap, args);
        abortOnError(jsap, jsapResult);
        return jsapResult;
    }

    /**
     * Print out specific error messages describing the problems
     * with the command line, THEN print usage, THEN print full
     * help.  This is called "beating the user with a clue stick."
     *
     * @param jsap       The argument parser
     * @param jsapResult The results of the parse
     */
    private void abortOnError(final JSAP jsap, final JSAPResult jsapResult) {
        if (!jsapResult.success()) {
            System.err.println();
            for (final Iterator errs = jsapResult.getErrorMessageIterator(); errs.hasNext();) {
                System.err.println("Error: " + errs.next());
            }
            System.err.println();
            printUsage(jsap);
            System.exit(1);
        }
    }
}
