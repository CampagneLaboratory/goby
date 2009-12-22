/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.goby.modes;

import com.martiansoftware.jsap.IDMap;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import org.apache.commons.lang.ClassUtils;
import org.apache.commons.lang.StringUtils;

import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

/**
 * Base abstract class for MaqToolsMode's.
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
     * Returns the mode name defined by subclasses, OR null
     */
    public abstract String getModeName();

    /**
     * Returns the mode description defined by subclasses, OR null
     */
    public abstract String getModeDescription();

    /**
     * Print Usage for all of the Cascading Mode (except CascadingDriver).
     *
     * @param jsap     the configured JSAP object
     */
    public void printUsage(final JSAP jsap) {
        final String modeName = getModeName();
        final String modeDescription = getModeDescription();
        final Parameter modeId = jsap.getByID("mode");
        if (modeId != null) {
            jsap.unregisterParameter(modeId);
        }
        System.err.println("Usage: ");
        System.err.println("java -jar " + jarFilename + " --mode "
                + modeName + " " + jsap.getUsage());
        System.err.println();
        System.err.println("Description: ");
        System.err.println(modeDescription);
        System.err.println();
        System.err.println("Options: ");
        System.err.println(jsap.getHelp());
    }

    /**
     * Load JSAP from resource with no modifications to help / defaults.
     *
     * @return the configured JSAP object
     * @throws IOException error reading or configuring
     * @throws JSAPException
     *                             error reading or configuring
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
     * @throws IOException error reading or configuring
     * @throws JSAPException
     *                             error reading or configuring
     */
    @SuppressWarnings("unchecked")
    public JSAP loadJsapFromResource(final Map<String, String> helpValues)
            throws IOException, JSAPException {
        final Class thisClass = this.getClass();
        final String className = ClassUtils.getShortCanonicalName(thisClass);
        JSAP jsap;
        try {
            jsap = new JSAP(thisClass.getResource(className + ".jsap"));
        } catch (NullPointerException e) {
            // No JSAP for "this", no additional args for the specific driver
            jsap = new JSAP();
        }

        // add options from driver to mode specific JSAP:
        final JSAP jsapSource = new JSAP(GenericToolsDriver.class.getResource(
                ClassUtils.getShortCanonicalName(GenericToolsDriver.class) + ".jsap"));
        final Iterator<String> idsIt = jsapSource.getIDMap().idIterator();
        while (idsIt.hasNext()) {
            final String id = idsIt.next();
            //remove parameters if they already exist in the destination.
            // This makes sure we use the source parameters, such as mode, help, which contain variables to be replaced.
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

    private JSAP findJSAP(final Class thisClass) throws IOException, JSAPException {
        final String className = ClassUtils.getShortCanonicalName(thisClass);
        final JSAP jsap;
        try {
            jsap = new JSAP(thisClass.getResource(className + ".jsap"));
            return jsap;
        } catch (NullPointerException e) {
            // No JSAP for "this", no additional args for the specific driver
            return null;
        }

    }

    /**
     * Parse the JSAP to JSAPResult. Handle if "--help" was requested.
     *
     * @param jsap     the JSAP to use to parse args
     * @param args     the args
     * @return the JSAPResult
     */
    public JSAPResult parseJsap(final JSAP jsap, final String[] args) {
        final JSAPResult jsapResult = jsap.parse(args);
        if (getModeName() != null) {
            if (jsap.getByID("help") != null && jsapResult.getBoolean("help")) {
                printUsage(jsap);
                System.exit(1);
            }
        }
        return jsapResult;
    }

    /**
     * Configure the mode via command line arguments.
     *
     * @param args command line arguments
     * @return this object for chaining
     * @throws IOException error configuring
     * @throws JSAPException
     *                             error configuring
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
     * each mode, see   loadJsapFromResource
     *
     * @param args command line arguments.
     * @param helpValues
     * @return Parsed arguments.
     * @throws IOException
     * @throws JSAPException
     */
    protected JSAPResult parseJsapArguments(final String[] args, final Map<String, String> helpValues) throws IOException, JSAPException {
        final JSAP jsap = loadJsapFromResource(helpValues);
        final JSAPResult jsapResult = parseJsap(jsap, args);
        abortOnError(jsap, jsapResult);
        return jsapResult;
    }

    /**
     * Parse the JSAP arguments defined for the mode. Different arguments can be defined for
     * each mode, see   loadJsapFromResource
     *
     * @param args command line arguments.
     * @return Parsed arguments.
     * @throws IOException
     * @throws JSAPException
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
     */
    private void abortOnError(final JSAP jsap, final JSAPResult jsapResult) {
        if (!jsapResult.success()) {
            System.err.println();
            for (Iterator errs = jsapResult.getErrorMessageIterator();
                 errs.hasNext();) {
                System.err.println("Error: " + errs.next());
            }
            System.err.println();
            printUsage(jsap);
            System.exit(1);
        }
    }
}
