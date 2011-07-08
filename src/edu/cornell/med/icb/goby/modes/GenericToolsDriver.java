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

import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import edu.cornell.med.icb.io.ResourceFinder;
import org.apache.commons.lang.ClassUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.lang.reflect.Modifier;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * When running java -jar goby.jar this is the class that is called.
 * Using the --mode argument, this calls one of several classes and passes the
 * command line arguments to those classes.
 *
 * @author Kevin Dorff
 */
public class GenericToolsDriver extends AbstractCommandLineMode {
    private static final Log LOG = LogFactory.getLog(GenericToolsDriver.class);

    /**
     * The mode to run.
     */
    private String mode;

    /**
     * The JSAP configuration.
     */
    private JSAP jsap;

    /**
     * The command line arguments.
     */
    private String[] args;

    /**
     * Mode name to Class map.
     */
    private final Map<String, Class> MODES_MAP;

    /**
     * Map to override help / default values.
     */
    protected final Map<String, String> HELP_VALUES;

    @Override
    public String getModeName() {
        return null;
    }

    @Override
    public String getModeDescription() {
        return null;
    }

    public GenericToolsDriver(final String jarFilename) {
        super(jarFilename);
        MODES_MAP = loadModeMap();

        HELP_VALUES = new HashMap<String, String>();
        final StringBuilder modesList = new StringBuilder();
        final List<String> modes = new ArrayList<String>(MODES_MAP.keySet());
        Collections.sort(modes);
        for (final String mode : modes) {
            modesList.append("* ").append(mode).append('\n');
        }
        HELP_VALUES.put("[MODES_LIST]", modesList.toString());
    }

    /**
     * Display JSAP usage information (overrides the base version).
     *
     * @param jsapVal the JSAP as configured
     */
    @Override
    public void printUsage(final JSAP jsapVal) {
        try {
            System.err.println("java -jar " + getJarFilename() + " " + jsapVal.getUsage());
            System.err.println();
            System.err.println(jsapVal.getHelp());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Create the GenericToolsDriver object, parse argument(s). This intentionaly doesn't check
     * for errors as we only care for "mode" and we KNOW there will be other command line
     * parameters. All other arguments are specifed by the specific mode being called
     * (including things like --help).
     *
     * @param argsVal the arguments
     * @return this object for chaining
     * @throws IOException   error parsing arguments
     * @throws JSAPException error parsing arguments
     */
    @Override
    public GenericToolsDriver configure(final String[] argsVal) throws IOException, JSAPException {
        jsap = loadJsapFromResource(HELP_VALUES);
        args = argsVal;
        final JSAPResult jsapResult = parseJsap(jsap, args);
        mode = jsapResult.getString("mode");
        if (mode == null) {
            unregisterFormattedHelp(jsap);
            printUsage(jsap);
            System.exit(1);
        }
        return this;
    }

    /**
     * Removes any "special" options from being displayed.  In this case,
     * these are help strings that are only useful with modes.
     *
     * @param jsapVal the JSAP as configured
     */
    private void unregisterFormattedHelp(final JSAP jsapVal) {
        // unregister help parameters for modes
        final Parameter htmlHelpId = jsapVal.getByID("htmlhelp");
        if (htmlHelpId != null) {
            jsapVal.unregisterParameter(htmlHelpId);
        }
        final Parameter wikiHelpId = jsapVal.getByID("wikihelp");
        if (wikiHelpId != null) {
            jsapVal.unregisterParameter(wikiHelpId);
        }
    }

    /**
     * Execute the appropriate mode.
     *
     * @throws IOException error parsing arguments (in the mode), or other exception
     */
    @Override
    public void execute() throws IOException {
        final Class modeClass = MODES_MAP.get(mode);
        Exception caught = null;
        try {
            if (modeClass != null) {
                final AbstractCommandLineMode modeObject = (AbstractCommandLineMode) modeClass.newInstance();
                modeObject.configure(args);
                modeObject.execute();
            }
        } catch (JSAPException e) {
            caught = e;
        } catch (IllegalAccessException e) {
            caught = e;
        } catch (InstantiationException e) {
            caught = e;
        } catch (OutOfMemoryError e) {
            e.printStackTrace();
            System.err.println(
                    "\n!!\n"
                            + "!! An out of memory exception was thrown."
                            + "!! Try running with more memory such as -Xmx1200m -XX:PermSize=384m\n"
                            + "!!\n");
            throw e;
        }

        if (modeClass == null || caught != null) {
            if (modeClass == null) {
                System.err.println("Unrecognized mode: '" + mode + "'");
                unregisterFormattedHelp(jsap);
            }
            printUsage(jsap);
            if (caught != null) {
                System.err.println();
                System.err.println("Exception caught type=" + caught.getClass().toString());
                System.err.println("   Message=" + caught.getMessage());
                caught.printStackTrace(System.err);

            }
            System.exit(1);
        }
        System.exit(0);
    }

    /**
     * Load the list of concrete modes.
     *
     * @return the mode map
     */
    private Map<String, Class> loadModeMap() {
        final Map<String, Class> modeMap = new HashMap<String, Class>();
        for (final String modeClassName : modesClassNamesList()) {
            try {
                LOG.debug("About to process modeClassName: " + modeClassName);
                final Class modeClass = Class.forName(modeClassName);
                if (Modifier.isAbstract(modeClass.getModifiers())) {
                    // Ignore abstract classes
                    continue;
                }
                // Since we're not abstract, we can make an instance
                final Object modeInstance = modeClass.newInstance();
                if (modeInstance instanceof AbstractCommandLineMode) {
                    final String modeName = ((AbstractCommandLineMode) modeInstance).getModeName();
                    // Make sure the class has a static field "MODE_NAME"
                    if (modeName != null) {
                        modeMap.put(modeName, modeClass);
                    }
                }
            } catch (ClassNotFoundException e) {
                System.err.println(
                        "Could find a class for " + modeClassName
                                + " ClassNotFoundException: " + e.getMessage());
            } catch (IllegalAccessException e) {
                System.err.println(
                        "Could not find MODE_NAME for class " + modeClassName
                                + " IllegalAccessException: " + e.getMessage());
            } catch (InstantiationException e) {
                System.err.println(
                        "Could not find MODE_NAME for class " + modeClassName
                                + " InstantiationException: " + e.getMessage());
            }
        }
        return modeMap;
    }

    /**
     * Retrieve the class names of the likely modes.
     *
     * @return the list of class names
     */
    private List<String> modesClassNamesList() {
        final List<String> modesList = new LinkedList<String>();
        try {
            final String modesPackage = ClassUtils.getPackageName(getClass());
            final String[] files = new ResourceFinder().getResourceListing(getClass());
            for (final String file : files) {
                if (file.endsWith("Mode.class")) {
                    // Keep only classes whose class name end in "Mode"
                    final String modeClassName = modesPackage + "."
                            + file.substring(0, file.length() - 6);
                    modesList.add(modeClassName);
                }
            }
        } catch (URISyntaxException e) {
            System.err.println("Could not list mode class names URISyntaxException: "
                    + e.getMessage());
            e.printStackTrace();
        } catch (IOException e) {
            System.err.println("Could not list mode class names IOException: "
                    + e.getMessage());
            e.printStackTrace();
        }
        return modesList;
    }
}
