package edu.cornell.med.icb.goby.alignments;

import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.lang.MutableString;
import org.apache.commons.io.FilenameUtils;

import java.util.HashMap;
import java.util.Map;

/**
 * Class to help parse read group information.
 *
 * @author Fabien Campagne
 *         Date: 4/3/13
 *         Time: 4:43 PM
 */
public class ReadGroupHelper {

    public boolean isOverrideReadGroups() {
        return overrideReadGroups;
    }

    private IndexedIdentifier basenames = new IndexedIdentifier();

    public void setOverrideReadGroups(boolean overrideReadGroups) {
        this.overrideReadGroups = overrideReadGroups;
    }

    private boolean overrideReadGroups = false;
    /**
     * A map with read group associations: PL->platform.
     */
    private Map<String, String> readGroupMap = new HashMap<String, String>();


    /**
     * Parse the read group info associated with a set of samples from the command line options. The JSAP
     * argument include-reference-names must be defined.
     *
     * @param jsapResult     The jsapResult available to the mode.
     * @param inputFilenames Input files for alignments associated with this slice.
     */
    public void parseReadGroupOptions(final JSAPResult jsapResult, String[] inputFilenames) {
        if (jsapResult.userSpecified("platform")) {
            readGroupMap.put("PL", jsapResult.getString("platform"));
        }
        for (String inputFilename : inputFilenames) {
            basenames.registerIdentifier(new MutableString(inputFilename));
        }

    }

    public Map<String, String> getReadGroupMap() {
        return readGroupMap;
    }

    /**
     * Not that this method currently ignores the input filename argument.
     *
     * @param inputFilename
     * @return
     */
    public String getPlatform(String inputFilename) {
        if (readGroupMap.containsKey("PL")) {
            return readGroupMap.get("PL");
        } else {
            return "platform";
        }
    }

    public String getId(String inputFilename) {
        basenames.registerIdentifier(new MutableString(inputFilename));
        return Integer.toString(basenames.getInt(new MutableString(inputFilename)));
    }

    public String getSample(String basename) {
        // remove the path from basename:
        return FilenameUtils.getBaseName(basename);
    }
}
