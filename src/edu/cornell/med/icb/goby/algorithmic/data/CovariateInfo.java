package edu.cornell.med.icb.goby.algorithmic.data;

import edu.cornell.med.icb.io.TSVReader;
import edu.cornell.med.icb.io.TsvToFromMap;
import edu.cornell.med.icb.iterators.TsvLineIterator;
import edu.cornell.med.icb.maps.LinkedHashToMultiTypeMap;
import it.unimi.dsi.fastutil.objects.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

/**
 * Store and query information about sample covariates. This class can load covariate information from files.
 * These files must be formatted as TSV, with one column called sample-id,
 * and other columns named after keys. Elements of the columns are the values associated to the key for the sample-id.
 * Below is an example of a valid covariate file with four samples:
 * <p/>
 * <table border="1" style="border-collapse:collapse">
 * <tr><td>sample-id</td><td>Gender</td><td>Type</td><td>Kind-Of-Sample</td><td>Tissue</td></tr>
 * <tr><td>S1</td><td>Male</td><td>Father</td><td>Germline</td><td>Blood</td></tr>
 * <tr><td>S2</td><td>Female</td><td>Mother</td><td>Germline</td><td>Blood</td></tr>
 * <tr><td>S3</td><td>Male</td><td>Patient</td><td>Somatic</td><td>Blood</td></tr>
 * <tr><td>S4</td><td>Male</td><td>Patient</td><td>Germline</td><td>Skin</td></tr></table>
 * In this example, sample S3 will be associated with Kind-Of-Sample=Somatic, while samples S1,S2 and S4 will be
 * associated with Kind-Of-Sample=Germline.
 * <p/>
 * The covariate info data structure makes it possible to query for samples with Kind-Of-Sample=Somatic (see method
 * samplesWithExactCovariate(key,value).
 *
 * @author Fabien Campagne
 *         Date: 3/9/13
 *         Time: 10:26 AM
 */
public class CovariateInfo {

    Object2ObjectMap<String, Object2ObjectMap<String, String>> map = new Object2ObjectAVLTreeMap<String, Object2ObjectMap<String, String>>();


    public void defineSample(String sampleId) {
        if (!map.containsKey(sampleId)) {
            map.put(sampleId, new Object2ObjectAVLTreeMap<String, String>());
        }
    }

    /**
     * Retrieve the set of samples that satisfy the covariate relation key=value.
     *
     * @param key   Covariate/Key that sample must have to be returned.
     * @param value CovariateValue that sample must have to be returned.
     */
    public ObjectArraySet<String> samplesWithExactCovariate(String key, String value) {
        ObjectArraySet<String> result = new ObjectArraySet<String>();
        for (String sampleId : map.keySet()) {
            Object2ObjectMap<String, String> objectMap = map.get(sampleId);
            if (objectMap.containsKey(key)) {
                if (objectMap.get(key).equals(value)) {
                    result.add(sampleId);
                }
            }
        }
        return result;
    }

    /**
     * Retrieve the set of samples that satisfy the covariate relation key.contains(value). In contrast to
     * samplesWithExactCovariate, this method looks for samples whose covariate string contains the value.
     * This means that you can encode covariates with multiple values separated by some character and lookup
     * all samples that include the covariate element. For instance, in the following table, it is possible
     * to lookup samples that have S1 as a parent by calling samplesContainCovariate("Parents",S1): {S3,S4}:
     * Note that lookups are case sensitive.
     * <p/>
     * <table border="1" style="border-collapse:collapse">
     * <tr><td>sample-id</td><td>Gender</td><td>Type</td><td>Kind-Of-Sample</td><td>Tissue</td><td>Parents</td><td>Offspring</td></tr>
     * <tr><td>S1</td><td>Male</td><td>Father</td><td>Germline</td><td>Blood</td><td>N/A</td><td>S3|S4</td></tr>
     * <tr><td>S2</td><td>Female</td><td>Mother</td><td>Germline</td><td>Blood</td><td>N/A</td><td>S3|S4</td></tr>
     * <tr><td>S3</td><td>Male</td><td>Patient</td><td>Somatic</td><td>Blood</td><td>S1|S2</td><td>N/A</td></tr>
     * <tr><td>S4</td><td>Male</td><td>Patient</td><td>Germline</td><td>Skin</td><td>S1|S2</td><td>N/A</td></tr>
     * </table>
     *
     * @param key   Covariate/Key that sample must have to be returned.
     * @param value CovariateValue that sample must include to be returned.
     */
    public ObjectArraySet<String> samplesContainCovariate(String key, String value) {
        ObjectArraySet<String> result = new ObjectArraySet<String>();
        for (String sampleId : map.keySet()) {
            Object2ObjectMap<String, String> objectMap = map.get(sampleId);
            if (objectMap.containsKey(key)) {
                if (objectMap.get(key).contains(value)) {
                    result.add(sampleId);
                }
            }
        }
        return result;
    }

    /**
     * Indicate that the random access cache was created.
     */
    private boolean randomAccessCacheCreated;

    private Object2ObjectAVLTreeMap<String, ObjectArraySet<String>> covariateToSampleSetCache =
            new Object2ObjectAVLTreeMap<String, ObjectArraySet<String>>();
    ObjectSet<String> allCovariates;
    private void refresh() {

        if (!randomAccessCacheCreated) {
            allCovariates = new ObjectArraySet<String>();
            for (Map.Entry<String, Object2ObjectMap<String, String>> entry : map.entrySet()) {
                ObjectSet<String> covariates = entry.getValue().keySet();
                allCovariates.addAll(covariates);
            }

        }
        randomAccessCacheCreated = true;
    }


    /**
     * Parse a covariate information file.
     *
     * @param covInfoFilename
     * @return
     * @throws IOException
     */
    public static CovariateInfo parse(String covInfoFilename) throws IOException {
        File dataFile = new File(covInfoFilename);
        CovariateInfo result = new CovariateInfo();

        final TsvToFromMap tsvReader = TsvToFromMap.createFromTsvFile(dataFile);
        for (final LinkedHashToMultiTypeMap<String> dataline : new TsvLineIterator(dataFile, tsvReader)) {
            for (String key : tsvReader.getColumnHeaders()) {
                if (!"sample-id".equals(key)) {
                    String sampleId = dataline.getString("sample-id");
                    String value = dataline.getString(key);
                    result.defineSample(sampleId);
                    result.associate(sampleId, key, value);
                }
            }
        }
        return result;

    }

    private void associate(String sampleId, String key, String value) {
        Object2ObjectMap<String, String> keyValuePairs = map.get(sampleId);
        if (keyValuePairs == null) {
            throw new IllegalArgumentException("SampleId is not defined in the covariates map, cannot associate key value pair to it.");

        }
        keyValuePairs.put(key, value);
    }

    /**
     * Return the value associated to a sample for a given covariate.
     *
     * @param sampleId     sample of interest.
     * @param covariateKey Key of the covariate for which the value is sought.
     * @return null when either the sampleId or covariate key is not defined, or the value of the associated covariate.
     */
    public String getCovariateValue(String sampleId, String covariateKey) {

        Object2ObjectMap<String, String> sampleMap = map.get(sampleId);
        if (sampleMap == null) {
            return null;
        }
        return sampleMap.get(covariateKey);
    }

    /**
     * Returns true iff the sample has the specific value of the covariate.
     *
     * @param sampleId
     * @param covariateKey
     * @param value
     * @return
     */
    public boolean hasCovariateValue(String sampleId, String covariateKey, String value) {
        String v = getCovariateValue(sampleId, covariateKey);
        if (v == null) return false;
        return (value.equals(v));
    }

    public ObjectSet<String> getCovariateKeys() {
        refresh();
        return allCovariates;
    }
}
