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

package edu.cornell.med.icb.goby.stats;

import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.algorithmic.data.GroupComparison;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import it.unimi.dsi.fastutil.objects.*;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.ArrayList;
import java.util.ServiceLoader;

/**
 * Created by IntelliJ IDEA.
 * User: nyasha
 * Date: Mar 24, 2010
 * Time: 3:29:26 PM
 * Class that organizes group comparisons and facillitates differential expression
 * statistical testing.
 */
public class DifferentialExpressionAnalysis {
    /**
     * Used to log debug and informational messages.
     */
    private static final Log LOG = LogFactory.getLog(DifferentialExpressionAnalysis.class);

    /**
     * The groups that should be compared, order matters.
     */
    private String[] groupComparison;

    private String[] groups;
    private final Object2ObjectMap<String, Integer> groupSizes =
            new Object2ObjectOpenHashMap<String, Integer>();

    /**
     * Flag for an invalid t-test due to any group size < 2.
     */
    private boolean ttestflag;

    /**
     * The set of normalization methods to use for the comparison.
     */

    private ObjectArraySet<NormalizationMethod> normalizationMethods;


    private static final ServiceLoader<NormalizationMethod> normalizationMethodLoader
            = ServiceLoader.load(NormalizationMethod.class);
    private boolean runInParalell;

    public DifferentialExpressionAnalysis() {
        super();
    }

    public void parseGroupsDefinition(final String groupsDefinition,
                                      final DifferentialExpressionCalculator deCalculator,
                                      final String[] inputFilenames) {
        if (groupsDefinition == null) {
            // no groups definition to parse.
            return;
        }

        final String[] groupsTmp = groupsDefinition.split("/");
        ObjectSet<String> groupSet = new ObjectArraySet<String>();
        for (final String group : groupsTmp) {
            final String[] groupTokens = group.split("=");
            if (groupTokens.length < 2) {
                System.err.println("The --group argument must have the syntax groupId=basename");
                System.exit(1);
            }
            final String groupId = groupTokens[0];
            final String groupBasenames = groupTokens[1];
            assert groupTokens.length == 2 : "group definition must have only two elements separated by an equal sign.";
            deCalculator.defineGroup(groupId);
            groupSet.add(groupId);
            for (final String groupString : groupBasenames.split(",")) {

                final String sampleBasename = FilenameUtils.getBaseName(AlignmentReaderImpl.getBasename(groupString));
                if (!isInputBasename(sampleBasename, inputFilenames)) {
                    System.err.printf("The group basename %s is not a valid input basename.%n", sampleBasename);
                    System.exit(1);
                }
                System.out.println("Associating basename: " + sampleBasename + " to group: " + groupId);
                deCalculator.associateSampleToGroup(sampleBasename, groupId);
            }
            final int groupSize = (groupBasenames.split(",")).length;
            groupSizes.put(groupId, groupSize);
        }
        groups = groupSet.toArray(new String[groupSet.size()]);
    }

    public String[] getGroups() {
        return groups;
    }

    /**
     * Return true if the basename is on the command line as an input basename.
     *
     * @param basename
     * @return
     */
    private boolean isInputBasename(final String basename, final String[] inputFilenames) {
        for (final String inputFilename : inputFilenames) {
            if (FilenameUtils.getBaseName(AlignmentReaderImpl.getBasename(inputFilename)).equals(basename)) {
                return true;
            }
        }
        return false;
    }

    ArrayList<GroupComparison> groupComparisons;

    public ArrayList<GroupComparison> parseCompare(final String compare) {
        groupComparisons = new ArrayList<GroupComparison>();
        if (compare == null) {
            return groupComparisons;
        }
        String comparisonPairs[] = compare.split(",");
        ObjectSet<String> groupSet = new ObjectArraySet<String>();
        groupSet.addAll(ObjectArrayList.wrap(groups));

        int comparisonIndex = 0;
        for (String pair : comparisonPairs) {
            final String[] groupLanguageText = pair.split("/");
            if (groupLanguageText.length > 2) {
                System.err.println("Group comparison description" + pair + " used in --compare must name at most two groups. "
                        + "Please see the --groups option to define groups.");
                System.exit(1);
            }
            if (groupLanguageText.length == 1) {
                // comparisons with only one group are ignored.
                continue;
            }
            if (groupLanguageText.length < 1) {
                System.err.println("Group comparison description" + pair + " used in --compare must name at least one groups. "
                        + "Please see the --groups option to define groups.");
                System.exit(1);
            }
            if (groupLanguageText[0].equals(groupLanguageText[1])) {
                System.err.printf("Group names must be different when specifying two groups %s/%s. Use a single group %s or leave out the --compare option to disable group comparisons.%n",
                        groupLanguageText[0],
                        groupLanguageText[0], groupLanguageText[0]);
                System.exit(1);
            } else {
                System.err.printf("Ignoring group comparison %s. Use A/B to compare two groups.%n", groupLanguageText[0]);
            }
            int groupIndex = -1;
            for (final String groupId : groupLanguageText) {
                if (!groupSet.contains(groupId)) {
                    System.err.println("Group " + groupId + " used in --compare must be defined. "
                            + "Please see the --groups option to define groups.");
                    System.exit(1);
                }
            }
            groupComparison = groupLanguageText;
            int firstGroupIndex = 0;
            int secondGroupIndex = groupLanguageText.length - 1;
            groupComparisons.add(new GroupComparison(groupLanguageText[firstGroupIndex],
                    groupLanguageText[secondGroupIndex],
                    groupIndex(groupLanguageText[firstGroupIndex]),
                    groupIndex(groupLanguageText[secondGroupIndex]),
                    comparisonIndex++));
        }
        return groupComparisons;

    }

    private int groupIndex(final String groupName) {
        int index = 0;
        for (final String s : getGroups()) {
            if (groupName.equals(s)) {
                return index;
            }
            index++;
        }
        assert false : "group " + groupName + " must be found.";
        return -1;
    }

    public ObjectArraySet<NormalizationMethod> parseNormalization(final JSAPResult jsapResult) {
        final String normalizationMethodNames = jsapResult.getString("normalization-methods");
        final String[] methodIds = normalizationMethodNames.split(",");
        this.normalizationMethods = new ObjectArraySet<NormalizationMethod>();
        LOG.info("Looking up services");
        for (final String methodId : methodIds) {
            for (final NormalizationMethod aMethod : normalizationMethodLoader) {
                if (aMethod.getIdentifier().equals(methodId)) {
                    LOG.info("Adding " + aMethod.getIdentifier());
                    this.normalizationMethods.add(aMethod);
                }
            }
        }
        if (this.normalizationMethods.size() == 0) {
            LOG.error("Could not locate any normalization method with the names provided: " + normalizationMethodNames);
            System.exit(1);
        }
        return normalizationMethods;
    }

    public boolean checkTtest() {
        boolean flag = true;
        for (final int size : groupSizes.values()) {
            if (size < 2) {
                flag = false;
                System.out.println("Insufficient data for t-test: need at least 2 samples per group.");
            }
        }
        return flag;
    }

    public DifferentialExpressionResults evaluateDifferentialExpressionStatistics(
            final DifferentialExpressionCalculator deCalculator, final boolean doComparison,
            final ObjectArraySet<NormalizationMethod> normalizationMethods) {
        DifferentialExpressionResults results = null;

        if (!doComparison) {
            deCalculator.createDefaultGroup();
            // make groupComparison be all samples
            groupComparison = new String[]{"all-samples"};
        }
        results = null;
        deCalculator.setRunInParallel(runInParalell);
        LOG.info("Evaluating statistics..");
        if (eval("raw-counts")) {
            results = deCalculator.compare(results, null, new CountRawSampleIdsCalculator());
        }
        for (final NormalizationMethod method : normalizationMethods) {
            method.normalize(deCalculator, getGroupComparison());

            // evaluate per-sample statistics:

            if (eval("samples")) {
                results = deCalculator.compare(results, method, new SampleCountCalculator());
            }

            if (eval("counts")) {
                results = deCalculator.compare(results, method, new CountCalculator());
            }

            if (doComparison) {
                // evaluate differences between groups:
                if (eval("fold-change")) {
                    results = deCalculator.compare(results, method, new FoldChangeCalculator(), groupComparison);
                }
                //results.setOmitNonInformativeColumns(omitNonInformativeColumns);
                if (eval("fold-change-magnitude")) {
                    results = deCalculator.compare(results, method, new FoldChangeMagnitudeCalculator(), groupComparison);
                }
                if (eval("log2-fold-change")) {
                    results = deCalculator.compare(results, method, new Log2FoldChangeCalculator(), groupComparison);
                }

                if (eval("group-averages")) {
                    results = deCalculator.compare(results, method, new AverageCalculator(), groupComparison);
                }

                ttestflag = checkTtest();
                if (ttestflag) {
                    if (eval("t-test")) {
                        results = deCalculator.compare(results, method, new TTestCalculator(), groupComparison);
                    }
                }
                if (eval("fisher")) {
                    results = deCalculator.compare(results, method, new FisherExactTestCalculator(), groupComparison);
                }
                if (eval("fisher-r")) {
                    results = deCalculator.compare(results, method, new FisherExactRCalculator(), groupComparison);
                }
                if (eval("chi-square")) {
                    results = deCalculator.compare(results, method, new ChiSquareTestCalculator(), groupComparison);
                }
            }
            final BenjaminiHochbergAdjustment benjaminiHochbergAdjustment = new BenjaminiHochbergAdjustment();
            final BonferroniAdjustment bonferroniAdjustment = new BonferroniAdjustment();

            if (eval("Bonferroni")) {
                results = bonferroniAdjustment.adjust(results, method, "t-test", "fisher-exact-test", "fisher-exact-R", "chi-square-test");
            }
            if (eval("BH")) {
                results = benjaminiHochbergAdjustment.adjust(results, method, "t-test", "fisher-exact-test", "fisher-exact-R", "chi-square-test");
            }

        }

        return results;
    }


    public boolean eval(final String evalName) {
        return evalSet.contains(evalName.toLowerCase());
    }

    private ObjectSet<String> evalSet;

    public void setEvalNames(final ObjectSet<String> evalSet) {
        this.evalSet = evalSet;
    }

    public void setRunInParallel(boolean parallel) {
        this.runInParalell = parallel;
    }

    public String[] getGroupComparison() {
        String[] result = new String[2];
        GroupComparison val = groupComparisons.get(0);
        result[0] = val.nameGroup1;
        result[1] = val.nameGroup1;
        return result;
    }
}

