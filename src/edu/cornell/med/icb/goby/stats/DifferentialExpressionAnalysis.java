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

package edu.cornell.med.icb.goby.stats;

import com.martiansoftware.jsap.JSAPResult;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.modes.CompactAlignmentToAnnotationCountsMode;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.ServiceLoader;

/**
 * Created by IntelliJ IDEA.
 * User: nyasha
 * Date: Mar 24, 2010
 * Time: 3:29:26 PM
 * Class that organizes group comparisons and facillitates differential expression statistical testing
 */
public class DifferentialExpressionAnalysis {

    public DifferentialExpressionAnalysis() {
    }

    /* Used to log debug and informational messages.
         */
    private static final Log LOG = LogFactory.getLog(CompactAlignmentToAnnotationCountsMode.class);

    /**
     * The groups that should be compared, order matters.
     */
    private String[] groupComparison;
    private final ObjectSet<String> groups = new ObjectArraySet<String>();
    private final Object2ObjectMap<String, Integer> groupSizes = new Object2ObjectOpenHashMap<String, Integer>();
    /**
     * Flag for an invalid t-test due to any group size < 2
     *
     * @return
     */
    private boolean ttestflag;
    /**
     * The set of normalization methods to use for the comparison.
     */

    private ObjectArraySet<NormalizationMethod> normalizationMethods;


    private static final ServiceLoader<NormalizationMethod> normalizationMethodLoader
            = ServiceLoader.load(NormalizationMethod.class);

    public void parseGroupsDefinition(final String groupsDefinition, final DifferentialExpressionCalculator deCalculator,
                                      final String[] inputFilenames) {
        if (groupsDefinition == null) {
            // no groups definition to parse.
            return;
        }

        final String[] groupsTmp = groupsDefinition.split("/");
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
            groups.add(groupId);
            for (final String groupString : groupBasenames.split(",")) {

                final String groupBasename = FilenameUtils.getBaseName(AlignmentReader.getBasename(groupString));
                if (!isInputBasename(groupBasename, inputFilenames)) {
                    System.err.printf("The group basename %s is not a valid input basename.%n", groupBasename);
                    System.exit(1);
                }
                System.out.println("Associating basename: " + groupBasename + " to group: " + groupId);
                deCalculator.associateSampleToGroup(groupBasename, groupId);
            }
            final int groupSize = (groupBasenames.split(",")).length;
            groupSizes.put(groupId, groupSize);
        }
    }

    /**
     * Return true if the basename is on the command line as an input basename.
     *
     * @param basename
     * @return
     */
    private boolean isInputBasename(final String basename, final String[] inputFilenames) {
        for (final String inputFilename : inputFilenames) {
            if (FilenameUtils.getBaseName(AlignmentReader.getBasename(inputFilename)).equals(basename)) {
                return true;
            }
        }
        return false;
    }

    public void parseCompare(final String compare) {
        final String[] groupLanguageText = compare.split("/");
        for (final String groupId : groupLanguageText) {
            if (!groups.contains(groupId)) {
                System.err.println("Group " + groupId + " used in --compare must be defined. "
                        + "Please see the --groups option to define groups.");
                System.exit(1);
            }
        }
        groupComparison = groupLanguageText;
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
        if (doComparison) {
            results = null;

            for (final NormalizationMethod method : normalizationMethods) {
                method.normalize(deCalculator, groupComparison);
                // evaluate differences between groups:
                results = deCalculator.compare(results, method, new FoldChangeCalculator(), groupComparison);
                //results.setOmitNonInformativeColumns(omitNonInformativeColumns);
                results = deCalculator.compare(results, method, new FoldChangeMagnitudeCalculator(), groupComparison);
                results = deCalculator.compare(results, method, new Log2FoldChangeCalculator(), groupComparison);

                results = deCalculator.compare(results, method, new AverageCalculator(), groupComparison);
                ttestflag = checkTtest();
                if (ttestflag) {
                    results = deCalculator.compare(results, method, new TTestCalculator(), groupComparison);
                }
                results = deCalculator.compare(results, method, new FisherExactTestCalculator(), groupComparison);
                results = deCalculator.compare(results, method, new FisherExactRCalculator(), groupComparison);
                results = deCalculator.compare(results, method, new ChiSquareTestCalculator(), groupComparison);

                final BenjaminiHochbergAdjustment benjaminiHochbergAdjustment = new BenjaminiHochbergAdjustment();
                final BonferroniAdjustment bonferroniAdjustment = new BonferroniAdjustment();

                results = bonferroniAdjustment.adjust(results, method, "t-test", "fisher-exact-test", "fisher-exact-R", "chi-square-test");
                results = benjaminiHochbergAdjustment.adjust(results, method, "t-test", "fisher-exact-test", "fisher-exact-R", "chi-square-test");
            }
        }
        return results;
    }
}

