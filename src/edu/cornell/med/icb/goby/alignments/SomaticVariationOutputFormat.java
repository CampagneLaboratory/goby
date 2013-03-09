package edu.cornell.med.icb.goby.alignments;

import edu.cornell.med.icb.goby.algorithmic.data.CovariateInfo;
import edu.cornell.med.icb.goby.modes.DiscoverSequenceVariantsMode;
import edu.cornell.med.icb.goby.readers.vcf.ColumnType;
import edu.cornell.med.icb.goby.util.OutputInfo;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;

/**
 * @author Fabien Campagne
 *         Date: 3/9/13
 *         Time: 10:13 AM
 */
public class SomaticVariationOutputFormat extends GenotypesOutputFormat {
    private int somaticPValueIndex;

    public void defineColumns(OutputInfo writer, DiscoverSequenceVariantsMode mode) {
        // define columns for genotype format
        super.defineColumns(writer, mode);
        CovariateInfo covInfo = mode.getCovariateInfo();
        ObjectArraySet<String> somaticSampleIds = covInfo.samplesWithExactCovariate("kind-of-sample", "Somatic");
        // add column(s) for p-values of somatic variation:
        for (String sample : somaticSampleIds) {
            somaticPValueIndex = statsWriter.defineField("INFO",
                    String.format("Somatic-P-value[%s]", sample),
                    1, ColumnType.String,
                    "P-value that the variation is somatic in this particular sample, compared to other germline samples (e.g., germline skin, or mother/father).");
        }

    }
}
