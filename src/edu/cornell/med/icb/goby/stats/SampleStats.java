/*
 * Copyright (C) 2009-2012 Institute for Computational Biomedicine,
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

/**
 * @author Fabien Campagne
 *         Date: 1/9/12
 *         Time: 2:08 PM
 */

import edu.cornell.med.icb.goby.modes.VCFCompareMode;
import edu.cornell.med.icb.goby.util.LongNamedCounter;
import edu.cornell.med.icb.goby.util.NamedCounters;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.objects.*;

import java.io.PrintWriter;
import java.util.Collections;

/**
 * Keeps information about sample variation statistics. Used by vcf-compare.
 *
 * @author Fabien Campagne
 *         Date: 1/9/12
 *         Time: 2:00 PM
 */
public class SampleStats {
    public int numGenotypeDisagreements;
    public int numMissedVariantCalls;
    public int numGenotypeAgreements;
    public String sampleId;
    /**
     * Matrix of counters. One counter for each type of error, times the number of files.
     */
    NamedCounters namedCounters;
    public int numHadDifferentAllele;
    private int numFiles;

    public SampleStats(int numFiles) {
        this.numFiles = numFiles;
        transitionCount = new int[numFiles];
        transversionCount = new int[numFiles];
        namedCounters = new NamedCounters(numFiles);
        namedCounters.register("numGenotypeNotInFile", numFiles);
        namedCounters.register("missedOneAllele", numFiles);
        namedCounters.register("missedTwoAlleles", numFiles);
        namedCounters.register("missedMoreThanTwoAlleles", numFiles);
        namedCounters.register("numGenotypeDisagreements", 1);
        namedCounters.register("numGenotypeAgreements", 1);
        namedCounters.register("otherDifferencesInGenotype", 1);

    }

    public NamedCounters counters() {
        return namedCounters;
    }


    @Override
    public String toString() {
        final StringBuffer titvBuffer = new StringBuffer();
        for (int i = 0; i < numFiles; i++) {
            if (i != 0) {
                titvBuffer.append("\t");
            }
            titvBuffer.append(getTransversionToTransitionRatio(i));
        }
        return String.format("STATS\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s%n", sampleId, numGenotypeAgreements,
                namedCounters.get("numGenotypeDisagreements", 0).getCount(),
                numMissedVariantCalls,
                formatSamples(namedCounters.getArray("numGenotypeNotInFile")),
                formatSamples(namedCounters.getArray("missedOneAlleles")),
                formatSamples(namedCounters.getArray("missedTwoAlleles")),
                formatSamples(namedCounters.getArray("missedMoreThanTwoAlleles")),
                numHadDifferentAllele, titvBuffer
        );
    }


    private double fraction(double a, double b) {
        return (a / (b)) * 100;
    }

    private String formatSamples(LongNamedCounter[] counters) {
        StringBuffer buffer = new StringBuffer();
        int i = 0;
        for (int val : LongNamedCounter.valuesInt(counters)) {
            if (i != 0) {
                buffer.append("\t");
            }
            buffer.append(val);
            i++;
        }
        return buffer.toString();
    }


    public void analyze(ObjectSet<String> distinctGenotypes, ObjectList<String> sampleGenotypes, VCFCompareMode.VCFPosition position) {
        int numAllelesAgreed = 0;
        int numAllelesMissed = 0;
        int numAlleleDifference = 0;
        ObjectIterator<String> iterator = distinctGenotypes.iterator();
        String first = iterator.next();
        String second = iterator.next();
        String[] tokenFirst = first.split("/");
        String[] tokenSecond = second.split("/");
        int fileIndex = -1;
        for (int i = 0; i < Math.min(tokenFirst.length, tokenSecond.length); i++) {
            if (tokenFirst[i].equals(tokenSecond[i])) {
                numAllelesAgreed++;
            } else {
                //  System.out.println(i+" differ: "+tokenFirst[i] +" "+isRef(tokenSecond[i]));
                if (isRef(tokenFirst[i]) || isRef(tokenSecond[i])) {
                    numAllelesMissed++;
                    fileIndex = sampleGenotypes.indexOf("ref/ref/");
                } else {
                    numAlleleDifference++;
                }
            }
        }
        if (tokenFirst.length != tokenSecond.length) {
            numAllelesMissed += Math.max(tokenFirst.length, tokenSecond.length) - Math.min(tokenFirst.length, tokenSecond.length);
        }
        if (numAllelesMissed != 0) {
            fileIndex = getIndexWithRef(sampleGenotypes);
            assert fileIndex >= 0 : "fileIndex must have been set";

            if (numAllelesMissed == 2) {
                namedCounters.get("missedTwoAlleles", fileIndex).increment(position);
            }
            if (numAllelesMissed == 1) {
                namedCounters.get("missedOneAllele", fileIndex).increment(position);
            }
            if (numAllelesMissed > 2) {
                namedCounters.get("missedMoreThanTwoAlleles", fileIndex).increment(position);
            }
            if (numAlleleDifference != 0) {
                System.out.println(sampleGenotypes);
                namedCounters.get("otherDifferencesInGenotype", 0).increment(position);
            }
        } else {
            System.out.println(sampleGenotypes);
            namedCounters.get("otherDifferencesInGenotype", 0).increment(position);

        }
    }

    private int getIndexWithRef(ObjectList<String> sampleGenotypes) {
        int sampleIndex = 0;
        int indexhasRef = -1;
        for (String value : sampleGenotypes) {
            ObjectSet<String> unique = new ObjectOpenHashSet<String>();
            for (String token : value.split("/")) {
                unique.add(token);
            }
            if (unique.size() == 1 && "ref".equals(unique.iterator().next())) {
                return sampleIndex;
            } else if (unique.size() >= 1 && value.contains("ref")) {
                indexhasRef = sampleIndex;
            }
            sampleIndex++;
        }
        if (indexhasRef != -1) {
            return indexhasRef;
        }
        return -1;
    }

    private boolean isRef(String allele) {
        return "ref".equals(allele);
    }

    int transitionCount[];
    int transversionCount[];

    public void observeTransitionToTransversions(int fileIndex, ObjectList<String> sampleGenotypes, String ref) {
        String genotype = sampleGenotypes.get(fileIndex);
        int isTransition;
        int isTransversion;
        String alleles[] = genotype.split("/");
        String aBase = "A";
        String gBase = "G";
        transitionCount[fileIndex] += variantWithBases(ref, alleles, aBase, gBase) ? 1 : 0;
        String cBase = "C";
        String tBase = "T";
        // two possible transitions:
        transitionCount[fileIndex] += variantWithBases(ref, alleles, cBase, tBase) ? 1 : 0;
        transversionCount[fileIndex] += variantWithBases(ref, alleles, aBase, cBase) ? 1 : 0;
        // four possible transversions:
        transversionCount[fileIndex] += variantWithBases(ref, alleles, gBase, tBase) ? 1 : 0;
        transversionCount[fileIndex] += variantWithBases(ref, alleles, aBase, tBase) ? 1 : 0;
        transversionCount[fileIndex] += variantWithBases(ref, alleles, cBase, gBase) ? 1 : 0;
    }

    private boolean variantWithBases(String ref, String[] alleles, String base1, String base2) {
        return base1.equals(ref) && contains(alleles, base2, ref) || base2.equals(ref) && contains(alleles, base1, ref);
    }

    public float getTransversionToTransitionRatio(int fileIndex) {
        return ratio(transitionCount[fileIndex], transversionCount[fileIndex]);
    }

    private float ratio(int a, int b) {
        return ((float) a) / ((float) b);
    }

    private boolean contains(String[] alleles, String allele, String ref) {
        for (String v : alleles) {
            if ("ref".equals(v)) {
                //ignore reference allele since this is not a variant.
                continue;
            }
            if (v.indexOf(allele) >= 0) return true;
        }
        return false;
    }

    /**
     * Return text describing the random sample of positions for each counter.
     *
     * @param fileIndex
     * @param reverseIdentifiers
     * @return
     */
    public String toStringExamples(int fileIndex, DoubleIndexedIdentifier reverseIdentifiers) {
        StringBuffer randomSampleText = new StringBuffer();
        String counterNamesByFile[] = {
                "numGenotypeNotInFile", "missedOneAllele", "missedTwoAlleles", "missedMoreThanTwoAlleles", "missedTwoAlleles",
        };
        for (String counterName : counterNamesByFile) {

            printRandomSampleForOneCounter(fileIndex, reverseIdentifiers, randomSampleText, counterName);
        }
        String counterNamesConstant[] = {
                "numGenotypeDisagreements", "otherDifferencesInGenotype"
        };
        for (String counterName : counterNamesConstant) {


            printRandomSampleForOneCounter(0, reverseIdentifiers, randomSampleText, counterName);
        }
        return randomSampleText.toString();
    }

    private void printRandomSampleForOneCounter(int fileIndex, DoubleIndexedIdentifier reverseIdentifiers, StringBuffer randomSampleText, String counterName) {
        LongNamedCounter counter = namedCounters.get(counterName, fileIndex);
        ObjectArrayList<Object> randomSample = counter.getRandomSample();
        long count = counter.getCount();
        if (count != 0) {
            randomSampleText.append("File " + fileIndex + ", out of a total of " + count + ", random sample of positions for counter " + counterName + " (chromosome tab position): \n");
            for (Object o : randomSample) {
                VCFCompareMode.VCFPosition pos = (VCFCompareMode.VCFPosition) o;
                randomSampleText.append(pos.toString(reverseIdentifiers));
                randomSampleText.append("\n");
            }
        }
    }

    public String toStringSample(int fileIndex) {
        String a =
                String.format("COUNT_STATS\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g%n",
                        fileIndex, sampleId,
                        namedCounters.get("numGenotypeAgreements", 0).getCount(),
                        namedCounters.get("numGenotypeDisagreements", 0).getCount(),
                        namedCounters.get("numGenotypeNotInFile", fileIndex).getCount(),
                        namedCounters.get("missedOneAllele", fileIndex).getCount(),
                        namedCounters.get("missedTwoAlleles", fileIndex).getCount(),
                        namedCounters.get("missedMoreThanTwoAlleles", fileIndex).getCount(),
                        namedCounters.get("otherDifferencesInGenotype", 0).getCount(),
                        getTransversionToTransitionRatio(fileIndex)
                );
        int sumErrors = 0;


        sumErrors += namedCounters.get("numGenotypeNotInFile", fileIndex).getCount();
        sumErrors += namedCounters.get("missedOneAllele", fileIndex).getCount();
        sumErrors += namedCounters.get("missedTwoAlleles", fileIndex).getCount();
        sumErrors += namedCounters.get("missedMoreThanTwoAlleles", fileIndex).getCount();
        sumErrors += namedCounters.get("otherDifferencesInGenotype", 0).getCount();

        long disagreementDenominator = sumErrors;
        long observedDenominator = namedCounters.get("numGenotypeDisagreements", 0).getCount() +
                namedCounters.get("numGenotypeAgreements", 0).getCount();

        String b = String.format("FREQ_STATS\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g%n",
                fileIndex, sampleId,
                fraction(namedCounters.get("numGenotypeAgreements", 0).getCount(), observedDenominator),
                fraction(namedCounters.get("numGenotypeDisagreements", 0).getCount(), observedDenominator),
                fraction(namedCounters.get("numGenotypeNotInFile", fileIndex).getCount(), disagreementDenominator),
                fraction(namedCounters.get("missedOneAllele", fileIndex).getCount(), disagreementDenominator),
                fraction(namedCounters.get("missedTwoAlleles", fileIndex).getCount(), disagreementDenominator),
                fraction(namedCounters.get("missedMoreThanTwoAlleles", fileIndex).getCount(), disagreementDenominator),
                fraction(namedCounters.get("otherDifferencesInGenotype", 0).getCount(), disagreementDenominator),
                getTransversionToTransitionRatio(fileIndex)
        );
        return a + b;
    }
}