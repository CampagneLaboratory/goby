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

import it.unimi.dsi.fastutil.objects.ObjectIterator;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectSet;

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
    public int numGenotypeNotCalled;
    public int missedTwoAlleles;
    public int missedOneAlleles;
    public int missedMoreThanTwoAlleles;
    public int numHadDifferentAllele;
    private int numFiles;

    public SampleStats(int numFiles) {
        this.numFiles = numFiles;
        transitionCount = new int[numFiles];
        transversionCount = new int[numFiles];
    }

    public static String header(int numFiles) {
        final StringBuffer titvheaders = new StringBuffer();
        for (int i = 0; i < numFiles; i++) {
            if (i != 0) {
                titvheaders.append("\t");
            }
            titvheaders.append("ti/tv_ratio_file_");
            titvheaders.append(i);
        }

        return "sampleId,numGenotypeAgreements, numGenotypeDisagreements,numMissedVariantCalls,numGenotypeNotCalled,missedOneAlleles,missedTwoAlleles,missedMoreThanTwoAlleles,numHadDifferentAllele".replace(',', '\t') + titvheaders + "\n";
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
        return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s%n", sampleId, numGenotypeAgreements, numGenotypeDisagreements,
                numMissedVariantCalls, numGenotypeNotCalled,
                missedOneAlleles, missedTwoAlleles, missedMoreThanTwoAlleles, numHadDifferentAllele, titvBuffer
        );

    }


    public void analyze(ObjectSet<String> distinctGenotypes, ObjectList<String> sampleGenotypes) {
        int numAllelesAgreed = 0;
        int numAllelesMissed = 0;
        int numAlleleDifference = 0;
        ObjectIterator<String> iterator = distinctGenotypes.iterator();
        String first = iterator.next();
        String second = iterator.next();
        String[] tokenFirst = first.split("/");
        String[] tokenSecond = second.split("/");

        for (int i = 0; i < Math.min(tokenFirst.length, tokenSecond.length); i++) {
            if (tokenFirst[i].equals(tokenSecond[i])) {
                numAllelesAgreed++;
            } else {
                //  System.out.println(i+" differ: "+tokenFirst[i] +" "+isRef(tokenSecond[i]));
                if (isRef(tokenFirst[i]) || isRef(tokenSecond[i])) {
                    numAllelesMissed++;
                } else {
                    numAlleleDifference++;
                }
            }
        }
        if (tokenFirst.length != tokenSecond.length) {
            numAllelesMissed += Math.max(tokenFirst.length, tokenSecond.length) - Math.min(tokenFirst.length, tokenSecond.length);
        }
        if (numAllelesMissed != 0) {

            if (numAllelesMissed == 2) {
                missedTwoAlleles += 1;
            }
            if (numAllelesMissed == 1) {
                missedOneAlleles += 1;
            }
            if (numAllelesMissed > 2) {
                missedMoreThanTwoAlleles += 1;
            }
            if (numAlleleDifference != 0) {
                numHadDifferentAllele++;
            }
        }
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
}