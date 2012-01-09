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

    public static String header() {
        return "sampleId,numGenotypeAgreements, numGenotypeDisagreements,numMissedVariantCalls,numGenotypeNotCalled,missedOneAlleles,missedTwoAlleles,missedMoreThanTwoAlleles,numHadDifferentAllele\n".replace(',', '\t');
    }

    @Override
    public String toString() {
        return String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n", sampleId, numGenotypeAgreements, numGenotypeDisagreements,
                numMissedVariantCalls, numGenotypeNotCalled,
                missedOneAlleles, missedTwoAlleles, missedMoreThanTwoAlleles,numHadDifferentAllele
        );
    }


    public void analyze(ObjectSet<String> distinctGenotypes) {
        int numAllelesAgreed = 0;
        int numAllelesMissed = 0;
        int numAlleleDifference=0;
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
            if (numAlleleDifference!=0) {
                numHadDifferentAllele++;
            }
        }
    }

    private boolean isRef(String allele) {
        return "ref".equals(allele);
    }
}