/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.methylation;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;


/**
 *  Stores information about differentially methylated regions (DMRs) across groups
 *
 * @author Nyasha Chambwe
 * Date: 10/3/11
 * Time: 12:07 PM
 */
public class DifferentiallyMethylatedRegion {

    public String chromosome;
    public int start;
    public int end;
    public String strand;
    public ObjectArrayList<MethylationRegion> [] methylatedRegionPerSample;
    public double [] meanMethylationRatePerSample;
    public double foldChange;
    public double pValue;
    public double qValue;


    public DifferentiallyMethylatedRegion(String chromosome, int start, int end, String strand) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    public DifferentiallyMethylatedRegion(String chromosome, int start, int end, String strand, ObjectArrayList<MethylationRegion>[] methylatedRegionPerSample) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.methylatedRegionPerSample = methylatedRegionPerSample;
    }
}
