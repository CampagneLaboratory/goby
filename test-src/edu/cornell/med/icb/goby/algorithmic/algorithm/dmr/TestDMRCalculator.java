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

package edu.cornell.med.icb.goby.algorithmic.algorithm.dmr;

import com.sun.tools.javac.util.Position;
import com.sun.xml.internal.bind.marshaller.MinimumEscapeHandler;
import edu.cornell.med.icb.goby.algorithmic.data.MethylRateInfo;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.easymock.internal.matchers.Null;
import org.junit.Test;

import java.util.List;

/**
 * User: nyasha
 * Date: 6/14/11
 * Time: 4:03 PM
 */
public class TestDMRCalculator {

    @Test

    public void testCase1() {

        DMRCalculator calculator = new DMRCalculator();

        int position = 0;

        ObjectArrayList<MethylRateInfo> ratesAtPosition = new ObjectArrayList<MethylRateInfo>();

        calculator.observe(position, ratesAtPosition);
        calculator.observe(1, new ObjectArrayList<MethylRateInfo>());
        calculator.observe(2, new ObjectArrayList<MethylRateInfo>());

        // sampleIndex:groupIndex pos,methylationRate
        String m = "0:0 1,0.9 2,0.9 3,0.9 4,0.9 5,0.9";
        String m1 = "1:1 1,0.1 2,0.1 3,0.1 4,0.1 5,0.1";

        String[] data = {m, m1};

        ObjectArrayList<MethylRateInfo> result = getPositionData(1, data);



    }

    //take  string and convert to methylrateinfo list
    private ObjectArrayList<MethylRateInfo> getPositionData(int position, String[] data) {
        int numSamples = data.length;

        ObjectArrayList<MethylRateInfo> positionList = new ObjectArrayList<MethylRateInfo>();
        int sampleIndex, groupIndex, genomicPosition;
        float methylationRate;
        String toParse;
        String[] tokens, sampleInfo, methylInfo;
        for (int x = 0; x < numSamples; x++) {
            toParse = data[x];
            tokens = toParse.split("\\s");
            sampleInfo = tokens[0].split(":");
            methylInfo = tokens[position].split(",");
            sampleIndex = Integer.parseInt(sampleInfo[0]);
            groupIndex = Integer.parseInt(sampleInfo[1]);
            genomicPosition = Integer.parseInt(methylInfo[0]);
            methylationRate = Float.valueOf(methylInfo[1]);
            MethylRateInfo addToList = new MethylRateInfo(sampleIndex, groupIndex, genomicPosition, methylationRate);
            positionList.add(addToList);
        }
        return positionList;  
    }


}
