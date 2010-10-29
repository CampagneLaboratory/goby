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

package edu.cornell.med.icb.goby.util;

import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;

/**
 * @author Fabien Campagne
 *         Date: Mar 26, 2010
 *         Time: 5:25:57 PM
 */
class BasenameParallelRegion extends ParallelRegion {
    private final String[] inputFilenames;
    private final DoInParallel loop;

    BasenameParallelRegion(final DoInParallel loop, final String[] inputFilenames) {
        super();
        this.loop = loop;
        this.inputFilenames = inputFilenames;

    }

    @Override
    public void run() throws Exception {
        execute(0, inputFilenames.length - 1 /* end index must be inclusive. This is counter-intuitive */, new IntegerForLoop() {

            @Override
            public void run(final int startIndex, final int endIndex) {
                //   System.out.println(String.format("executing start= %d end=%d ",startIndex, endIndex));
                for (int i = startIndex; i <= endIndex; ++i) {
                    if (i >= 0 && i < inputFilenames.length) {

                        final String inputBasename = AlignmentReader.getBasename(inputFilenames[i]);
                        loop.action(loop, inputBasename, i);

                    }
                }
            }
        });
    }
}
