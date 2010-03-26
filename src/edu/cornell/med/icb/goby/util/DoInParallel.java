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

package edu.cornell.med.icb.goby.util;

import edu.rit.pj.ParallelTeam;
import org.apache.log4j.Logger;

/**
 * Helper class to support parallel processing of a set of files.
 *
 * @author Fabien Campagne
 *         Date: Mar 26, 2010
 *         Time: 1:55:21 PM
 */
public abstract class DoInParallel {
    private ParallelTeam team;
    private boolean parallel;
    protected static final Logger LOG = Logger.getLogger(DoInParallel.class);

    public DoInParallel() {
    }

    protected synchronized ParallelTeam getParallelTeam() {
        if (team == null) {
            if (!parallel) {
                // 1 thread only is sequential
                team = new ParallelTeam(1);
            } else {
                // as many threads as configured with -Dpj.nt or default.
                team = new ParallelTeam();
            }
        }
        LOG.info("Executing on " + team.getThreadCount() + " threads.");
        return team;
    }

   public abstract void action(DoInParallel forDataAccess, String inputBasename, int loopIndex);

    public void execute(boolean parallel, String[] inputFilenames) throws Exception {

        final BasenameParallelRegion region = new BasenameParallelRegion(this,inputFilenames);
        this.parallel=parallel;
        getParallelTeam().execute(region);


    }
}
