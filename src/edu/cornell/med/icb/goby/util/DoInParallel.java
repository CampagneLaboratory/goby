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
    private int currentlyActive;

    protected DoInParallel() {
        super();
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

    protected void debugEnd(final String inputBasename) {
        synchronized (this) {
            currentlyActive--;
            LOG.debug("                                               Stopped " + inputBasename);
        }
    }

    protected void debugStart(final String inputBasename) {
        LOG.debug("                                               Starting " + inputBasename);
        synchronized (this) {
            currentlyActive++;
            LOG.debug("== Currently active: " + currentlyActive);
            System.out.flush();
        }
    }

    public void execute(final boolean parallel, final String[] inputFilenames) throws Exception {
        final BasenameParallelRegion region = new BasenameParallelRegion(this, inputFilenames);
        this.parallel = parallel;
        getParallelTeam().execute(region);
    }

    public static void main(final String[] args) throws Exception {
        final DoInParallel loop = new DoInParallel() {

            @Override
            public void action(final DoInParallel forDataAccess, final String inputBasename, final int loopIndex) {
                try {
                    debugStart(inputBasename);

                    Thread.sleep(10000);
                    debugEnd(inputBasename);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        };
        final String[] inputs = new String[100];
        for (int i = 0; i < inputs.length; i++) {

            inputs[i] = Integer.toString(i + 1);
        }
        loop.execute(true, inputs);
    }
}
