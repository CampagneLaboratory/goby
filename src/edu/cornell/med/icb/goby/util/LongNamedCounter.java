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

package edu.cornell.med.icb.goby.util;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.Random;

/**
 * A counter to keep track of the number of objects of a certain kind. Each counter can collect
 * a random sample of objects that contribute to the count.
 *
 * @author Fabien Campagne
 *         Date: 1/13/12
 *         Time: 2:37 PM
 */
public class LongNamedCounter {
    public String getName() {
        return name;
    }

    private String name;
    private int index;
    private long count;
    private static final int MAX_RANDOM = 5;

    public LongNamedCounter(final String name) {
        this.name = name;
    }

    public void increment() {
        count += 1;
    }

    public void increment(final Object value) {
        count += 1;
        randomSample(value);

    }

    Random randomEngine = new Random();
    private ObjectArrayList<Object> randomSample = new ObjectArrayList<Object>();

    public ObjectArrayList<Object> getRandomSample() {
        return randomSample;
    }

    public void randomSample(Object value) {
        if (randomSample.size() < MAX_RANDOM) {
            randomSample.add(value);
        } else {

            double v = randomEngine.nextDouble();
            // as more elements are seen, make it increasingly hard to replace an element seen earlier:

            double threshold = 1.0 / (double) Math.log(count);
            if (v < threshold) {
                double l=v/threshold;
                final int index = (int) Math.round(l  * ((double)randomSample.size()-1f));
              //  System.out.printf("l=%g index=%d%n",l, index);
                // replace value at random position.
                randomSample.set(index, value);
            }
        }
    }

    public static long[] valuesLong(LongNamedCounter[] counters) {
        long[] result = new long[counters.length];
        int i = 0;
        for (LongNamedCounter counter : counters) {
            result[i++] = counter.getCount();
        }
        return result;
    }
     @Override
    public String toString() {
        return String.format("[%s: %d]%n", name, count);
    }

    public static int[] valuesInt(LongNamedCounter[] counters) {
        int[] result = new int[counters.length];
        int i = 0;
        for (LongNamedCounter counter : counters) {
            result[i++] = (int) counter.getCount();
        }
        return result;
    }

    public long getCount() {
        return count;
    }
}
