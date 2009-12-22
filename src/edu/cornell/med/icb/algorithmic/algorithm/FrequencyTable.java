/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.algorithmic.algorithm;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import org.apache.commons.io.IOUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class FrequencyTable {
    Int2IntMap a;
    Int2IntMap c;
    Int2IntMap t;
    Int2IntMap g;
    Int2IntMap n; //position to count

    public FrequencyTable() {
        a = new Int2IntOpenHashMap();
        c = new Int2IntOpenHashMap();
        t = new Int2IntOpenHashMap();
        g = new Int2IntOpenHashMap();
        n = new Int2IntOpenHashMap();
        a.defaultReturnValue(0);
        c.defaultReturnValue(0);
        t.defaultReturnValue(0);
        g.defaultReturnValue(0);
        n.defaultReturnValue(0);
    }

    public void put(final int position, final char base) {
        if (base == 'A') {
            int count = a.get(position);
            a.put(position, ++count);
        } else if (base == 'C') {
            int count = c.get(position);
            c.put(position, ++count);
        } else if (base == 'T') {
            int count = t.get(position);
            t.put(position, ++count);
        } else if (base == 'G') {
            int count = g.get(position);
            g.put(position, ++count);
        } else {
            int count = n.get(position);
            n.put(position, ++count);
        }
    }

    public int get(final int position, final char base) {
        final int rval;
        if (base == 'A') {
            rval = a.get(position);
        } else if (base == 'C') {
            rval = c.get(position);
        } else if (base == 'T') {
            rval = t.get(position);
        } else if (base == 'G') {
            rval = g.get(position);
        } else {
            rval = n.get(position);
        }
        return rval;
    }

    public void addAll(final FrequencyTable right) {
        add(this.a, right.a);
        add(this.c, right.c);
        add(this.t, right.t);
        add(this.g, right.g);
        add(this.n, right.n);
    }

    public void add(final Int2IntMap x, final Int2IntMap y) {
        for (final Int2IntMap.Entry entry : y.int2IntEntrySet()) {
            final int key = entry.getKey();
            final int ycount = entry.getValue();
            final int xcount = x.get(key);
            x.put(key, xcount + ycount);
            //     entry.setValue(xcount + ycount);
        }
    }

    public void output(final String filename) throws IOException {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(filename));
            writer.write("\ta\tc\tt\tg\tn\n");
            for (int position = -2; position <= 2; position++) {
                writer.write(position + "\t");
                //Start writing to the output stream
                final int acount = get(position, 'A');
                final int ccount = get(position, 'C');
                final int tcount = get(position, 'T');
                final int gcount = get(position, 'G');
                final int ncount = get(position, 'N');
                writer.write(acount + "\t" + ccount + "\t" + tcount + "\t" + gcount + "\t" + ncount + "\n");
            }
        } finally {
            IOUtils.closeQuietly(writer);
        }
    }
}
