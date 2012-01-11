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

package edu.cornell.med.icb.goby.alignments;

/**
 * A codec that always fails. Useful to determine that codecs are wired-in.
 *  @author Fabien Campagne
 *         Date: 1/11/12
 *         Time: 5:02 PM
 */
public class AlignmentCodecDummy implements AlignmentCodec {
    @Override
    public Alignments.AlignmentEntry.Builder encode(Alignments.AlignmentEntry.Builder source) {
      throw new RuntimeException("bug!");
    }

    @Override
    public Alignments.AlignmentEntry.Builder decode(Alignments.AlignmentEntry source) {
       throw new RuntimeException("bug!");

    }

    @Override
    public String name() {
        return "bug";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void newChunk() {
      throw new RuntimeException("bug!");
    }

    @Override
    public byte registrationCode() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
