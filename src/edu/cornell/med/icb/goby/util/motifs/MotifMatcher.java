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

package edu.cornell.med.icb.goby.util.motifs;

/**
 * @author campagne
 *         Date: 9/28/11
 *         Time: 12:27 PM
 */
public interface MotifMatcher {
    /**
     * Accept a new character.
     * @param c character to accept.
     */
    public void accept(char c);

    /**
     * Indicate the start of a new sequence.
     */
    public void newSequence();
    /**
     * Determine if the motif was matched up to the last character accepted.
     * @return
     */
    public boolean matched();

}
