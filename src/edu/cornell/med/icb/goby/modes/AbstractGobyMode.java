/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

package edu.cornell.med.icb.goby.modes;

/**
 * Specifies the maqsupport jar file for all goby modes.
 *
 * @author Fabien Campagne
 * Date: Apr 7, 2009
 * Time: 3:14:02 PM
 */
public abstract class AbstractGobyMode extends AbstractCommandLineMode {
    protected AbstractGobyMode(final String jarFilename) {
        super(jarFilename);
    }

    protected AbstractGobyMode() {
        super("goby.jar");
    }
}
