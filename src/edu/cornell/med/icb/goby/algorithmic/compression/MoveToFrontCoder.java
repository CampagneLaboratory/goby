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

package edu.cornell.med.icb.goby.algorithmic.compression;

/**
 * @author Fabien Campagne
 *         Date: 1/15/12
 *         Time: 11:02 AM
 */
public class MoveToFrontCoder {
    int[] symbols;
    private int numSymbols;

    public MoveToFrontCoder(int numberOfSymbols) {
        numSymbols = numberOfSymbols;

        symbols = new int[numberOfSymbols ];
        for (int i=0;i<numberOfSymbols;i++) {
            symbols[i]=i;
        }

    }

    /**
     * Return the encoded index of the symbol.
     * @param symbol Symbol to encode with move to front (precondition: 0<= symbol <numberOfSymbols).
     * @return encoded symbol (post-condition: 0<= encoded-symbol <numberOfSymbols).
     */
    public int encode(int symbol) {
        if (symbol < 0) {
            throw new IllegalArgumentException("symbol must be positive.");
        }
        if (symbol > numSymbols) {
            throw new IllegalArgumentException(String.format("You have only %d symbols, symbol value is too large %d.",
                    numSymbols, symbol));
        }
      try {
        int carry = -1;
        for (int i = 0; i < numSymbols; i++) {

            if (symbols[i] != symbol) {
                int tmp=symbols[i];
                if (carry!=-1) {
                    symbols[i ] = carry;
                }
                carry = tmp;
            } else {
                if (carry!=-1) {
                    symbols[i]=carry;
                }
                symbols[0] = symbol;
                return i;
            }
        }
      } finally {
       //   System.out.printf("encode symbol %d then: %s%n", symbol, IntArrayList.wrap(symbols));
      }
        throw new IllegalArgumentException(String.format("Symbol %d must be found.",
                    symbol));
    }
}
