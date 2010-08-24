//
// Copyright (C) 2010 Institute for Computational Biomedicine,
//                    Weill Medical College of Cornell University
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <fstream>
#include <iostream>
#include <string>

#include "common.h"
#include "hash.h"

namespace goby {
  // load values from java style propertes file and store them into the provided map
  // copied from propertyutil.cpp - http://www.senzee5.com/2008/02/java-style-properties-files-in-c.html
  LIBGOBY_EXPORT void load_properties(const std::string& filename, LIBGOBY_HASH_MAP<std::string, std::string>& property_map) {
    std::ifstream is(filename.c_str());
    if (!is) {
      std::cerr << "Failed to open property file: " << filename << std::endl;
    }

    char next = 0;
    char ch = is.get();
    while (!is.eof())
    {
      switch (ch)
      {
      case '#' :
      case '!' : 
        do ch = is.get();
        while (!is.eof() && ch >= 0 && ch != '\n' && ch != '\r');               
        continue;
      case '\n':
      case '\r':
      case ' ' :
      case '\t': ch = is.get(); continue;
      }

      // Read the key
      std::ostringstream key, val;

      while (!is.eof() && ch >= 0 && ch != '=' && ch != ':' &&
        ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r')
      {
        key << ch;
        ch = is.get();
      }

      while (!is.eof() && (ch == ' ' || ch == '\t'))
        ch = is.get();

      if (!is.eof() && (ch == '=' || ch == ':'))
        ch = is.get();

      while (!is.eof() && (ch == ' ' || ch == '\t'))
        ch = is.get();

      // Read the value
      while (!is.eof() && ch >= 0 && ch != '\n' && ch != '\r')
      {
        int next = 0;
        if (ch == '\\')
        {
          ch = is.get();
          switch (ch)
          {
          case '\r':
            ch = is.get();
            if (ch != '\n' && ch != ' ' && ch != '\t')
              continue;
            // fall through
          case '\n':
            ch = is.get();
            while (!is.eof() && (ch == ' ' || ch == '\t')) { is >> ch; }
            continue;
          case 't': ch = '\t'; next = is.get(); break;
          case 'n': ch = '\n'; next = is.get(); break;
          case 'r': ch = '\r'; next = is.get(); break;
          case 'u':
            {
              ch = is.get();
              while (!is.eof() && ch == 'u') { is >> ch; }
              int d = 0;
loop:
              for (int i = 0 ; !is.eof() && i < 4 ; i++)
              {
                next = is.get();
                switch (ch)
                {
                case '0': case '1': case '2': case '3': case '4':
                case '5': case '6': case '7': case '8': case '9':           d = (d << 4) +      ch - '0'; break;
                case 'a': case 'b': case 'c': case 'd': case 'e': case 'f': d = (d << 4) + 10 + ch - 'a'; break;
                case 'A': case 'B': case 'C': case 'D': case 'E': case 'F': d = (d << 4) + 10 + ch - 'A'; break;
                default:                                                                                  goto loop;
                }
                ch = is.get();
              }
              ch = d;
              break;
            }
          default: next = is.get(); break;
          }
        }
        else
          next = is.get();

        val << ch;
        ch = next;
      }
      property_map[key.str()] = val.str();
    }
  }
}
