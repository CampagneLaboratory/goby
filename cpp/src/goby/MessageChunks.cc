/**
 * Copyright (C) 2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>
#include "Alignments.pb.h"
#include "MessageChunks.h"

using namespace std;

namespace goby {
  MessageChunksReader::MessageChunksReader(const string& filename) : DELIMITER_LENGTH(8) {
    this->filename = filename;
    this->currentIndex = 0;
    this->currentCollection = AlignmentCollection::default_instance();

    // get the positions each of the chunks in the file
    ifstream stream;
    stream.open(filename.c_str(), ios::in | ios::binary);
    while (stream.good()) {
      // each chunk is delimited by DELIMITER_LENGTH bytes
      stream.seekg(DELIMITER_LENGTH, ios::cur);

      // then the size of the next chunk follows
      int size = readInt(stream);
      cout << "size is " << size << endl;

      // the last chunk has a size of zero bytes
      if (!stream.eof() && size != 0) {
        const streampos position = stream.tellg();
        this->positions.push_back(position);
        this->lengths.push_back(size);
        stream.seekg(size, ios::cur);
      } else {
        break;
      }
    }
    stream.close();
  }
  
  MessageChunksReader::~MessageChunksReader(void) {
  }

  // Java DataInput.readInt()
  int MessageChunksReader::readInt(std::istream& stream) {
    const int ch1 = stream.get();
    const int ch2 = stream.get();
    const int ch3 = stream.get();
    const int ch4 = stream.get();
    return (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);
  }
}