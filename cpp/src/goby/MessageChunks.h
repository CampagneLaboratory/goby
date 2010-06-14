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

#pragma once

#ifndef GOBY_MESSAGE_CHUNKS_H
#define GOBY_MESSAGE_CHUNKS_H

#include <fstream>
#include <iostream>
#include <istream>
#include <iterator>
#include <string>
#include <vector>
#include "common.h"
#include "Alignments.pb.h"

namespace goby {
  class LIBGOBY_EXPORT MessageChunksReader : public std::iterator<std::input_iterator_tag, AlignmentCollection> {
    // the name of the chunked file
    std::string filename;

    // positions of compressed alignment collections within the entries file
    std::vector<std::streampos> positions;
    
    // positions of compressed alignment collections within the entries file    
    std::vector<int> lengths;

    unsigned currentIndex;
    AlignmentCollection currentCollection;

    static int readInt(std::istream& stream);

  public:
    // length of the delimiter tag (in bytes)
    const unsigned DELIMITER_LENGTH;
    
    MessageChunksReader(const std::string& filename);
    virtual ~MessageChunksReader(void);

    MessageChunksReader(const MessageChunksReader& reader);
    MessageChunksReader& operator++();
    MessageChunksReader& operator++(int);
    bool operator==(const MessageChunksReader& rhs);
    bool operator!=(const MessageChunksReader& rhs);
    AlignmentCollection& operator*() { return currentCollection; };
  };
}

#endif // GOBY_MESSAGE_CHUNKS_H
