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
#define GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH 8    // length of the delimiter tag (in bytes)

  struct MessageChunk {
    std::streampos position;
    size_t length;
  };

  template <class T> class LIBGOBY_EXPORT MessageChunksReader : public std::iterator<std::input_iterator_tag, T> {
    // the name of the chunked file
    std::string filename;

    // the underlying stream
    std::ifstream *stream;

    // positions of compressed alignment collections within the entries file
    std::vector<MessageChunk> chunks;

    // index of current chunk
    unsigned currentChunk;

    // the current processed chunk
    T current;

    // Java DataInput.readInt()
    static int readInt(std::istream &stream) {
      // TODO? Do we need to worry about endian order here?
      const int ch1 = stream.get();
      const int ch2 = stream.get();
      const int ch3 = stream.get();
      const int ch4 = stream.get();
      return (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);
    };

  public:
    MessageChunksReader(const std::string& filename) {
      this->filename = filename;
      this->currentChunk = 0;
      this->current = T::default_instance();
      this->stream = new std::ifstream(filename.c_str(), std::ios::in | std::ios::binary);

      int chunkNumber = 0;
      // get the positions each of the chunks in the file
      while (stream->good()) {
        // each chunk is delimited by DELIMITER_LENGTH bytes
        stream->seekg(GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH, std::ios::cur);

        // then the size of the next chunk follows
        const int size = readInt(*stream);
#if GOBY_DEBUG
        std::cout << "length of chunk #" << chunkNumber++ << " is " << size << std::endl;
#endif // GOBY_DEBUG

        // the last chunk has a size of zero bytes
        if (!stream->eof() && size != 0) {
          const std::streampos position = stream->tellg();
          MessageChunk chunk;
          chunk.position = position;
          chunk.length = size;
          chunks.push_back(chunk);
          stream->seekg(size, std::ios::cur);
        } else {
          break;
        }
      }
    };

    virtual ~MessageChunksReader(void) {
      stream->close();
      delete stream;
    };

    // TODO MessageChunksReader(const MessageChunksReader& reader);

    // Prefix increment operator
    MessageChunksReader& operator++() {
      std::cout << "Prefix operator++() " << this->currentChunk << std::endl;
      ++currentChunk;
      return *this;
    };

    // Postfix increment operator
    MessageChunksReader& operator++(int) {
      std::cout << "Postfix operator++(int) " << this->currentChunk << std::endl;
      currentChunk++;
      return *this;
    };

    bool operator==(const MessageChunksReader& rhs) {
      return currentChunk == rhs.currentChunk;  // TODO
    };

    bool operator!=(const MessageChunksReader& rhs) {
      return currentChunk != rhs.currentChunk;  // TODO
    };

    T& operator*() {
#if GOBY_DEBUG
        std::cout << "chunk number is " << currentChunk << std::endl;
#endif // GOBY_DEBUG
      return current;
    };

    friend std::ostream &operator<<(std::ostream &out, const MessageChunksReader& reader) {
      out << "ostream &operator<< " << reader.currentChunk;
      return out;
    };

    MessageChunksReader begin() {
      return(MessageChunksReader(NULL)); // TODO
    };

    MessageChunksReader end() {
      return(MessageChunksReader(NULL)); // TODO
    };
  };
}

#endif // GOBY_MESSAGE_CHUNKS_H
