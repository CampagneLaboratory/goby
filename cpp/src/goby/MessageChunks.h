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
#include <iterator>
#include <string>
#include <vector>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include "common.h"
#include "Alignments.pb.h"

namespace goby {
  #define GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH 8      // length of the delimiter tag (in bytes)

  template <typename T> class MessageChunksIterator : public std::iterator<std::input_iterator_tag, T> {
    // the name of the chunked file
    std::string filename;

    // The underlying stream
    std::ifstream stream;

    // The current position in the stream (we keep this since "tellg()" is not const)
    std::streampos current_position;

    // the current processed chunk
    T current_chunk;

    // populate T with the data from the given chunk
    // The assumption here is that the stream is positioned at the start of a chunk boundary
    // The stream will be positioned at the next chunk boundary or eof
    T populateChunk(std::ifstream& stream, T& chunk) {
      if (!stream.good()) {
        current_position = static_cast<std::streampos>(-1);
        chunk.Clear();
        return chunk;
      }

      // each chunk is delimited by DELIMITER_LENGTH bytes
      stream.seekg(GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH, std::ios::cur);
      current_position = stream.tellg();

      // Set up a stream that will only read up the current chunk length
      const size_t chunkLength = readInt(stream);
      current_position = stream.tellg();

#ifdef GOBY_DEBUG
      std::cout << "Chunk length: " << chunkLength << std::endl;
#endif

      // the last chunk has a length of zero bytes
      if (stream.eof() || chunkLength == 0) {
        current_position = static_cast<std::streampos>(-1);
        chunk.Clear();
        return chunk;
      }

      const std::streampos chunkStartPosition = current_position;
      const std::streampos nextChunkStartPosition = current_position + static_cast<std::streamoff>(chunkLength);

      // explicitly set the block size and limit it to the chunk length to prevent "over read"
      google::protobuf::io::IstreamInputStream istream(&stream, chunkLength);
      google::protobuf::io::LimitingInputStream rawChunkStream(&istream, chunkLength);

      // and handle the fact that each chunk is compressed with gzip
      google::protobuf::io::GzipInputStream gzipChunkStream(&rawChunkStream);

      // since the chunks may be large, we may need to increase the limit
      google::protobuf::io::CodedInputStream codedStream(&gzipChunkStream);
      codedStream.SetTotalBytesLimit(chunkLength, 0);

      // populate the current object from the compressed data
      if (!chunk.ParseFromCodedStream(&codedStream)) {
        std::cerr << "Failed to parse message chunk from " << filename << std::endl;
      }

      current_position = stream.tellg();

      // may need to adjust the position in the stream since the parsing above may have "over-read"
      // TODO - can probably remove this block since the block size is set above
      if (current_position != static_cast<std::streampos>(-1)) {
        const std::streamoff offset = nextChunkStartPosition - current_position;
        if (offset != 0) {
          stream.seekg(offset, std::ios::cur);
          current_position = stream.tellg();
          std::cerr << "WARNING: Stream position adjuted to: " << current_position << std::endl;
        }
      }

      // and retrun the processed chunk
      return chunk;
    };

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
    // TOOD: currently assumes that the position is at the start of chunk boundary
    MessageChunksIterator<T>(const std::string& filename, const std::streampos position = 0) :
        filename(filename),
        current_position(position),
        current_chunk(T::default_instance()) {
      stream.open(filename.c_str(), std::ios::in | std::ios::binary);
      stream.seekg(position);

      if (stream.good()) {
       populateChunk(stream, current_chunk);
      } else {
        std::cerr << "Failed to open " << filename << std::endl;
      }
    };

    virtual ~MessageChunksIterator(void) {
      stream.close();
    };

    MessageChunksIterator(const MessageChunksIterator<T>& that) :
        filename(that.filename),
        current_position(that.current_position),
        current_chunk(that.current_chunk) {
      // TODO: there is probably a better way to copy the stream
      stream.open(filename.c_str(), std::ios::in | std::ios::binary);
      stream.seekg(current_position);
      std::cout << "MessageChunksIterator Copy constructor" << std::endl;
    }

    MessageChunksIterator(const MessageChunksIterator<T>& that, std::streamoff off, std::ios_base::seekdir dir = std::ios_base::beg) :
      filename(that.filename) {
      // TODO: there is probably a better way to copy the stream
      stream.open(filename.c_str(), std::ios::in | std::ios::binary);
      stream.seekg(off, dir);
      current_position = stream.tellg();
      populateChunk(stream, current_chunk);
      std::cout << "MessageChunksIterator Copy constructor" << std::endl;
    }

    // TODO: Prefix and Postfix operators currently do the same thing!
    // Prefix increment operator
    MessageChunksIterator& operator++() {
      std::cout << "Prefix operator++() " << std::endl;
      populateChunk(stream, current_chunk);
      return *this;
    };

    // Postfix increment operator
    MessageChunksIterator& operator++(int) {
      std::cout << "Postfix operator++(int) " << std::endl;
      populateChunk(stream, current_chunk);
      return *this;
    };

    bool operator==(const MessageChunksIterator<T>& rhs) const {
      // the filenames and the stream positions must match
      return filename == rhs.filename && current_position == rhs.current_position;
    };

    bool operator!=(const MessageChunksIterator<T>& rhs) const {
      // the filenames and the stream positions must match
      return filename != rhs.filename || current_position != rhs.current_position;
    };

    // return the parsed results for the current chunk
    const T& operator*() {
      return current_chunk;
    };

    T* const operator->() {
      return &current_chunk;
    };

    // TODO - remove the operator<< - for testing only
    friend std::ostream &operator<<(std::ostream &out, MessageChunksIterator& iterator) {
      out << "ostream &operator<< " << iterator.current_position;
      return out;
    };

    MessageChunksIterator begin() const {
      return MessageChunksIterator(*this, static_cast<std::streamoff>(0), std::ios_base::beg);
    };

    MessageChunksIterator end() const {
      return MessageChunksIterator(*this, static_cast<std::streamoff>(0), std::ios_base::end);
    };

    /*
    MessageChunksIterator& operator=(const MessageChunksIterator<T>& that)  {
      if (this != &that) {
        // TODO
      }
      return *this;
    };
    */
  };
}

#endif // GOBY_MESSAGE_CHUNKS_H
