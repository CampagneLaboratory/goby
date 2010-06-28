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
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include "common.h"
#include "Alignments.pb.h"

namespace goby {
  #define GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH 8      // length of the delimiter tag (in bytes)

  // details of an individual chunk of protocol buffer messages
  struct MessageChunkInfo {
    // the position within the stream
    std::streampos position;
    // the length of the chunk
    size_t length;
  };

  template <typename T> class MessageChunksIterator : public std::iterator<std::input_iterator_tag, T> {
    // the name of the chunked file
    std::string filename;

    // positions of compressed alignment collections within the goby chunked file
    std::vector<MessageChunkInfo> chunks;

    // index to the current chunk details
    int chunkIndex;

    // the current processed chunk
    T currentChunk;

    // populate T with the data from the given chunk
    T populateChunk(T& chunk, const MessageChunkInfo& chunkInfo) {
       std::ifstream stream(filename.c_str(), std::ios::in | std::ios::binary);

      // position the stream to the current chunk location
      stream.seekg(chunkInfo.position, std::ios::beg);

      // Set up a stream that will only read up the current chunk length
      const size_t compressedChunkLength = chunkInfo.length;
      google::protobuf::io::IstreamInputStream istream(&stream);
      google::protobuf::io::LimitingInputStream rawChunkStream(&istream, compressedChunkLength);

      // uncompress the buffer so that it can be parsed
      google::protobuf::io::GzipInputStream chunkStream(&rawChunkStream);

      // the stream may not get all read in at once so we may need to copy in stages
      void* uncompressedChunk = NULL;
      size_t chunkSize = 0;

      const void* buffer;
      int bufferSize;
      while (chunkStream.Next(&buffer, &bufferSize)) {
        // store the end location of the chunk buffer
        int index = chunkSize;

        // resize the chunk buffer to fit the new data just read
        chunkSize += bufferSize;
        uncompressedChunk = (void *)realloc(uncompressedChunk, chunkSize);

        // and append the new data over to the end of the header buffer
        ::memcpy(reinterpret_cast<char*>(uncompressedChunk) + index, buffer, bufferSize);
      }

      // populate the current object from the uncompressed data
      if (!chunk.ParseFromArray(uncompressedChunk, chunkSize)) {
        std::cerr << "Failed to parse message chunk from " << filename << std::endl;
      }

      // free up the temporary buffers and close the stream
      ::free(uncompressedChunk);

      stream.close();

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
    // TODO: MessageChunksIterator(int fileDescriptor);
    // TODO: MessageChunksIterator(istream* stream);
    MessageChunksIterator(const std::string& filename) : filename(filename) {
      chunks.clear();

      std::ifstream stream(filename.c_str(), std::ios::in | std::ios::binary);

      int chunkNumber = 0;
      // get the positions each of the chunks in the file
      while (stream.good()) {
        // each chunk is delimited by DELIMITER_LENGTH bytes
        stream.seekg(GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH, std::ios::cur);

        // then the size of the next chunk follows
        const int size = readInt(stream);
#if GOBY_DEBUG
        std::cout << "length of chunk #" << chunkNumber++ << " is " << size << std::endl;
#endif // GOBY_DEBUG

        // the last chunk has a size of zero bytes
        if (!stream.eof() && size != 0) {
          const std::streampos position = stream.tellg();
          MessageChunkInfo chunkInfo;
          chunkInfo.position = position;
          chunkInfo.length = size;
          chunks.push_back(chunkInfo);
          stream.seekg(size, std::ios::cur);
        } else {
          break;
        }
      }

      stream.close();

      this->chunkIndex = 0;
      this->currentChunk = T::default_instance();
    };

    virtual ~MessageChunksIterator(void) {
    };

    MessageChunksIterator(const MessageChunksIterator& reader)
      : filename(reader.filename),
        chunks(reader.chunks),
        chunkIndex(reader.chunkIndex),
        currentChunk(reader.currentChunk) {
      std::cout << "MessageChunksIterator Copy constructor" << std::endl;
    }

    MessageChunksIterator(const MessageChunksIterator& reader, int chunkIndex)
      : filename(reader.filename),
        chunks(reader.chunks),
        chunkIndex(chunkIndex),
        currentChunk(reader.currentChunk) {
      std::cout << "MessageChunksIterator Copy constructor" << std::endl;
    }

    // Prefix increment operator
    MessageChunksIterator& operator++() {
      std::cout << "Prefix operator++() " << std::endl;
      ++chunkIndex;
      return *this;
    };

    // Postfix increment operator
    MessageChunksIterator& operator++(int) {
      std::cout << "Postfix operator++(int) " << std::endl;
      chunkIndex++;
      return *this;
    };

    bool operator==(const MessageChunksIterator<T>& rhs) const {
      // the filenames must match and either the chunks are
      // both empty or the chunk positions are the same
      return filename == rhs.filename && chunkIndex == rhs.chunkIndex;
    };

    bool operator!=(const MessageChunksIterator<T>& rhs) const {
      // the filenames must match and either the chunks are
      // both empty or the chunk positions are the same
      return filename != rhs.filename || chunkIndex != rhs.chunkIndex;
    };

    // return the parsed results for the current chunk
    // the contract of the input iterator is that this is only done once
    // so we uncompress at this time rather than during the increment
    const T& operator*() {
      populateChunk(currentChunk, chunks[chunkIndex]);
      return currentChunk;
    };

    T* const operator->() {
      populateChunk(currentChunk, chunks[chunkIndex]);
      return &currentChunk;
    };

    // TODO - remove the operator<< - for testing only
    friend std::ostream &operator<<(std::ostream &out, const MessageChunksIterator& reader) {
      out << "ostream &operator<< " << reader.chunkIndex;
      return out;
    };

    MessageChunksIterator begin() const {
      return MessageChunksIterator(*this, 0);
    };

    MessageChunksIterator end() const {
      return MessageChunksIterator(*this, chunks.size());
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
