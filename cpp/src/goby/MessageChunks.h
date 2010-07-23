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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include "common.h"
#include "Alignments.pb.h"

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

namespace goby {
  #define GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH 8           // length of the delimiter tag (in bytes)
  #define GOBY_MESSAGE_CHUNK_DELIMITER_CONTENT 0xFF       // value of each byte in the delimeter
  #define GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK 10000  // default number of entries per chunk

  template <typename T> class MessageChunksIterator : public std::iterator<std::input_iterator_tag, T> {
    // the name of the chunked file
    std::string filename;

    // The underlying stream
    std::ifstream stream;

    // The current position in the stream (we keep this since "tellg()" is not const)
    std::streampos current_position;

    // the current processed chunk
    T *current_chunk;

    // the length of the current raw chunk
    size_t current_chunk_length;

    // Move the stream poitner to the next chunk boundary or eof
    void advanceToNextChunk(std::ifstream& stream) {
      if (current_position != static_cast<std::streampos>(-1) && stream.good()) {
        // each chunk is delimited by DELIMITER_LENGTH bytes
        stream.seekg(GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH, std::ios::cur);
        current_position = stream.tellg();

        // Set up a stream that will only read up the current chunk length
        current_chunk_length = readInt(stream);
        current_position = stream.tellg();

#ifdef GOBY_DEBUG
        //std::cout << "Chunk length: " << current_chunk_length << std::endl;
#endif

        // the last chunk has a length of zero bytes
        if (stream.eof() || current_chunk_length == 0) {
          current_position = static_cast<std::streampos>(-1);
          current_chunk_length = 0;
        }
      } else {
        current_position = static_cast<std::streampos>(-1);
        current_chunk_length = 0;
      }
    }

    // populate T with the data from the given chunk
    // The assumption here is that the stream is positioned at the start of a chunk boundary
    T* populateChunk(std::ifstream& stream, T *chunk) {
      // if the stream is not valid just return an empty chunk
      if (current_position == static_cast<std::streampos>(-1) || !stream.good()) {
        current_position = static_cast<std::streampos>(-1);
        chunk->Clear();
        return chunk;
      }

      const std::streampos chunkStartPosition = current_position;
      const std::streampos nextChunkStartPosition = current_position + static_cast<std::streamoff>(current_chunk_length);

      // explicitly set the block size and limit it to the chunk length to prevent "over read"
      google::protobuf::io::IstreamInputStream istream(&stream, current_chunk_length);
      google::protobuf::io::LimitingInputStream rawChunkStream(&istream, current_chunk_length);

      // and handle the fact that each chunk is compressed with gzip
      google::protobuf::io::GzipInputStream gzipChunkStream(&rawChunkStream);

      // since the chunks may be large, we need to increase the limit
      // TODO: we may want to be smarter about the actual limit here
      google::protobuf::io::CodedInputStream codedStream(&gzipChunkStream);
      codedStream.SetTotalBytesLimit(INT_MAX, -1);

      // populate the current object from the compressed data
      if (!chunk->ParseFromCodedStream(&codedStream)) {
        std::cerr << "Failed to parse message chunk from " << filename << std::endl;
      }

      current_position = stream.tellg();

      // may need to adjust the position in the stream since the parsing above may have "over-read"
      // TODO - we can probably remove this check since the block size is set above
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
        current_chunk(new T),
        current_chunk_length(0) {
      stream.open(filename.c_str(), std::ios::in | std::ios::binary);
      stream.seekg(position);

      if (!stream.good()) {
        std::cerr << "Failed to open " << filename << std::endl;
      }

      advanceToNextChunk(stream);
      current_chunk->Clear();
    };

    MessageChunksIterator(const MessageChunksIterator<T>& that) :
        filename(that.filename),
        current_position(that.current_position),
        current_chunk(new T),
        current_chunk_length(that.current_chunk_length) {
      // TODO: there is probably a better way to copy the stream
      stream.open(filename.c_str(), std::ios::in | std::ios::binary);
      stream.seekg(current_position);

      current_chunk->Clear();
    }

    MessageChunksIterator(const MessageChunksIterator<T>& that, std::streamoff off, std::ios_base::seekdir dir = std::ios_base::beg) :
      filename(that.filename),
      current_chunk(new T) {
      // TODO: there is probably a better way to copy the stream
      stream.open(filename.c_str(), std::ios::in | std::ios::binary);
      stream.seekg(off, dir);
      current_position = stream.tellg();
      advanceToNextChunk(stream);
    }

    virtual ~MessageChunksIterator(void) {
      stream.close();
      delete current_chunk;
    };

    // TODO: Prefix and Postfix operators currently do the same thing!
    // Prefix increment operator
    MessageChunksIterator& operator++() {
      std::cout << "Prefix operator++() " << std::endl;
      advanceToNextChunk(stream);
      return *this;
    };

    // Postfix increment operator
    MessageChunksIterator& operator++(int) {
      advanceToNextChunk(stream);
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
      populateChunk(stream, current_chunk);
      return *current_chunk;
    };

    const T* const operator->() {
      populateChunk(stream, current_chunk);
      return current_chunk;
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

  template <typename T> class MessageChunksWriter {
    // the name of the chunked file
    std::string filename;

    // the file descriptor for the chunked file
    int fd;

    // The underlying streams
    google::protobuf::io::FileOutputStream *file_stream;
    google::protobuf::io::CodedOutputStream *coded_stream;

    // number of entries to store before serializing to disk
    unsigned number_of_entries_per_chunk;

    // The number of messages appended in a chunk.
    unsigned number_appended;

    // The total number of logical entries written to the output.
    unsigned long total_entries_written;

    // a temporary buffer for processing (class member to avoid creating eveytime)
    std::string buffer;

    // Java DataOutput.writeInt()
    static void writeInt(google::protobuf::io::CodedOutputStream *stream, int v) {
      unsigned char *tmp = new unsigned char[4];
      tmp[0] = (v >> 24) & 0xFF;
      tmp[1] = (v >> 16) & 0xFF;
      tmp[2] = (v >> 8) & 0xFF;
      tmp[3] = v & 0xFF;
      stream->WriteRaw(tmp, 4);
      delete tmp;
    }

    // content for delimiter between each chunk
    static std::string delimiter;

    // indicates the writer has been explicitly closed by the user
    bool closed;

    // write termination record to the file and close all the underlying streams
    void close() {
      // Write the final delimiter
      coded_stream->WriteString(delimiter);

      // and write zero for the final chunk size
      writeInt(coded_stream, 0);

      // close all the streams
      delete coded_stream;
      delete file_stream;
      ::close(fd);

      closed = true;
    }

  public:
    MessageChunksWriter(const std::string& filename, unsigned number_of_entries_per_chunk = GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK) :
        filename(filename),
        number_of_entries_per_chunk(number_of_entries_per_chunk),
        number_appended(0),
        total_entries_written(0) {
      fd = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0644);
      if (fd < 0) {
        std::cerr << "Error opening file: " << filename << std::endl;
        file_stream = NULL;
        coded_stream = NULL;
        closed = true;
      } else {
        file_stream = new google::protobuf::io::FileOutputStream(fd);
        coded_stream = new google::protobuf::io::CodedOutputStream(file_stream);
        closed = false;
      }
    }

    virtual ~MessageChunksWriter(void) {
      if (!closed) {
        close();
      }
    }

    // Write the collection as needed to the output stream. When the number of entries
    // per chunk is reached, the chunk is written to disk and the collection cleared. Clients
    // can just keep adding to the collection and call writeAsNeeded for every entry.
    void writeAsNeeded(T* const collection, int multiplicity = 1) {
        total_entries_written += multiplicity;
        if (++number_appended >= number_of_entries_per_chunk) {
            flush(collection);
        }
        // TODO: return currentChunkStartOffset
    }

    // force writing the collection to the output stream.
    void flush(T* const collection) {
      // Write the delimiter between two chunks
      coded_stream->WriteString(delimiter);

      google::protobuf::io::StringOutputStream bufferStream(&buffer);
      google::protobuf::io::GzipOutputStream compressedStream(&bufferStream);
      //std::cout << collection->DebugString() << std::endl;

      if (!collection->SerializeToZeroCopyStream(&compressedStream)) {
        std::cerr << "There was a problem compressing the collection" << std::endl;
      }
      compressedStream.Close();

      // write the length of the compressed chunk to the file
      const int bufferSize = buffer.size();
      writeInt(coded_stream, bufferSize);

      // then write the actual compressed chunk
      coded_stream->WriteString(buffer);

      number_appended = 0;
      buffer.clear();
      collection->Clear();
      // TODO: delete elements of the collection?

    /*
       totalBytesWritten += serializedSize + 4 + DELIMITER_LENGTH;
    */
    }

    // write any remaining items in the collection to disk and close up the file
    // no more data items can be written after this method is called
    void close(T* const collection) {
      // Write any remaining items in the collection
      flush(collection);
      close();
    }
  };

  // initialize string for delimter between chunks
  template <typename T> std::string MessageChunksWriter<T>::delimiter(GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH, GOBY_MESSAGE_CHUNK_DELIMITER_CONTENT);

}

#ifdef _MSC_VER
#pragma warning(pop)  // Restores the warning state.
#endif

#endif // GOBY_MESSAGE_CHUNKS_H
