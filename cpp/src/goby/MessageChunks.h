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

#pragma once

#ifndef GOBY_MESSAGE_CHUNKS_H
#define GOBY_MESSAGE_CHUNKS_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
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

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif  // _MSC_VER

namespace goby {
  #define GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH 8           // length of the delimiter tag (in bytes)
  #define GOBY_MESSAGE_CHUNK_DELIMITER_CONTENT 0xFF       // value of each byte in the delimeter
  #define GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK 10000  // default number of entries per chunk

  template <typename T> class MessageChunksIterator : public std::iterator<std::input_iterator_tag, T> {
    // the name of the chunked file
    std::string filename;

    // the file descriptor for the chunked file
    int fd;

    // Start offset
    std::streamoff startOffset;
    // End offset
    std::streamoff endOffset;

    // whether or not to close the file descriptor when the object is deleted
    bool close_fd_on_delete;

    // The underlying stream
    google::protobuf::io::ZeroCopyInputStream *input_stream;

    // whether or not this instance created the stream and therefore should delete it
    bool owns_input_stream;

    // the initial position in the file
    std::streampos initial_position;

    // The current chunk in the file
    int current_chunk_index;

    // the current processed chunk
    T *current_chunk;

    // the length of the current raw chunk
    size_t current_chunk_length;

    // True if the message chunks delimiter has already been read
    bool delimiterAlreadySkipped;

    /**
     * Move to startOffset and then to the start of the next chunk.
     */
    void initializePosition() {
      if (startOffset == 0) {
        return;
      }
      std::cerr << "skipping " << startOffset << " bytes in input file" << std::endl;
      const bool result = input_stream->Skip(startOffset);
      if (!result) {
        std::cerr << "failed to skip bytes." << std::endl;
        current_chunk_index = -1;
        current_chunk_length = 0;
        return;
      }
      // advance to the next chunk
      int b;
      int contiguousZeroBytes = 0;
      int skipped = 0;
      std::streamoff position = 0;

      // search though the input stream until a delimiter chunk or end of stream is reached
      while (true) {
        if ((endOffset != 0) && (position >= endOffset)) {
            break;
        }
        b = readByte(input_stream);
        if (b < 0) {
          break; // end of file
        }

        if (b == GOBY_MESSAGE_CHUNK_DELIMITER_CONTENT) {
          contiguousZeroBytes++;
        } else {
          contiguousZeroBytes = 0;
        }
        ++skipped;
        if (contiguousZeroBytes == GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH) {
          delimiterAlreadySkipped = true;
          break;
        }
        position = startOffset + skipped;
      }
      std::cerr << "New chunk at position " << input_stream->ByteCount() << std::endl;
    }

    // Move the stream poitner to the next chunk boundary or eof
    void advanceToNextChunk() {
      if (current_chunk_index != -1) {
        if (endOffset != 0 && input_stream->ByteCount() >= endOffset) {
            current_chunk_index = -1;
            current_chunk_length = 0;
        } else {
          // each chunk is delimited by DELIMITER_LENGTH bytes
          bool result;
          if (delimiterAlreadySkipped) {
              delimiterAlreadySkipped = false;
              result = true;
          } else {
              result = input_stream->Skip(GOBY_MESSAGE_CHUNK_DELIMITER_LENGTH);
          }
          if (!result) {
            current_chunk_index = -1;
            current_chunk_length = 0;
          } else {
            current_chunk_index++;

            // Set up a stream that will only read up the current chunk length
            current_chunk_length = readInt(input_stream);

#ifdef GOBY_DEBUG
            std::cout << "Chunk length: " << current_chunk_length << std::endl;
#endif

            // the last chunk has a length of zero bytes
            if (current_chunk_length == 0) {
              current_chunk_index = -1;
              current_chunk_length = 0;
            }
          }
        }
      } else {
        current_chunk_length = 0;
        std::cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of " << filename << std::endl;
      }
    }

    // populate T with the data from the given chunk
    // The assumption here is that the stream is positioned at the start of a chunk boundary
    void populateChunk(T *chunk) const {
      // if the stream is not valid just clear out the chunk
      if (current_chunk_index == -1) {
        chunk->Clear();
      } else {
        // explicitly set the block size and limit it to the chunk length to prevent "over read"
        google::protobuf::io::LimitingInputStream limiting_stream(input_stream, current_chunk_length);

        // and handle the fact that each chunk is compressed with gzip
        google::protobuf::io::GzipInputStream gzip_stream(&limiting_stream);

        // populate the current object from the compressed data
        if (!chunk->ParseFromZeroCopyStream(&gzip_stream)) {
          std::cerr << __FILE__ ":" << __LINE__ << " ParseFromZeroCopyStream() returned false " << filename << "(" << fd << ")" << std::endl;
        }
      }
    };

    // Java DataInput.readInt()
    static int readInt(google::protobuf::io::ZeroCopyInputStream *stream) {
      int bytes_needed = 4;           // the number of bytes needed to read an int value
      unsigned char ch[4];            // raw bytes (4) that represent the int value
      unsigned char *cp = &ch[0];

      const void *buffer;             // buffer from the input stream
      int size;                       // size of the buffer returned from the input stream

      // need to loop in case the stream does not return the number of bytes needed in a single pass
      while (bytes_needed > 0) {
        if (stream->Next(&buffer, &size)) {
          // the stream may return a size less than we need (including zero) and still have more data
          if (size > 0) {
            // only use what is needed (which may be less than what we got)
            const int bytes_used = std::min<int>(bytes_needed, size);

            // copy the data from the buffer to the local copy
            ::memcpy(cp, buffer, bytes_used);

            // reset the stream and local buffer pointers to just after the data read
            stream->BackUp(size - bytes_used);
            cp += bytes_used;

            // and update the number of bytes still needed
            bytes_needed -= bytes_used;
          }
        } else {
          // reached eof or there was an issue reading from the stream
          std::cerr << __FILE__ ":" << __LINE__ << " ZeroCopyInputStream::Next() returned false" << std::endl;
          return 0;
        }
      };

      return (ch[0] << 24) + (ch[1] << 16) + (ch[2] << 8) + (ch[3] << 0);
    };

    // Reads an unsigned byte from the file
    // Returns -1 at end of file (or error)
    static int readByte(google::protobuf::io::ZeroCopyInputStream *stream) {
      int bytes_needed = 1;           // the number of bytes needed to read an int value
      unsigned char ch[1];            // data read from the stream
      unsigned char *cp = &ch[0];

      const void *buffer;             // buffer from the input stream
      int size;                       // size of the buffer returned from the input stream

      // need to loop in case the stream does not return the number of bytes needed in a single pass
      while (bytes_needed > 0) {
        if (stream->Next(&buffer, &size)) {
          // the stream may return a size less than we need (including zero) and still have more data
          if (size > 0) {
            // only use what is needed (which may be less than what we got)
            const int bytes_used = std::min<int>(bytes_needed, size);

            // reset the stream and local buffer pointers to just after the data read
            stream->BackUp(size - bytes_used);

            // copy the data from the buffer to the local copy
            ::memcpy(cp, buffer, bytes_used);

            bytes_needed -= bytes_used;
          }
        } else {
          // reached eof or there was an issue reading from the stream
          return -1;
        }
      }

      return ch[0];
    }

  public:
    // TOOD: currently assumes that the position is at the start of chunk boundary
    MessageChunksIterator<T>(const int fd, std::streamoff startOffsetVal=0, std::streamoff endOffsetVal=0, std::ios_base::seekdir dir = std::ios_base::beg) :
        fd(fd),
        startOffset(startOffsetVal),
        endOffset(endOffsetVal),
        filename(""),
        current_chunk(new T),
        current_chunk_index(0),
        current_chunk_length(0),
        close_fd_on_delete(false),
        delimiterAlreadySkipped(false),
        owns_input_stream(true) {

      if ((endOffset != 0) && (endOffset < startOffset)) {
        // Swap start and end
        std::streamoff tempOffset = startOffset;
        startOffset = endOffset;
        endOffset = tempOffset;
      }

      current_chunk->Clear();
      if ((dir == std::ios_base::beg && startOffset < 0) || (dir == std::ios_base::end && startOffset >= 0)) {
        input_stream = NULL;
        current_chunk_index = -1;
      } else {
        input_stream = new google::protobuf::io::FileInputStream(fd);
        if (startOffset != 0) {
          initializePosition();
        }
        advanceToNextChunk();
      }
    }

    MessageChunksIterator<T>(const std::string& filename, std::streamoff startOffsetVal=0, std::streamoff endOffsetVal=0, std::ios_base::seekdir dir = std::ios_base::beg) :
        filename(filename),
        startOffset(startOffsetVal),
        endOffset(endOffsetVal),
        current_chunk(new T),
        current_chunk_index(0),
        current_chunk_length(0),
        close_fd_on_delete(true),
        delimiterAlreadySkipped(false),
        owns_input_stream(true) {

      if (endOffset < startOffset) {
        // Swap start and end
        std::streamoff tempOffset = startOffset;
        startOffset = endOffset;
        endOffset = tempOffset;
      }

      current_chunk->Clear();
      if ((dir == std::ios_base::beg && startOffset < 0) || (dir == std::ios_base::end && startOffset >= 0)) {
        input_stream = NULL;
        current_chunk_index = -1;
      } else {
        fd = ::open(filename.c_str(), O_RDONLY | O_BINARY);
        if (fd < 0) {
          std::cerr << "Error opening file: " << filename << std::endl;
          input_stream = NULL;
          current_chunk_index = -1;
        } else {
          input_stream = new google::protobuf::io::FileInputStream(fd);
          if (startOffset != 0) {
            initializePosition();
          }
          advanceToNextChunk();
        }
      }
    }

    MessageChunksIterator(const MessageChunksIterator<T>& that) :
        filename(that.filename),
        startOffset(that.startOffset),
        endOffset(that.endOffset),
        current_chunk(new T),
        current_chunk_index(that.current_chunk_index),
        current_chunk_length(that.current_chunk_length),
        fd(that.fd), close_fd_on_delete(false),
        input_stream(that.input_stream), owns_input_stream(false) {
      current_chunk->Clear();
    }

/*
REMOVE?
    MessageChunksIterator(const MessageChunksIterator<T>& that, std::streamoff off, std::ios_base::seekdir dir = std::ios_base::beg) :
      filename(that.filename),
      current_chunk(new T),
      current_chunk_index(-1),
      current_chunk_length(0),
      input_stream(NULL),
      fd(0),
      close_fd_on_delete(false) {
        // TODO
        current_chunk->Clear();
    }
*/

    virtual ~MessageChunksIterator(void) {
      delete current_chunk;

      if (owns_input_stream) {
        delete input_stream;
      }
      if (close_fd_on_delete) {
        ::close(fd);
      }
    };

    // TODO: Prefix and Postfix operators currently do the same thing!
    // Prefix increment operator
    MessageChunksIterator& operator++() {
      advanceToNextChunk();
      return *this;
    };

    // Postfix increment operator
    MessageChunksIterator& operator++(int) {
      advanceToNextChunk();
      return *this;
    };

    bool operator==(const MessageChunksIterator<T>& rhs) const {
      // the filenames and the stream positions must match
      return fd == rhs.fd && filename == rhs.filename && current_chunk_index == rhs.current_chunk_index;
    };

    bool operator!=(const MessageChunksIterator<T>& rhs) const {
      // the filenames and the stream positions must match
      return fd != rhs.fd || filename != rhs.filename || current_chunk_index != rhs.current_chunk_index;
    };

    // return the parsed results for the current chunk
    const T& operator*() const {
      populateChunk(current_chunk);
      return *current_chunk;
    };

    const T* const operator->() const {
      populateChunk(current_chunk);
      return current_chunk;
    };
  };

  template <typename T> class MessageChunksReader {
    // the name of the chunked file
    std::string filename;

    // the file descriptor for the chunked file
    int fd;

    // whether or not to close the file descriptor when the object is deleted
    bool close_fd_on_delete;

  public:
    MessageChunksReader(const std::string& file_name) :
      filename(file_name),
      close_fd_on_delete(true) {
      fd = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0644);
      if (fd < 0) {
        std::cerr << "Error opening file: " << filename << std::endl;
      }
    }

    MessageChunksReader(const int fd) :
      filename(""),
      fd(fd),
      close_fd_on_delete(false) {
    }

    ~MessageChunksReader() {
      if (close_fd_on_delete) {
        ::close(fd);
      }
    }

    MessageChunksIterator<T> begin() const {
      return MessageChunksIterator<T>(fd);
    };

    MessageChunksIterator<T> end() const {
      return MessageChunksIterator<T>(fd, 0, std::ios_base::end);
    };

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
      unsigned char tmp[4];
      tmp[0] = (v >> 24) & 0xFF;
      tmp[1] = (v >> 16) & 0xFF;
      tmp[2] = (v >> 8) & 0xFF;
      tmp[3] = v & 0xFF;
      stream->WriteRaw(tmp, 4);
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

      if (total_entries_written == 0 || number_appended > 0) {
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
