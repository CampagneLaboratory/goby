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

#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include "common.h"
#include "MessageChunks.h"
#include "Reads.h"

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

using namespace std;

namespace goby {
  ReadsIterator::ReadsIterator(const int fd) :
    fd(fd),
    filename(""),
    close_on_delete(false),
    message_chunks_iterator(MessageChunksIterator<ReadCollection>(fd)),
    read_collection(new ReadCollection),
    current_read_index(0) {
  }

  ReadsIterator::ReadsIterator(const int fd, const std::streampos position, std::ios_base::seekdir dir = std::ios_base::beg) :
    fd(fd),
    filename(""),
    close_on_delete(false),
    message_chunks_iterator(MessageChunksIterator<ReadCollection>(fd, position, dir)),
    read_collection(new ReadCollection),
    current_read_index(-1) {
  }

  ReadsIterator::ReadsIterator(const string& filename, const std::streampos position = 0, std::ios_base::seekdir dir = std::ios_base::beg) :
    filename(filename),
    fd(::open(filename.c_str(), O_RDONLY | O_BINARY)),
    message_chunks_iterator(MessageChunksIterator<ReadCollection>(fd, position, dir)),
    close_on_delete(true),
    read_collection(new ReadCollection),
    current_read_index(0) {

    if (fd < 0) {
      std::cerr << "Error opening file: " << filename << std::endl;
    }
  }

  ReadsIterator::ReadsIterator(const ReadsIterator& that) :
    filename(that.filename),
    fd(that.fd),
    close_on_delete(false),
    message_chunks_iterator(that.message_chunks_iterator),
    read_collection(new ReadCollection),
    current_read_index(that.current_read_index) {
  }

  ReadsIterator::~ReadsIterator() {
    delete read_collection;
    if (close_on_delete) {
      ::close(fd);
    }
  }

  // Prefix increment operator
  ReadsIterator& ReadsIterator::operator++() {
    ++current_read_index;
    // if we're at the end of the current chunk, move on to the next
    if (current_read_index >= read_collection->reads_size()) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator.end()) {
        message_chunks_iterator++;
        current_read_index = 0;
      } else {
        read_collection->Clear();
        current_read_index = -1;
      }
    }
    return *this;
  };

  // Postfix increment operator
  ReadsIterator& ReadsIterator::operator++(int) {
    current_read_index++;
    if (current_read_index >= read_collection->reads_size()) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator.end()) {
        message_chunks_iterator++;
        current_read_index = 0;
      } else {
        read_collection->Clear();
        current_read_index = -1;
      }
    }
    return *this;
  };

  bool ReadsIterator::operator==(const ReadsIterator& rhs) const {
    // the filenames must match and the chunk/read indicies must be the same
    return filename == rhs.filename && current_read_index == rhs.current_read_index && message_chunks_iterator == rhs.message_chunks_iterator;
  };

  bool ReadsIterator::operator!=(const ReadsIterator& rhs) const {
    // if the filenames or the chunk/read indicies don't match, the reader is different.
    return filename != rhs.filename || current_read_index != rhs.current_read_index || message_chunks_iterator != rhs.message_chunks_iterator;
  };

  // return the parsed results for the current chunk
  const ReadEntry& ReadsIterator::operator*() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_read_index >= read_collection->reads_size() || current_read_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator.end()) {
        *read_collection = *message_chunks_iterator;
      } else {
        read_collection->Clear();
      }
    }

    if (read_collection->reads_size() > current_read_index) {
      return read_collection->reads().Get(current_read_index);
    } else {
      return ReadEntry::default_instance();
    }
  };

  const ReadEntry* const ReadsIterator::operator->() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_read_index >= read_collection->reads_size() || current_read_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator.end()) {
        *read_collection = *message_chunks_iterator;
      } else {
        read_collection->Clear();
      }
    }

    if (read_collection->reads_size() > current_read_index) {
      return &read_collection->reads().Get(current_read_index);
    } else {
      return &ReadEntry::default_instance();
    }
  };


  Reads::Reads(const std::string& basename) : filename(getBasename(basename) + ".compact-reads") {
  }

  Reads::~Reads(void) {
  }

  string Reads::getBasename(const char* filename) {
    return getBasename(string(filename));
  }

  string Reads::getBasename(const string& filename) {
    const string COMPACT_READS_FILE_EXTS[] = {
      ".compact-reads"
    };

    const size_t len = sizeof(COMPACT_READS_FILE_EXTS) / sizeof(COMPACT_READS_FILE_EXTS[0]);
    const size_t dotindex = filename.find_last_of(".");
    if (dotindex != string::npos) {
      const string extension = filename.substr(dotindex);
      for (int i = 0; i < len; i++) {
        if (extension.compare(COMPACT_READS_FILE_EXTS[i]) == 0) {
          return filename.substr(0, dotindex);
        }
      }
    }

    return filename;
  }

  ReadsReader::ReadsReader(const string& filename) : Reads(filename),
    fd(::open(filename.c_str(), O_RDONLY | O_BINARY)) {
    if (fd < 0) {
      std::cerr << "Error opening file: " << filename << std::endl;
    }
  }

  ReadsReader::~ReadsReader(void) {
    if (fd >= 0) {
      ::close(fd);
    }
  }

  ReadsIterator ReadsReader::begin() const {
    return ReadsIterator(fd);
  };

  ReadsIterator ReadsReader::end() const {
    return ReadsIterator(fd, static_cast<std::streamoff>(0), std::ios_base::end);
  };

  ReadsWriter::ReadsWriter(const string& filename, unsigned number_of_entries_per_chunk) : Reads(getBasename(filename)),
    message_chunks_writer(new MessageChunksWriter<ReadCollection>(filename, number_of_entries_per_chunk)),
    read_collection(ReadCollection::default_instance()),
    current_read_index(0),
    sequence(NULL),
    description(NULL),
    identifier(NULL),
    quality_scores(NULL) {
  }

  ReadsWriter::ReadsWriter(const Reads& reads) : Reads(reads) {
    // TODO: testing only
    this->filename = "foo.compact-reads";
  }

  ReadsWriter::~ReadsWriter(void) {
    read_collection.Clear();
    delete message_chunks_writer;
  }

  void ReadsWriter::appendEntry() {
    // set fields in the new read entry
    goby::ReadEntry *entry = read_collection.add_reads();
    entry->set_read_index(current_read_index++);
    if (sequence != NULL) {
      entry->set_sequence(sequence);
      entry->set_read_length(strlen(sequence));
      sequence = NULL;
    } else {
      entry->set_read_length(0);
    }

    if (description != NULL) {
      entry->set_description(description);
      description = NULL;
    }

    if (identifier != NULL) {
      entry->set_read_identifier(identifier);
      identifier = NULL;
    }

    if (quality_scores != NULL) {
      entry->set_quality_scores(quality_scores);
      quality_scores = NULL;
    }

    //cout << entry->DebugString() << endl;

    // and pass it along to the chunk writer
    message_chunks_writer->writeAsNeeded(&read_collection);
  }

  // flush any remaining items to the file and close the underlying streams
  void ReadsWriter::close() {
    message_chunks_writer->close(&read_collection);
  }
}

#ifdef _MSC_VER
#pragma warning(pop)  // Restores the warning state.
#endif
