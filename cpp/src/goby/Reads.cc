//
// Copyright (C) 2009-2012 Institute for Computational Biomedicine,
//                         Weill Medical College of Cornell University
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

#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include "common.h"
#include "MessageChunks.h"
#include "Reads.h"
#include "Reads.pb.h"

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

using namespace std;

namespace goby {
  ReadEntryIterator::ReadEntryIterator(const int fd, std::streamoff startOffset=0, std::streamoff endOffset=0, ios_base::seekdir dir=ios_base::beg) :
    fd(fd),
    message_chunks_iterator(MessageChunksIterator<ReadCollection>(fd, startOffset, endOffset, dir)),
    message_chunks_iterator_end(MessageChunksIterator<ReadCollection>(fd, 0, 0, ios_base::end)),
    read_collection(new ReadCollection),
    current_read_index(0) {
  }

  ReadEntryIterator::ReadEntryIterator(const ReadEntryIterator& that) :
    fd(that.fd),
    message_chunks_iterator(that.message_chunks_iterator),
    message_chunks_iterator_end(that.message_chunks_iterator_end),
    read_collection(new ReadCollection),
    current_read_index(that.current_read_index) {
  }

  ReadEntryIterator::~ReadEntryIterator() {
    delete read_collection;
  }

  // Prefix increment operator
  ReadEntryIterator& ReadEntryIterator::operator++() {
    if (current_read_index != -1) {
      ++current_read_index;
      // if we're at the end of the current chunk, move on to the next
      if (current_read_index >= read_collection->reads_size()) {
        // if there is another chunk, get it otherwise set defaults
        if (message_chunks_iterator != message_chunks_iterator_end) {
          message_chunks_iterator++;
          current_read_index = 0;
        } else {
          read_collection->Clear();
          current_read_index = -1;
        }
      }
    } else {
        cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of fd " << fd << endl;
    }
    return *this;
  };

  // Postfix increment operator
  ReadEntryIterator& ReadEntryIterator::operator++(int) {
    if (current_read_index != -1) {
      current_read_index++;
      if (current_read_index >= read_collection->reads_size()) {
        // if there is another chunk, get it otherwise set defaults
        if (message_chunks_iterator != message_chunks_iterator_end) {
          message_chunks_iterator++;
          current_read_index = 0;
        } else {
          read_collection->Clear();
          current_read_index = -1;
        }
      }
    } else {
        cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of fd " << fd << endl;
    }
    return *this;
  };

  bool ReadEntryIterator::operator==(const ReadEntryIterator& rhs) const {
    // the filenames must match and the chunk/read indicies must be the same
    return current_read_index == rhs.current_read_index && message_chunks_iterator == rhs.message_chunks_iterator;
  };

  bool ReadEntryIterator::operator!=(const ReadEntryIterator& rhs) const {
    // if the filenames or the chunk/read indicies don't match, the reader is different.
    return current_read_index != rhs.current_read_index || message_chunks_iterator != rhs.message_chunks_iterator;
  };

  // return the parsed results for the current chunk
  const ReadEntry& ReadEntryIterator::operator*() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_read_index >= read_collection->reads_size() || current_read_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator_end) {
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

  const ReadEntry* const ReadEntryIterator::operator->() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_read_index >= read_collection->reads_size() || current_read_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator_end) {
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


  Reads::Reads(const string& basename) : filename(getBasename(basename) + ".compact-reads") {
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
    fd(::open(filename.c_str(), O_RDONLY | O_BINARY)),
    read_entry_iterator_end(new ReadEntryIterator(fd, 0, 0, ios_base::end)) {
    if (fd < 0) {
      cerr << "Error opening file: " << filename << endl;
    }
  }

  ReadsReader::~ReadsReader(void) {
    delete read_entry_iterator_end;
    if (fd >= 0) {
      ::close(fd);
    }
  }

  ReadEntryIterator ReadsReader::begin() const {
    return ReadEntryIterator(fd, 0, 0);
  }

  ReadEntryIterator ReadsReader::begin(std::streamoff startPosition=0, std::streamoff endPosition=0) const {
    return ReadEntryIterator(fd, startPosition, endPosition);
  };

  ReadEntryIterator* ReadsReader::beginPointer(std::streamoff startPosition=0, std::streamoff endPosition=0) const {
    return new ReadEntryIterator(fd, startPosition, endPosition);
  };

  ReadEntryIterator ReadsReader::end() const {
    return *read_entry_iterator_end;
  };

  const ReadEntryIterator* ReadsReader::endPointer() const {
    return read_entry_iterator_end;
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

  ReadsWriter::~ReadsWriter(void) {
    read_collection.Clear();
    delete message_chunks_writer;
  }

  void ReadsWriter::appendEntry() {
    // set fields in the new read entry
    ReadEntry *entry = read_collection.add_reads();
    entry->set_read_index(current_read_index++);
    if (sequence != NULL) {
      entry->set_sequence(sequence);
      entry->set_read_length(strlen(sequence));
      sequence = NULL;
    } else {
      entry->set_read_length(0);
    }

    if (description != NULL) {
      if (strlen(description) > 0) {
        entry->set_description(description);
      }
      description = NULL;
    }

    if (identifier != NULL) {
      if (strlen(identifier) > 0) {
        entry->set_read_identifier(identifier);
      }
      identifier = NULL;
    }

    if (quality_scores != NULL) {
      if (strlen(quality_scores) > 0) {
        entry->set_quality_scores(quality_scores);
      }
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
