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
#include "Alignments.h"

#ifdef _MSC_VER
// Disable Microsoft deprecation warnings for POSIX functions called from this class (open, close)
#pragma warning(push)
#pragma warning(disable:4996)
#endif

using namespace std;

namespace goby {
  AlignmentEntryIterator::AlignmentEntryIterator(const int fd, streamoff off = 0, ios_base::seekdir dir = ios_base::beg) :
    fd(fd),
    message_chunks_iterator(MessageChunksIterator<AlignmentCollection>(fd, off, dir)),
    message_chunks_iterator_end(MessageChunksIterator<AlignmentCollection>(fd, 0, ::ios_base::end)),
    alignment_collection(new AlignmentCollection),
    current_alignment_entry_index(0) {
  }

  AlignmentEntryIterator::AlignmentEntryIterator(const AlignmentEntryIterator& that) :
    fd(that.fd),
    message_chunks_iterator(that.message_chunks_iterator),
    message_chunks_iterator_end(that.message_chunks_iterator_end),
    alignment_collection(new AlignmentCollection),
    current_alignment_entry_index(that.current_alignment_entry_index) {
  }

  AlignmentEntryIterator::~AlignmentEntryIterator() {
    delete alignment_collection;
  }

  // Prefix increment operator
  AlignmentEntryIterator& AlignmentEntryIterator::operator++() {
    if (current_alignment_entry_index != -1) {
      ++current_alignment_entry_index;
      // if we're at the end of the current chunk, move on to the next
      if (current_alignment_entry_index >= alignment_collection->alignment_entries_size()) {
        // if there is another chunk, get it otherwise set defaults
        if (message_chunks_iterator != message_chunks_iterator_end) {
          message_chunks_iterator++;
          current_alignment_entry_index = 0;
        } else {
          alignment_collection->Clear();
          current_alignment_entry_index = -1;
        }
      }
    } else {
        std::cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of fd " << fd << std::endl;
    }
    return *this;
  };

  // Postfix increment operator
  AlignmentEntryIterator& AlignmentEntryIterator::operator++(int) {
    if (current_alignment_entry_index != -1) {
      current_alignment_entry_index++;
      if (current_alignment_entry_index >= alignment_collection->alignment_entries_size()) {
        // if there is another chunk, get it otherwise set defaults
        if (message_chunks_iterator != message_chunks_iterator_end) {
          message_chunks_iterator++;
          current_alignment_entry_index = 0;
        } else {
          alignment_collection->Clear();
          current_alignment_entry_index = -1;
        }
      }
    } else {
        std::cerr << __FILE__ ":" << __LINE__ << " - Attempt to advance past end of fd " << fd << std::endl;
    }
    return *this;
  };

  bool AlignmentEntryIterator::operator==(const AlignmentEntryIterator& rhs) const {
    // the filenames must match and the chunk/read indicies must be the same
    return current_alignment_entry_index == rhs.current_alignment_entry_index && message_chunks_iterator == rhs.message_chunks_iterator;
  };

  bool AlignmentEntryIterator::operator!=(const AlignmentEntryIterator& rhs) const {
    // if the filenames or the chunk/read indicies don't match, the reader is different.
    return current_alignment_entry_index != rhs.current_alignment_entry_index || message_chunks_iterator != rhs.message_chunks_iterator;
  };

  // return the parsed results for the current chunk
  const AlignmentEntry& AlignmentEntryIterator::operator*() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_alignment_entry_index >= alignment_collection->alignment_entries_size() || current_alignment_entry_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator_end) {
        *alignment_collection = *message_chunks_iterator;
      } else {
        alignment_collection->Clear();
      }
    }

    if (alignment_collection->alignment_entries_size() > current_alignment_entry_index) {
      return alignment_collection->alignment_entries().Get(current_alignment_entry_index);
    } else {
      return AlignmentEntry::default_instance();
    }
  };

  const AlignmentEntry* const AlignmentEntryIterator::operator->() {
    // if we're at the end of the current chunk or at the beginning of a new one
    if (current_alignment_entry_index >= alignment_collection->alignment_entries_size() || current_alignment_entry_index == 0) {
      // if there is another chunk, get it otherwise set defaults
      if (message_chunks_iterator != message_chunks_iterator_end) {
        *alignment_collection = *message_chunks_iterator;
      } else {
        alignment_collection->Clear();
      }
    }

    if (alignment_collection->alignment_entries_size() > current_alignment_entry_index) {
      return &alignment_collection->alignment_entries().Get(current_alignment_entry_index);
    } else {
      return &AlignmentEntry::default_instance();
    }
  };

  Alignment::Alignment(const string& basename) {
    this->basename = basename;
    this->header = AlignmentHeader::default_instance();
  }

  Alignment::~Alignment(void) {
  }

  string Alignment::getBasename(const char* filename) {
    return getBasename(string(filename));
  }

  string Alignment::getBasename(const string& filename) {
    const string COMPACT_ALIGNMENT_FILE_EXTS[] = {
      ".entries", ".header", ".tmh", ".stats", ".counts", ".index"
    };

    const size_t len = sizeof(COMPACT_ALIGNMENT_FILE_EXTS) / sizeof(COMPACT_ALIGNMENT_FILE_EXTS[0]);
    const size_t dotindex = filename.find_last_of(".");
    if (dotindex != string::npos) {
      const string extension = filename.substr(dotindex);
      for (int i = 0; i < len; i++) {
        if (extension.compare(COMPACT_ALIGNMENT_FILE_EXTS[i]) == 0) {
          return filename.substr(0, dotindex);
        }
      }
    }

    return filename;
  }

  vector<unsigned> Alignment::getQueryLengths() const {
    vector<unsigned> query_lengths(header.number_of_queries());
    if (hasConstantQueryLength()) {
      // The alignment has constant query lengths
      query_lengths.assign(header.number_of_queries(), getConstantQuerylength());
    } else {
      query_lengths.assign(header.query_length().begin(), header.query_length().end());
    }
  
    return query_lengths;
  }

  vector<unsigned> Alignment::getTargetLengths() const {
    vector<unsigned> target_lengths(header.target_length().begin(), header.target_length().end());
    return target_lengths;
  }
  
  AlignmentReader::AlignmentReader(const string& basename) : Alignment(basename), alignment_entry_iterator_end(NULL) {
    // open the "header" file
    const string header_filename = basename + ".header";
    const int fd = ::open(header_filename.c_str(), O_RDONLY | O_BINARY);

    if (fd > 0) {
      // uncompress file into memory so that it can be parsed
      google::protobuf::io::FileInputStream header_file_stream(fd);
      google::protobuf::io::GzipInputStream gzipHeaderStream(&header_file_stream);

      // populate the alignment header object from the uncompressed data
      if (!header.ParseFromZeroCopyStream(&gzipHeaderStream)) {
        cerr << "Failed to parse alignment header file: " << header_filename << endl;
      }

      // close the streams and files
      header_file_stream.Close();
    } else {
      cerr << "Failed to open alignment header file: " << header_filename << endl;
    }

    // populate the target identifiers
    google::protobuf::RepeatedPtrField<const goby::IdentifierInfo>::const_iterator target_mapping_iterator;
    for (target_mapping_iterator = header.target_name_mapping().mappings().begin(); target_mapping_iterator != header.target_name_mapping().mappings().end(); target_mapping_iterator++) {
      const string target_name = target_mapping_iterator->name();
      const unsigned target_index = target_mapping_iterator->index();
      target_identifiers.insert(pair<string,unsigned>(target_name, target_index));
    }

    // populate the query identifiers
    google::protobuf::RepeatedPtrField<const goby::IdentifierInfo>::const_iterator query_mapping_iterator;
    for (query_mapping_iterator = header.query_name_mapping().mappings().begin(); query_mapping_iterator != header.query_name_mapping().mappings().end(); query_mapping_iterator++) {
      string query_name = query_mapping_iterator->name();
      const unsigned query_index = query_mapping_iterator->index();
      query_identifiers.insert(pair<string,unsigned>(query_name, query_index));
    }

    // open the "entries" file
    const string entries_filename = basename + ".entries";
    entries_fd = ::open(entries_filename.c_str(), O_RDONLY | O_BINARY);
    if (entries_fd < 0) {
      cerr << "Error opening file: " << entries_filename << endl;
    }

    alignment_entry_iterator_end = new AlignmentEntryIterator(entries_fd, static_cast<streamoff>(0), ios_base::end);
  }

  AlignmentReader::~AlignmentReader(void) {
    delete alignment_entry_iterator_end;
  }

  AlignmentEntryIterator AlignmentReader::begin() const {
    return AlignmentEntryIterator(entries_fd);
  };

  AlignmentEntryIterator AlignmentReader::end() const {
    return *alignment_entry_iterator_end;
  };

  AlignmentWriter::AlignmentWriter(const string& basename) : Alignment(basename) {
  }

  AlignmentWriter::AlignmentWriter(const Alignment& alignment) : Alignment(alignment) {
    // TODO: testing only
    this->basename = "foo";
  }

  AlignmentWriter::~AlignmentWriter(void) {
  }

  void AlignmentWriter::write() {
    // Write to the "header" file
    const string headerFilename = basename + ".header";
    cout << "Writing file: " << headerFilename << endl;
    int fd = ::open(headerFilename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0644);

    // set up a gzip output stream to compress the header
    google::protobuf::io::FileOutputStream headerFileStream(fd);
    google::protobuf::io::GzipOutputStream gzipHeaderStream(&headerFileStream);

    if (!header.SerializeToZeroCopyStream(&gzipHeaderStream)) {
      cerr << "Failed to write alignment header file: " << headerFilename << endl;
    }

    gzipHeaderStream.Close();
    // TODO? ::close(fd);
  }

#ifdef _MSC_VER
#pragma warning(pop)  // Restores the warning state.
#endif
}
