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
  Alignment::Alignment(const std::string& basename) {
    this->basename = basename;
    this->pbHeader = AlignmentHeader::default_instance();
  }

  Alignment::~Alignment(void) {
  }

  string Alignment::getBasename(const char* filename) {
    return getBasename(string(filename));
  }

  string Alignment::getBasename(const string& filename) {
    const string COMPACT_ALIGNMENT_FILE_EXTS[] = {
      ".entries", ".header", ".tmh", ".stats", ".counts"
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
    vector<unsigned> queryLengths(pbHeader.number_of_queries());
    if (hasConstantQueryLength()) {
      // The alignment has constant query lengths
      queryLengths.assign(pbHeader.number_of_queries(), getConstantQuerylength());
    } else {
      queryLengths.assign(pbHeader.query_length().begin(), pbHeader.query_length().end());
    }
  
    return queryLengths;
  }

  vector<unsigned> Alignment::getTargetLengths() const {
    vector<unsigned> targetLengths(pbHeader.target_length().begin(), pbHeader.target_length().end());
    return targetLengths;
  }
  
  AlignmentReader::AlignmentReader(const string& basename) : Alignment(basename) {
    // open the "header" file
    const string headerFilename = basename + ".header";
    const int fd = ::open(headerFilename.c_str(), O_RDONLY | O_BINARY);

    if (fd > 0) {
      // uncompress file into memory so that it can be parsed
      google::protobuf::io::FileInputStream headerFileStream(fd);
      google::protobuf::io::GzipInputStream gzipHeaderStream(&headerFileStream);

      // the stream may not get all read in at once so we may need to copy in chunks
      void* header = NULL;
      size_t headerSize = 0;

      const void* buffer;
      int bufferSize;
      while (gzipHeaderStream.Next(&buffer, &bufferSize)) {
        // store the end location of the header buffer
        int index = headerSize;

        // resize the header buffer to fit the new data just read
        headerSize += bufferSize;
        header = (void *)realloc(header, headerSize);

        // and append the new data over to the end of the header buffer
        ::memcpy(reinterpret_cast<char*>(header) + index, buffer, bufferSize);
      }

      // populate the alignment header object from the uncompressed data
      if (!pbHeader.ParseFromArray(header, headerSize)) {
        cerr << "Failed to parse alignment header file: " << headerFilename << endl;
      }

      // free up the temporary buffers
      ::free(header);

      // close the streams and files
      headerFileStream.Close();
    } else {
      cerr << "Failed to open alignment header file: " << headerFilename << endl;
    }

    // populate the target identifiers
    google::protobuf::RepeatedPtrField<const goby::IdentifierInfo>::const_iterator targetMappingIterator;
    for (targetMappingIterator = pbHeader.target_name_mapping().mappings().begin(); targetMappingIterator != pbHeader.target_name_mapping().mappings().end(); targetMappingIterator++) {
      const string targetName = targetMappingIterator->name();
      const unsigned targetIndex = targetMappingIterator->index();
      targetIdentifiers.insert(pair<string,unsigned>(targetName, targetIndex));
    }

    // populate the query identifiers
    google::protobuf::RepeatedPtrField<const goby::IdentifierInfo>::const_iterator queryMappingIterator;
    for (queryMappingIterator = pbHeader.query_name_mapping().mappings().begin(); queryMappingIterator != pbHeader.query_name_mapping().mappings().end(); queryMappingIterator++) {
      string queryName = queryMappingIterator->name();
      const unsigned queryIndex = queryMappingIterator->index();
      queryIdentifiers.insert(pair<string,unsigned>(queryName, queryIndex));
    }
    
    this->messageChunksIterator = new MessageChunksIterator<AlignmentCollection>(basename + ".entries");
  }

  AlignmentReader::~AlignmentReader(void) {
    delete messageChunksIterator;
  }

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
    int fd = ::open(headerFilename.c_str(), O_RDWR | O_CREAT | O_TRUNC | O_BINARY, 0644);

    // set up a gzip output stream to compress the header
    google::protobuf::io::FileOutputStream headerFileStream(fd);
    google::protobuf::io::GzipOutputStream gzipHeaderStream(&headerFileStream);

    if (!pbHeader.SerializeToZeroCopyStream(&gzipHeaderStream)) {
      cerr << "Failed to write alignment header file: " << headerFilename << endl;
    }

    gzipHeaderStream.Close();
    // TODO? ::close(fd);
  }

#ifdef _MSC_VER
#pragma warning(pop)  // Restores the warning state.
#endif
}
