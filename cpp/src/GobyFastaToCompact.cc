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

#include <iostream>
#include <string>
#include <zlib.h>

#include "goby/common.h"
#include "goby/MessageChunks.h"
#include "goby/Reads.pb.h"
#include "goby/Reads.h"

#include "kseq.h"     // see http://lh3lh3.users.sourceforge.net/parsefastq.shtml

using namespace std;

// declare the type of file handler and the read() function for the fasta parser
KSEQ_INIT(gzFile, gzread);

int main (int argc, const char *const argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <filename>" << endl;
    return -1;
  }
  
  const string fastaFilename = argv[1];

  // open the fasta file handler and initialize kseq
  gzFile fastaFp = gzopen(fastaFilename.c_str(), "r");
  kseq_t *seq = kseq_init(fastaFp);

  // set up a goby writer
  goby::ReadsWriter writer("foobar.txt");

  int l;
  while ((l = kseq_read(seq)) >= 0) {       // read sequences
    // the sequence itself
    writer.setSequence(seq->seq.s);

    // the description line
    string description(seq->name.s);
    if (seq->comment.l) {
      description.append(" ");
      description.append(seq->comment.s);
      writer.setDescription(description.c_str());
    }

    // the quality scores
    if (seq->qual.l) {
      writer.setQualityScores(seq->qual.s);
    }

    writer.appendEntry();
  }

  writer.close();                   // write the any remaining data
  kseq_destroy(seq);                // destroy seq  
  gzclose(fastaFp);                 // close the fasta file handler  

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

	return 0;
}

