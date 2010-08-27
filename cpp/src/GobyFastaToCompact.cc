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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#ifdef HAVE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

#include "goby/Reads.pb.h"
#include "goby/Reads.h"

#include "kseq.h"     // see http://lh3lh3.users.sourceforge.net/parsefastq.shtml

using namespace std;

// declare the type of file handler and the read() function for the fasta parser
KSEQ_INIT(gzFile, gzread);

string stripExtension(const string& filename) {
  // Probable file extensions for FASTA/FASTQ files.
  const string PROBABLE_FILE_EXTS[] = {
    ".fa.gz", ".fna.gz", ".fasta.gz", ".fq.gz", ".fnq.gz", ".fastq.gz",
    ".csfasta", ".csfasta.gz", ".csfastq", ".csfastq.gz",
    ".csfa", ".csfa.gz", ".csfq", ".csfq.gz",
    ".fa", ".fna", ".fasta", ".fq", ".fnq", ".fastq",
    ".txt", ".txt.gz", ".gz"
  };

  const size_t len = sizeof(PROBABLE_FILE_EXTS) / sizeof(PROBABLE_FILE_EXTS[0]);
  const size_t dotindex = filename.find_last_of(".");
  if (dotindex != string::npos) {
    const string extension = filename.substr(dotindex);
    for (int i = 0; i < len; i++) {
      if (extension.compare(PROBABLE_FILE_EXTS[i]) == 0) {
        return filename.substr(0, dotindex);
      }
    }
  }

  return filename;
}

int main(int argc, char **argv) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

#ifdef HAVE_BOOST_PROGRAM_OPTIONS
  // Declare the supported options.
  boost::program_options::options_description desc("Converts FASTA/FASTQ files to the Goby \"compact-reads\" file format");
  desc.add_options()
    ("help,h", "Displays usage help information")
    ("input", boost::program_options::value< vector<string> >(), "The input fasta files to convert to compact reads. The output files will have the same filename but end in .compact-reads. If the input file ends in .gz it will be decompressed on the fly.")
    ("include-descriptions,d", "When this switch is provided, include description lines into the compact output. By default, ignore description lines.")
    ("include-identifiers,x", "When this switch is provided, include identifiers into the compact output. By default, ignore identifiers. Identifiers are parsed out of description lines as the token before the first space or tab character.")
    ("exclude-sequences", "When this switch is provided, exclude sequences. This results in not writing sequences to the compact file. This can be useful to keep only an association between sequence index and identifier.")
    ("exclude-quality", "When this switch is provided, exclude quality scores. This results in not writing quality scores to the compact file.")
    ("sequence-per-chunk,n", boost::program_options::value<unsigned>()->default_value(GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK), "The number of sequences that will be written in each compressed chunk. Default is suitable for very many short sequences. Reduce to a few sequences per chunk if each sequence is very large.")
    ("output,o", boost::program_options::value<string>(), "If there is only one read file, this will force the output file to this specific filename. If there is more than one input file, the output filename will always be the input filename without the .fasta, .gz, etc.  extensions with an extension of .compact-reads. You should generally use an extension of .compact-reads when writing a compact reads file.")
    ("verbose-quality-scores", "Print quality scores to the console as they are read and converted to Phred score. Useful for testing with a small number of reads.");
    
  // it's ok to specify the input files as the last argument
  boost::program_options::positional_options_description p;
  p.add("input", -1);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    cout << "Usage: " << argv[0] << " input1 input2 ... inputN" << endl;
    cout << desc;
    return 0;
  }

  if (!vm.count("input")) {
    cerr << "Error: Parameter 'input' is required." << endl;
    cout << "Usage: " << argv[0] << " input1 input2 ... inputN" << endl;
    cout << desc;
    return -1;
  }

  // get options from command line
  const bool include_descriptions = vm.count("include-descriptions") > 0;
  const bool include_identifiers = vm.count("include-identifiers") > 0;
  const bool exclude_sequences = vm.count("exclude-sequences") > 0;
  const bool exclude_quality = vm.count("exclude-quality") > 0;
  const bool verbose_quality_scores = vm.count("verbose-quality-scores") > 0;
  const unsigned sequence_per_chunk = vm["sequence-per-chunk"].as<unsigned>();

  // iterate over the input filenames
  const vector<string> input_filenames = vm["input"].as< vector<string> >();
#else
  // not using boost program options so get filename from command line and default everything else
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <filename>" << endl;
    return -1;
  }

  // use default option values
  const bool include_descriptions = false;
  const bool include_identifiers = false;
  const bool exclude_sequences = false;
  const bool exclude_quality = false;
  const bool verbose_quality_scores = false;
  const unsigned sequence_per_chunk = GOBY_DEFAULT_NUMBER_OF_ENTRIES_PER_CHUNK;

  // iterate over the input filenames
  const vector<string> input_filenames(1, argv[1]);
#endif // HAVE_BOOST_PROGRAM_OPTIONS

  for (vector<string>::const_iterator it = input_filenames.begin() ; it != input_filenames.end(); it++) {
    const string input_filename = *it;

    string output_filename;
#ifdef HAVE_BOOST_PROGRAM_OPTIONS
    // if there is only one file to process and the user specified an output file
    if (input_filenames.size() == 1 && vm.count("output")) {
      // use the output file name
      output_filename = vm["output"].as<string>();
    } else {
      // create the output file name based on the input file name
#endif // HAVE_BOOST_PROGRAM_OPTIONS
      output_filename = stripExtension(input_filename) + ".compact-reads";
#ifdef HAVE_BOOST_PROGRAM_OPTIONS
    }
#endif // HAVE_BOOST_PROGRAM_OPTIONS

    cout << "Creating file " << output_filename << endl;

    // open the fasta file handler and initialize kseq
    gzFile input_fp = gzopen(input_filename.c_str(), "r");
    kseq_t *seq = kseq_init(input_fp);

    // set up a goby writer
    goby::ReadsWriter writer(output_filename, sequence_per_chunk);

    int l;
    while ((l = kseq_read(seq)) >= 0) {       // read sequences
      if (!exclude_sequences) {
        // the sequence itself
        writer.setSequence(seq->seq.s);
      }

      if (include_identifiers) {
        writer.setIdentifier(seq->name.s);
      }

      if (include_descriptions) {
        // the description line
        string description(seq->name.s);
        if (seq->comment.l) {
          description.append(" ");
          description.append(seq->comment.s);
          writer.setDescription(description.c_str());
        }
      }

      if (!exclude_quality) {
        // the quality scores
        if (seq->qual.l) {
          if (verbose_quality_scores) {
            cout << seq->qual.s << endl;
          }
          // TODO: convert quality scores based on user options
          writer.setQualityScores(seq->qual.s);
        }
      }

      writer.appendEntry();
    }

    writer.close();                   // write the any remaining data
    kseq_destroy(seq);                // destroy seq
    gzclose(input_fp);                // close the fasta file handler
  }

  // Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}

