//
// Copyright (C) 2009-2012 Institute for Computational Biomedicine,
//                         Weill Medical College of Cornell University
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#if	HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#include "goby/C_Gsnap.h"
#include "goby/C_Alignments.h"
#include "goby/C_Reads.h"
#include <pcrecpp.h>

using namespace	std;

char *read_file(char *filename)	{
	cout <<	filename << endl;
	FILE *pFile	= fopen(filename, "r" );
	if (pFile == NULL) {
		fputs("File	error",	stderr);
		exit(1);
	}

	// obtain file size:
	fseek(pFile, 0,	SEEK_END);
	long lSize = ftell(pFile);
	rewind(pFile);

	// allocate	memory to contain the whole	file + null	terminator for the string
	char *buffer = (char*) malloc((sizeof(char)	* lSize) + sizeof(char));
	if (buffer == NULL)	{
		fputs("Memory error", stderr);
		exit(1);
	}

	// copy	the	file into the buffer:
	size_t result =	fread(buffer, 1, lSize,	pFile);
	if (result != lSize) {
		fputs("Reading error", stderr);
		exit(1);
	}
	// null	terminate the buffer
	buffer[lSize] =	'\0';

	fclose (pFile);
	return buffer;
}

void registerChromosomes(CAlignmentsWriterHelper *writerHelper)	{
	char *test1	= read_file("test-data/zebrafish.chromosome.list.txt");
	gobyGsnap_test_registerTargets(writerHelper, test1);
	free(test1);
}

void targetIdentiferTest(CAlignmentsWriterHelper *writerHelper)	{
	cout <<	"MT?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "MT") <<	" ";
	cout <<	"Zv9...?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Zv9_scaffold3564") << "	";
	cout <<	"25?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "25") <<	" ";
	cout <<	"z?" <<	gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Z") << endl;

	unsigned t_mt =	gobyAlignments_targetIndexForIdentifier(writerHelper, "MT");
	unsigned t_zv9 = gobyAlignments_targetIndexForIdentifier(writerHelper, "Zv9_scaffold3564");
	unsigned t_25 =	gobyAlignments_targetIndexForIdentifier(writerHelper, "25");
	unsigned t_z = gobyAlignments_targetIndexForIdentifier(writerHelper, "Z");
	cout <<	"MT=" << t_mt << " " <<	"Zv9...=" << t_zv9 << "	" << "25=" << t_25 << "	" << "z=" << t_z <<	endl;

	// Verify that getting targetIndex didn't effect the value of the index
	cout <<	"MT?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "MT") <<	" ";
	cout <<	"Zv9...?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Zv9_scaffold3564") << "	";
	cout <<	"25?" << gobyAlignments_isTargetIdentifierRegistered(writerHelper, "25") <<	" ";
	cout <<	"z?" <<	gobyAlignments_isTargetIdentifierRegistered(writerHelper, "Z") << endl;
}

void testSingleForwardDel(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/gsnap-output-single-qual-del.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testPairedEndSpliced(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/gsnap-paired-spliced-multi.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}


void testPairedEnd1(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/gsnap-output-pair-test-1.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testPairedEnd2(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/gsnap-paired-multi.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testUnpaired1(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/unpaired-keep-left.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testUnpaired2(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/unpaired-keep-right.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}


void testSingleEndNoQual(CAlignmentsWriterHelper *writerHelper)	{
	char *test1	= read_file("test-data/gsnap-output-nonpair-noqual-test-1.gnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testSingleReverseInsQual(CAlignmentsWriterHelper *writerHelper) {
	char *test1	= read_file("test-data/gsnap-output-single-reverse-ins.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testSeqvarX(CAlignmentsWriterHelper *writerHelper,	int	which) {
	char buffer[256];
	sprintf(buffer,	"test-data/%02d-seq-var-reads-gsnap.gsnap",	which);
	char *test1	= read_file(buffer);
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
}

void testSeqvar() {
	CAlignmentsWriterHelper	*writerHelper;
	gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk("deleteme-seqvar", &writerHelper);
    gobyAlignments_setQualityAdjustment(writerHelper, -33);

	char *test1	= read_file("test-data/synth.chromosome.list.txt");
	gobyGsnap_test_registerTargets(writerHelper, test1);
	free(test1);

	int	x =	0;
	testSeqvarX(writerHelper, 2); x++;
	testSeqvarX(writerHelper, 4); x++;
	testSeqvarX(writerHelper, 7); x++;
	testSeqvarX(writerHelper, 12); x++;
	testSeqvarX(writerHelper, 13); x++;
	testSeqvarX(writerHelper, 14); x++;
	testSeqvarX(writerHelper, 24); x++;
	testSeqvarX(writerHelper, 25); x++;
	testSeqvarX(writerHelper, 26); x++;

	gobyAlignments_finished(writerHelper, x);
	goby_shutdownProtobuf();
}

void pcreTest()	{
	char *a	= "s-reads:3.2..b:5.1,l_reads:3,c:3,d:5,a:0..b:2.8";
	string temp1;
	string da;	   string daClip;	string daProb;
	string temp2;
	string std;	   string stdClip;

	pcrecpp::StringPiece input(a);	// Wrap	in a StringPiece
	pcrecpp::RE	re("(([A-Za-z_\\-]+):(\\d+)\\.(\\d))|(([A-Za-z_\\-]+):(\\d))");

	printf("pcre start\n");
	while (re.FindAndConsume(&input,
			&temp1,
			&da, &daClip, &daProb,
			&temp2,
			&std, &stdClip)) {

		if (da.length()	> 0) {
			cout <<	"da=" << da	<< " ";
			cout <<	"daClip=" << daClip	<< " ";
			cout <<	"daProb=" << daProb	<< endl;
		} else {
			cout <<	"std=" << std << " ";
			cout <<	"stdClip=" << stdClip << endl;
		}

	}
	printf("done\n");
}

void spliceTests(CAlignmentsWriterHelper *writerHelper)	{
	char *test1	= read_file("test-data/spliced-pair.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);

	test1 =	read_file("test-data/shortexon-forward.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);

	test1 =	read_file("test-data/shortexon-reverse.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);
	
	test1 =	read_file("test-data/splice-two-hits-reverse.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);

	test1 =	read_file("test-data/splice-sense.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);

	test1 =	read_file("test-data/splice-antisense.gsnap");
	gobyGsnap_parse(writerHelper, test1);
	free(test1);

}

void sequenceTests() {
	CAlignmentsWriterHelper	*writerHelper;
	gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk("deleteme-seq", &writerHelper);
    gobyAlignments_setQualityAdjustment(writerHelper, -33);
	registerChromosomes(writerHelper);
    testUnpaired1(writerHelper);
    testUnpaired2(writerHelper);
    testPairedEndSpliced(writerHelper);
	testPairedEnd1(writerHelper);
    testPairedEnd2(writerHelper);
	spliceTests(writerHelper);
	targetIdentiferTest(writerHelper);
	testSingleForwardDel(writerHelper);
	testSingleReverseInsQual(writerHelper);
	testSingleEndNoQual(writerHelper);
	gobyAlignments_finished(writerHelper, 1);
    goby_shutdownProtobuf();
}

int	main(int argc, const char *const argv[]) {
	pcreTest();
	sequenceTests();
	testSeqvar();
	return 0;
}
