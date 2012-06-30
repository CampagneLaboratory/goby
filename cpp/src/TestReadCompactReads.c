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

/*
 * TestReadCompactReads.c
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */
#include <stdio.h>
#include "goby/GsnapSequence.h"
#include "goby/GsnapReads.h"

#define GSNAP
#undef PMAP

static CReadsHelper *readsHelper = (CReadsHelper *) NULL;

void gmap_dump_sequence(Sequence_T seq);
void gmap_dump_str(char *name, char *strval);
void gmap_dump_int(char *name, int intval);

int main(int argc, char **argv) {
	if (argc == 0) {
		printf("You must specify a compact-reads file to read");
	}
	argv++;
	argc--;
	printf("%d %s\n", argc, argv[0]);
	readsHelper = gobyReads_openReadsReader(argv, argc, 0);
	while (gobyReads_hasNext(readsHelper) == 1) {
		Sequence_T query = gobyReads_next(readsHelper);
		gmap_dump_sequence(query);
	}
	printf("Iterated %d times\n", readsHelper->numRead);
	printf("Calling finished\n");
	gobyReads_finished(readsHelper);
	return 0;
}

void gmap_dump_sequence(Sequence_T seq) {
	printf("\n\nDumping New Sequence  [[");
	gmap_dump_str("acc", seq.acc);
	gmap_dump_str("restofheader", seq.restofheader);
	gmap_dump_str("contents", seq.contents);
	gmap_dump_str("contents_alloc", seq.contents_alloc);
	gmap_dump_int("fulllength", seq.fulllength);
#ifdef GSNAP
	printf("  <GSNAP Structure>\n");
	gmap_dump_str("contents_uc", seq.contents_uc);
	gmap_dump_str("contents_uc_alloc", seq.contents_uc_alloc);
	gmap_dump_str("chop", seq.chop);
	gmap_dump_int("choplength", seq.choplength);
	gmap_dump_str("quality", seq.quality);
	gmap_dump_str("quality_alloc", seq.quality_alloc);
#endif
	gmap_dump_int("trimstart", seq.trimstart);
	gmap_dump_int("trimend", seq.trimend);
#ifdef PMAP
	printf("  <PMAP Structure>\n");
	gmap_dump_int("fulllength_given", seq.fulllength_given);
#endif
	gmap_dump_int("subseq_offset", seq.subseq_offset);
	gmap_dump_int("skiplength", seq.skiplength);
	printf("]]\n");
}

void gmap_dump_str(char *name, char *strval) {
	if (strval == NULL) {
		printf("  %s == NULL\n", name);
	} else {
		printf("  %s[l=%d]%s\n", name, strlen(strval), strval);
	}
}

void gmap_dump_int(char *name, int intval) {
	printf("  %s == %d\n", name, intval);
}
