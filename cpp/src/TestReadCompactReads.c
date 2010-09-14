/*
 * TestReadCompactReads.c
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */
#include <stdio.h>
#include "goby/GsnapSequence.h"
#include "goby/GsnapReads.h"

int main(int argc, char **argv) {
	// char *filename = "C:\\Users\\kdorff\\Desktop\\NGA\\small.compact-reads";
	// char *filename = "/desktop/NGA/small.compact-reads";
	char *filename = "C:\\Users\\kdorff\\Desktop\\NGA\\color-HBR.73_20080930_1_Ambion_HBR_HBR_F3-1.compact-reads";
	printf("Calling openReadsReader with %s\n", filename);
	CReadsHelper *readsHelper = gobyReads_openReadsReader(filename);
	while (gobyReads_hasNext(readsHelper) == 1) {
		Sequence_T query = gobyReads_next(readsHelper);
		printf(">%d\n%s\n", (readsHelper->numRead - 1), query.contents);
	}
	printf("Iterated %d times\n", readsHelper->numRead);
	printf("Calling finished\n");
	gobyReads_finished(readsHelper);
	return 0;
}
