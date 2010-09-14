/*
 * TestReadCompactReads.c
 *
 *  Created on: Sep 13, 2010
 *      Author: kdorff
 */
#include <stdio.h>
#include "goby/GsnapReads.h"

int main(int argc, char **argv) {
	// char *filename = "C:\\Users\\kdorff\\Desktop\\NGA\\small.compact-reads";
	// char *filename = "/desktop/NGA/small.compact-reads";
	char *filename = "C:\\Users\\kdorff\\Desktop\\NGA\\color-HBR.73_20080930_1_Ambion_HBR_HBR_F3-1.compact-reads";
	printf("Calling openReadsReader with %s\n", filename);
	CReadsHelper *readsHelper = openReadsReader(filename);
	int i = 0;
	while (hasNext(readsHelper) == 1) {
		next(readsHelper);
		i++;
	}
	printf("Iterated %d/%d times\n", i, readsHelper->numRead);
	printf("Calling finished\n");
	finished();
	return 0;
}
