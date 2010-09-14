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
	void *reader = openReadsReader(filename);
	void *it = getReadsIterator(reader);
	int i = 0;
	while (hasNext(reader, it, i) == 1) {
		next(it);
		i++;
	}
	printf("Iterated %d times\n", i);
	printf("Calling finished\n");
	void finished();
}
