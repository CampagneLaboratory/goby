/*
 * CReadsHelper.h
 *
 *  Created on: Sep 14, 2010
 *      Author: kdorff
 */

#ifndef CREADSHELPER_H_
#define CREADSHELPER_H_

#ifdef __cplusplus
	#include <queue>
	#include <string>
	#include "Reads.h"
	// More complex structure for C++
	struct CReadsHelper {
		goby::ReadsReader *readsReader;
		goby::ReadEntryIterator *it;
		goby::ReadEntryIterator *end;
		std::queue<std::string> *unopenedFiles;
		unsigned char circular;
		unsigned int numRead;
	};
#else
	// Opaque structure for C
	typedef struct {
		void *readsReader;
		void *it;
		void *end;
		void *unopenedFiles;
		unsigned char circular;
		unsigned int numRead;
	} CReadsHelper;
#endif

#endif /* CREADSHELPER_H_ */
