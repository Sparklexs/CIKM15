/*
 * rankinfo.h
 *
 *  Created on: Oct 9, 2014
 *      Author: Qi
 */

#ifndef RANKINFO_H_
#define RANKINFO_H_

#include "PostingOriented_BMW.h"

class rankinfo{
private:
	unsigned int* pages;

public:
	rankinfo(unsigned int* pgs) : pages(pgs) {}
	void operator()(int qn, lptrArray& lps);
};

#endif /* RANKINFO_H_ */
