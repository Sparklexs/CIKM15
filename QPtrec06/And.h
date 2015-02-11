/*
 * And.cpp
 *
 *  Created on: Oct 16th, 2014
 *      Author: Qi
 */

#ifndef AND_H_
#define AND_H_

#include "PostingOriented_BMW.h"
#include "CluewebReader.h"
#include <assert.h>

class And{
private:
	unsigned int* pages;

public:
	And(unsigned int* pgs) : pages(pgs) {}
	void operator()(CluewebReader* Reader, lptrArray& lps);
};

#endif /* AND_H_ */