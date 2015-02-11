/*
 * pairrank.h
 *
 *  Created on: Oct 21, 2014
 *      Author: Qi
 */

#ifndef PAIRRANK_H_
#define PAIRRANK_H_

#include "PostingOriented_BMW.h"

class pairrank{
private:
	unsigned int* pages;

public:
	pairrank(unsigned int* pgs) : pages(pgs) {}
	void operator()(int qn, pairlists& pls);
};

#endif /* PAIRRANK_H_ */
