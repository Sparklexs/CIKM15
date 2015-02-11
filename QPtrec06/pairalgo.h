/*
 * pairalgo.h
 *
 *  Created on: July 17, 2014
 *      Author: Qi
 */

#ifndef PAIRALGO_H_
#define PAIRALGO_H_

#include "PostingOriented_BMW.h"
#include "CluewebReader.h"

class pairalgo{
private:
	unsigned int* pages;

public:
	pairalgo(unsigned int* pgs) : pages(pgs) {}
	void operator()(CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topk, QpResult* r);
};

#endif /* PAIRALGO_H_ */