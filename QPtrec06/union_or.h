/*
 * union_or.cpp
 *
 *  Created on: Oct 16th, 2014
 *      Author: Qi
 */

#ifndef union_or_H_
#define union_or_H_

#include "PostingOriented_BMW.h"
#include "CluewebReader.h"
#include <assert.h>

struct pnode{
	unsigned int did;
	unsigned int freq;
	float score;
	float score1;
	float score2;
};

using namespace std;

class union_or{
private:
	unsigned int* pages;

public:
	union_or(unsigned int* pgs) : pages(pgs) {}
	void operator()(CluewebReader* Reader, lptrArray& lps);
	vector<pnode> findunion_or(vector<pnode> A, vector<pnode> B); //get the union_or
	vector<pnode> getpnodelist(string term, CluewebReader* Reader);   //get the impact-sorted list

};

#endif /* union_or_H_ */