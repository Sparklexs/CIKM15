/*
 * intersection.cpp
 *
 *  Created on: Oct 16th, 2014
 *      Author: Qi
 */

#ifndef INTERSECTION_H_
#define INTERSECTION_H_

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

class intersection{
private:
	unsigned int* pages;

public:
	intersection(unsigned int* pgs) : pages(pgs) {}
	// void operator()(CluewebReader* Reader, lptrArray& lps);
	void operator()(CluewebReader* Reader, 	const std::vector<std::string>& word_l);
	vector<pnode> findintersection(vector<pnode> A, vector<pnode> B); //get the intersection
	vector<pnode> getpnodelist(string term, CluewebReader* Reader);   //get the impact-sorted list

};

#endif /* INTERSECTION_H_ */