/*
 * algo_init.h
 *
 *  Created on: Aug 13, 2014
 *      Author: Qi
 */

#ifndef ALGO_INIT_H_
#define ALGO_INIT_H_

#include "PostingOriented_BMW.h"
#include "CluewebReader.h"
#include "hash.h"

 using namespace std;

 const static int slices = 6;

 typedef pair<string, int> PAIR;  

 struct pinfo{
 	unsigned int did;
	// float s;
 	float s1;
 	float s2;
 };

 struct sinfo{
 	unsigned int did;
	// unsigned int freq;
 	float score;
 };

 struct scores{
 	float s1;
 	float s2;
 	int f1;
 	int f2;
 };

 struct fullinfo{
 	int did;
 	float score;
 	short kbits;
 };

 bool myfunc (const sinfo& a, const sinfo& b);

 bool sortcdt (const fullinfo& a, const fullinfo& b);

 bool cmp_by_value(const PAIR& lhs, const PAIR& rhs); 

class algo_init{
private:
	unsigned int* pages;

	/*singleinfos has the pinfo vector stores all the dids in slices seq, while in each slice ordered by inmpact*/
	/*singleslicesizes map has #=slice elem, each elem is the # of dids within each slice*/
 	map<string, vector<sinfo>> singleinfos;
 	map<string, vector<size_t>> singleslicesizes;

 	/*pairinfos has the pinfo vector stores all the dids in slices seq, while in each slice ordered by inmpact*/
	/*pairslicesizes map has #=slice elems, each elem is the # of dids within each slice*/
 	map<string, vector<pinfo>> pairinfos;
 	map<string, vector<size_t>> pairslicesizes;

    /*The offsets array for each slice*/
 	int slice_offsets[slices+1];

	/*total # of docids, used for determine the hashtable size*/
 	size_t totaldids;

    /*this vector is to store the term and its order based on their listlengths*/
 	vector<PAIR> term_orders;

    /*each term has a # with a binary bits representation, 0 indicates existing, sorted by impact, the lower the bit is, the higher the impact
	say cat has the shortest listlen, then the corresponding bits is 11111110*/
 	map<string, int> termbits;

    //number of pairs we have
    int number_of_pairs;
	//number of single terms we have
    int number_of_singles;
	//total structures we have
    int total_number_of_structures;

    //For a 4 term query, the filter is 00001111 
    short filter;

    /*fresults is to store the final merge results*/
	vector<fullinfo> fresults;

	/*hash table*/
	hashTable *ht;

	/*elem(offset) to insert into the hashtable*/
 	unsigned short elem;

 	/*to parse the pairs*/
 	string dem;
 	string term1, term2;
 	size_t position;


public:
	algo_toplayer(unsigned int* pgs);
	void operator()(CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res, profilerC& p);
	void slice_docidspace();
	void load_singlelists_tomap(CluewebReader* Reader, lptrArray& lps);
	void load_pairlists_tomap(pairlists& pls);
	void decide_hashtablesize();
	void decide_termbits();
	void taat_merging();
	void writeout_results();
	~algo_toplayer();
};


#endif /* ALGO_INIT_H_ */
