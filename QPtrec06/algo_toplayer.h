/*
 * algo_toplayer.h
 *
 *  Created on: Nov 8, 2014
 *      Author: Qi
 */

#ifndef algo_toplayer_H_
#define algo_toplayer_H_

#include "PostingOriented_BMW.h"
#include "CluewebReader.h"
#include "hash.h"

 using namespace std;

 const static int slices = 200;

 typedef pair<string, int> PAIR;  

 struct pinfo{
 	unsigned int did;
 	float s1;
 	float s2;
 };

 struct sinfo{
 	unsigned int did;
	// unsigned int freq;
 	float score;
 };

 struct fullinfo{
 	int did;
 	float score;
 	short kbits;
 	short coverage;
 };

 bool myfunc (const sinfo& a, const sinfo& b);

 bool sortcdt (const fullinfo& a, const fullinfo& b);

 bool sort_by_coverage(const fullinfo& a, const fullinfo& b);

 bool cmp_by_value(const PAIR& lhs, const PAIR& rhs); 

class algo_toplayer{
private:
	unsigned int* pages;

	/*Assign each did an index, for the did -> index mapping*/
	int num_of_did_tracker;

	/*Num of lookups*/
	int num_of_lookups;

	/*Map for all the top layers*/
 	map<string, vector<sinfo> > TopLayer_Map;
 	map<string, int>TopLayer_Order_Map;
 	map<string, vector<sinfo> >::iterator iter;

 	/*Map for all the pairs*/
 	map<string, vector<pinfo> > TermPair_Map;
 	map<string, int>TermPair_Order_Map;
 	map<string, vector<pinfo> >::iterator iter_pair;


 	/*For each did there is an entry to access the vector, to make it faster, using a hashtable here*/
 	map<int, int> Did_to_Index_Map;

 	/*This is the accumulator for TAAT*/
 	vector<fullinfo> Accumulator_Vector;

	/*total # of docids from layers, need to be processed, may have dups*/
 	int totaldids_in_layers;

 	/*total # of docids from pairs, need to be processed, may have dups*/
	int totaldids_in_pairs;

  /*this vector is to store the term and its order based on their listlengths*/
 	vector<PAIR> term_orders;

  /*each term has a # with a binary bits representation, 0 indicates existing, sorted by impact, the lower the bit is, the higher the impact
	say cat has the shortest listlen, then the corresponding bits is 11111110*/
 	map<string, int> termbits_Map;

  //number of pairs we have
  int number_of_pairs;
	//number of single terms we have
  int number_of_singles;
	//total structures we have
  int total_number_of_structures;

  //11.9
  vector<uint> top_layer_listlengths;

  vector<uint> pairs_listlengths;  

  vector<uint> depths;

  vector<uint> depths_pairs;

  vector<uint> cutoffs;

  vector<uint> cutoffs_pairs;

  short filter_global;

 	/*to parse the pairs*/
 	string dem;
 	string term1, term2;
 	size_t position;
 	short pre_compute_kbits;

  /*hashtable*/
  hashTable *ht;


public:
	algo_toplayer(unsigned int* pgs);
	void operator()(const vector<vector<float> > &top_layer_model,
					const vector<vector<float> > &term_pair_model,
					const vector<uint> &posting_block_boundry,
					const vector<uint> &posting_block_sizes,
					const vector<uint> &listlen_block_boundry,
          const vector<uint> &intersection_block_boundry,
          CluewebReader* Reader, int qn, toplayers& tls,toplayers& otls, pairlists& pls, pairlists& opls,lptrArray& lps, const int topK, profilerC& p, int qid);
	void load_singlelists_tomap(toplayers& tls, toplayers& otls);
	void load_pairlists_tomap(pairlists& pls, pairlists& opls);
	void decide_termbits();
	void adjust_cutoffs();
	void adjust_cutoffs_pairs();
	void TAAT_merge();
	void Do_Lookups(lptrArray& lps);
	void writeout_results(int qid);
	short Get_Coverage_From_Kbits(short kbits);
	vector<short> Get_indexes_of_termlists_to_do_lookup(short kbits);

};


/*From Costas*/

class entry {
  public:
    // members
    uint idx;
    uint space;
    float score;
    
    // methods
    entry(uint, uint, float);
    entry();
    void setEntry(uint, uint, float);
};

bool compareScore(const entry& i, const entry& j);

uint findIdxInBucket(const vector<uint> vec, const uint input);

inline void onlineGreedyDepthSelectionAlgorithm(const uint space_budget,
                                        const vector<vector<float> > &model,
                                        const vector<uint> &row_block_boundary,
                                        const vector<uint> &col_block_boundary,
                                        const vector<uint> &posting_block_sizes,
                                        const vector<uint> &depths,
                                        vector<uint> &cutoffs) {

	// for(int i=0; i<model.size(); ++i){
	// 	for(int j=0; j<model[i].size(); ++j){
	// 		cout<<i<<" "<<j<<" "<<model[i][j]<<endl;
	// 	}
	// }

  vector<entry> entries (posting_block_sizes.size()*depths.size(), entry());
  // cout<<"depths size "<<depths.size()<<endl;

  // Compute score and add to heap.
  int cnt = 0;
  for (size_t i = 0; i < depths.size(); ++i) {

    if(depths[i]==0) //the corner case for toplayers with 0 cutoffs
      continue;
  	
    size_t bucket_idx = findIdxInBucket(row_block_boundary , depths[i]);
    size_t cur_col = 0;

    while(cur_col < col_block_boundary.size() && model[bucket_idx][cur_col] != 0){
    	 // cout<<"adding term and position "<<i<<" with size "<<posting_block_sizes[cur_col]<<
    	 // " score "<< model[bucket_idx][cur_col]<<endl;
    	 entries[cnt++].setEntry(i, posting_block_sizes[cur_col], model[bucket_idx][cur_col]);
    	 // cout<<"bucket_idx: "<<bucket_idx<<" term_indx: "<<i<<" cur_col "<<cur_col<<endl;
    	 ++cur_col;

    }

    // for (size_t j = 0; j < col_block_boundary.size() && model[bucket_idx][j] != 0; ++j) {
    //   // Note the leftovers, but we avoid ifs.
    //   entries[cnt++].setEntry(i, posting_block_sizes[j], model[bucket_idx][j]);
    //   // cout<<i<<" "<<j<<" "<<model[bucket_idx][j]<<endl;
    // }
  }
  // cout<<"after loop "<<cnt<<endl;
  // cout<<"entry size: "<<entries.size()<<endl;
  entries.resize(cnt);
  if(cnt==0){ //the corner case for all the toplayers with 0 cutoffs, cutoffs are initialized to 0
    return;
  }
  // Sort entries by score.
  // cout<<"after resize size is : " << entries.size()<<endl;
  sort(entries.begin(), entries.end(), compareScore); 
  // cout<<"after sort"<<endl;
  // Execute the greedy algorithm. Note fix leftovers ?
  uint current_budget = 0;
  uint idx = 0;
  while (current_budget < space_budget && idx < entries.size() ) {
  	// cout<<"current_budget: "<<current_budget<<" idx: "<<idx<<endl;
  	uint new_budget = current_budget + entries[idx].space;
  	// cout<<"new_budget: "<<new_budget<<" idx: "<<idx<<endl;
  	if(new_budget <= space_budget){
  		// cout<<"new_budget2: "<<new_budget<<" idx: "<<idx<<endl;

  		current_budget = new_budget;
  		cutoffs[entries[idx].idx] += entries[idx].space;
  		// cout<<"here"<<endl;
  	}else{
  		// cout<<"new_budget3: "<<new_budget<<" idx: "<<idx<<endl;
  		current_budget += space_budget - current_budget; 
  		cutoffs[entries[idx].idx] += space_budget - current_budget;
  	}
  	++idx;
  }

  // for(int i = 0; i<cutoffs.size(); i++){
  // 	cout<<i<<" cutoffs "<<cutoffs[i]<<" < "<<space_budget<<endl;
  // }
}
// Assumption: the cutoffs include the toplayer followed by the term pair ones.
//         	so the size of cutoff vector must be depths.size() + depths_pairs.size().
 inline void onlineGreedyDepthSelectionAlgorithmUnify(const uint space_budget,
                                           	const vector<vector<float> > &top_layer_model,
                                          	const vector<uint> &listlen_block_boundary,
                                           	const vector<uint> &posting_block_boundary,
                                          	const vector<uint> &posting_block_sizes,
                                          	const vector<uint> &depths,
						const vector<uint> &top_layer_listlengths,
                                          	const vector<vector<float> > &term_pair_model,
                                          	const vector<uint> &intersection_block_boundary,
                                          	const vector<uint> &depths_pairs,
						const vector<uint> &pairs_listlengths,
                                          	vector<uint> &cutoffs_top_layer,
                                                vector<uint> &cutoffs_pairs) {
   // Initialize the entries to the maximum pair and toplayer structures we could have, multiplied
   // by the total chunks of ranks we have.
  vector<entry> entries (posting_block_sizes.size()*(depths.size() + depths_pairs.size()), entry());

   // Initialize temporary cutoff vector
  vector<uint> cutoffs( depths.size() + depths_pairs.size(), 0);
   // Compute score and add to heap for top layer.
  int cnt = 0;
  for (size_t i = 0; i < depths.size(); ++i) {
    if (depths[i] == 0) continue; // corner case for 0 cutoffs
    size_t bucket_idx = findIdxInBucket(listlen_block_boundary , top_layer_listlengths[i]);
    size_t cur_col = 0;
    while (cur_col < posting_block_boundary.size() &&
           top_layer_model[bucket_idx][cur_col] != 0 &&
           depths[i] > posting_block_boundary[cur_col]) {
      entries[cnt++].setEntry(i, posting_block_sizes[cur_col], top_layer_model[bucket_idx][cur_col]);
      ++cur_col;
     }
    }
    
// Compute score and add to heap for term-pairs.
   for (size_t i = 0; i < depths_pairs.size(); ++i) {
     if (depths_pairs[i] == 0) continue; // corner case for 0 cutoffs
       size_t bucket_idx = findIdxInBucket(intersection_block_boundary , pairs_listlengths[i]);
       size_t cur_col = 0;
	 
 	while (cur_col < posting_block_boundary.size() &&
               term_pair_model[bucket_idx][cur_col] != 0 &&
               depths_pairs[i] > posting_block_boundary[cur_col]) {
   	  entries[cnt++].setEntry(depths.size() + i, posting_block_sizes[cur_col], 5*term_pair_model[bucket_idx][cur_col]);
   	  ++cur_col;
 	}
  }

  entries.resize(cnt);
  if (cnt==0) return; // corner case for all toplayers with 0 cutoffs.

  // Sort entries by score.
  sort(entries.begin(), entries.end(), compareScore);

  uint current_budget = 0;
  uint idx = 0;
	
  // Execute the greedy algorithm.
  while (current_budget < space_budget && idx < entries.size() ) {
	uint new_budget = current_budget + entries[idx].space;
	if(new_budget <= space_budget){
  	current_budget = new_budget;
  	cutoffs[entries[idx].idx] += entries[idx].space;
	} else {
  	current_budget += space_budget - current_budget;
  	cutoffs[entries[idx].idx] += space_budget - current_budget;
	}
	++idx;
  }


  if (depths.size() == 0 && depths_pairs.size() > 0) {
    cutoffs_pairs = cutoffs;
  } else if (depths.size() > 0 && depths_pairs.size() == 0) {
    cutoffs_top_layer = cutoffs;
  } else if (depths.size() > 0 && depths_pairs.size() > 0) {
    for (size_t i =0; i < depths.size(); ++i) {
      cutoffs_top_layer[i] = cutoffs[i];
    }
    for (size_t i = 0, j = depths.size(); i < depths_pairs.size(); ++i, ++j) {
      cutoffs_pairs[i] = cutoffs[j];
    }
  } //else both depths and depths_pairs are empty, do nothing.
}

#endif /* algo_toplayer_H_ */
