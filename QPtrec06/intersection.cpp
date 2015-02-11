/*
 * intersection.cpp
 * Generate index for the pairs, only two-term query (or one) will come in
 *  Created on: Oct 16th, 2014
 *      Author: Qi
 */

#include "intersection.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include "CluewebReader.h"
#include <fstream>
using namespace std;

bool myfctn (pnode a, pnode b){
	return (a.score > b.score);
}

vector<pnode> intersection::findintersection(vector<pnode> A, vector<pnode> B){
	vector<pnode> intersection_list;
	int n1 = A.size();
	int n2 = B.size();
	int i = 0, j = 0;
	while (i < n1 && j < n2) {
		if (A[i].did > B[j].did) {
		  j++;
		} else if (B[j].did > A[i].did) {
		  i++;
		} else {
		  A[i].score2 = B[j].score;
		  A[i].score = A[i].score + B[j].score;
		  // A[i].freq = A[i].freq + B[j].freq;
		  intersection_list.push_back(A[i]);
		  i++;
		  j++;
		}
	}
	return intersection_list;
}

vector<pnode> intersection::getpnodelist(string term, CluewebReader* Reader){

	RawIndexList Rlist = Reader->load_raw_list(term, Reader->term_map[term]);
    // cout<<term<<": "<< Reader->term_map[term]<<" "<<Reader->listLen_[Reader->term_map[term]]<<endl;
    // cout<<term<<": "<<Reader->listLen_[Reader->term_map[term]]<<endl;
    vector<pnode> t_list;

   	for(int h = 0; h < Rlist.unpadded_list_length; h++){
		pnode t_p;
		t_p.did  = Rlist.doc_ids.at(h);
		t_p.freq = Rlist.freq_s.at(h);
		t_p.score = Rlist.scores.at(h);
		t_p.score1 = Rlist.scores.at(h);
		t_p.score2 = 0;
		t_list.push_back(t_p);
		// cout<<t_p.did<<" "<<t_p.freq<<" "<<t_p.score<<endl;
		// int ts;
		// cin>>ts;
	}
	// sort(t_list.begin(), t_list.end(), myfctn);
	return t_list;
}

/*To generate pair lists in ASCII*/
// void intersection::operator() (CluewebReader* Reader, lptrArray& lps) {

// 	assert(("Only expect 2-term queries now", lps.size() == 2));	
// 	vector<pnode> listofterm1 = getpnodelist(lps[0]->term, Reader);
// 	vector<pnode> listofterm2 = getpnodelist(lps[1]->term, Reader);
// 	// vector<pnode> listofterm1 = getpnodelist("butterfly", Reader);
// 	// vector<pnode> listofterm2 = getpnodelist("pewter", Reader);
// 	vector<pnode> intersection_list = findintersection(listofterm1, listofterm2);
// 	sort(intersection_list.begin(), intersection_list.end(), myfctn);
// 	// cout<<"intersection size: "<<intersection_list.size()<<endl;

// 	if(intersection_list.size() > 0) {
// 		ofstream index;
// 		// index.open("/data/qw376/intersection_index/"+ lps[0]->term + "+" + lps[1]->term, ofstream::app);
// 		index.open("/home/vgc/qi/intersection_index/"+ lps[0]->term + "+" + lps[1]->term, ofstream::app);
// 		for(int j = 0; j < intersection_list.size(); j++){
// 			index<<intersection_list.at(j).did<<endl;
// 			// index<<intersection_list.at(j).did<<" "<<intersection_list.at(j).score<<endl;
// 			// index<<intersection_list.at(j).did<<" "<<intersection_list.at(j).freq<<" "<<intersection_list.at(j).score<<endl;
// 		}
// 		index.close();
// 	}
// }
/*To generate pair lists in ASCII*/

/*To generate pair lists in Binary Format*/
// void intersection::operator() (CluewebReader* Reader, lptrArray& lps) {
void intersection::operator() (CluewebReader* Reader, const std::vector<std::string>& word_l) {

	// const int buffer_size = 100*1024;
	// assert(("Only expect 2-term queries now", lps.size() == 2));
	assert(("Only expect 2-term queries now", word_l.size() == 2));	
	// vector<pnode> listofterm1 = getpnodelist(lps[0]->term, Reader);
	// vector<pnode> listofterm2 = getpnodelist(lps[1]->term, Reader);
	vector<pnode> listofterm1 = getpnodelist(word_l[0], Reader);
	vector<pnode> listofterm2 = getpnodelist(word_l[1], Reader);
	// vector<pnode> listofterm1 = getpnodelist("butterfly", Reader);  
	// vector<pnode> listofterm2 = getpnodelist("pewter", Reader);
	vector<pnode> intersection_list = findintersection(listofterm1, listofterm2);
	// sort(intersection_list.begin(), intersection_list.end(), myfctn);
	// cout<<"intersection size: "<<intersection_list.size()<<endl;

	ofstream index;
	index.open(CONSTS::pair_lex, ofstream::app);
	// index<<lps[0]->term + "+" + lps[1]->term<<" "<<intersection_list.size()<<endl;
	index<<word_l[0] + "+" + word_l[1]<<" "<<intersection_list.size()<<endl;
	index.close();
	
	// /*list less than a block*/
	// if(intersection_list.size()<=buffer_size){

	// 	const int size = intersection_list.size();
	// 	unsigned int t_a [size];
	// 	// float t_f1 [size];
	// 	// float t_f2 [size];

	// 	for(int j = 0; j < intersection_list.size(); j++){
	// 		t_a[j] = intersection_list[j].did;
	// 		// t_f1[j] = intersection_list[j].score1;
	// 		// t_f2[j] = intersection_list[j].score2;
	// 	}

	// 	// const string path = fd + lps[0]->term + "+" + lps[1]->term;
	// 	const string path_d = CONSTS::intersection_index_did + lps[0]->term + "+" + lps[1]->term;
	// 	// const string path_s1 = CONSTS::intersection_index_score_1 + lps[0]->term + "+" + lps[1]->term;
	// 	// const string path_s2 = CONSTS::intersection_index_score_2 + lps[0]->term + "+" + lps[1]->term;
	// 	FILE *f_d = fopen(path_d.c_str(),"wb");
	// 	// FILE *f_s1 = fopen(path_s1.c_str(),"wb");
	// 	// FILE *f_s2 = fopen(path_s2.c_str(),"wb");

	// 	// FILE *f_b = fopen(CONSTS::intersection_index.c_str(),"wb");
	// 	// cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), size, f_b)<<endl;
	// 	fwrite(t_a, sizeof(unsigned int), size, f_d);
	// 	fflush(f_d);
	// 	fclose(f_d);
	// 	// fwrite(t_f1, sizeof(float), size, f_s1);
	// 	// fflush(f_s1);
	// 	// fclose(f_s1);
	// 	// fwrite(t_f2, sizeof(float), size, f_s2);
	// 	// fflush(f_s2);
	// 	// fclose(f_s2);

	// }

 // 	/*list more than a block*/
	// if(intersection_list.size()>buffer_size){

	// 	int block = intersection_list.size() / buffer_size;
	// 	const int rem = intersection_list.size() % buffer_size;
	// 	// cout<<"block: "<<block<<" rem: "<<rem<<endl;

	// 	const string path_d = CONSTS::intersection_index_did + lps[0]->term + "+" + lps[1]->term;
	// 	// const string path_s1 = CONSTS::intersection_index_score_1 + lps[0]->term + "+" + lps[1]->term;
	// 	// const string path_s2 = CONSTS::intersection_index_score_2 + lps[0]->term + "+" + lps[1]->term;
	// 	FILE *f_d = fopen(path_d.c_str(),"wb");
	// 	// FILE *f_s1 = fopen(path_s1.c_str(),"wb");
	// 	// FILE *f_s2 = fopen(path_s2.c_str(),"wb");

	// 	for(int k = 0; k < block; k++){
	// 	  unsigned int t_a [buffer_size];
	// 	  // float t_f1 [buffer_size];
	// 	  // float t_f2 [buffer_size];

	// 	  for(int j = 0; j < buffer_size; j++){
	// 		t_a[j] = intersection_list[j+k*buffer_size].did;
	// 		// t_f1[j] = intersection_list[j+k*buffer_size].score1;
	// 		// t_f2[j] = intersection_list[j+k*buffer_size].score2;
	// 	  }
	// 	  // cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), buffer_size, f_b)<<endl;
	// 	  fwrite(t_a, sizeof(unsigned int), buffer_size, f_d);
	// 	  // fwrite(t_f1, sizeof(float), buffer_size, f_s1);
	// 	  // fwrite(t_f2, sizeof(float), buffer_size, f_s2);
	// 	}	

	// 	if(rem != 0){

	// 		unsigned int t_a [rem];
	// 		// float t_f1 [rem];
	// 		// float t_f2 [rem];

	// 			for(int j = 0; j < rem; j++){
	// 				t_a[j] = intersection_list[j+block*buffer_size].did;
	// 				// t_f1[j] = intersection_list[j+block*buffer_size].score1;
	// 				// t_f2[j] = intersection_list[j+block*buffer_size].score2;
	// 			}
	// 		// cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), rem, f_b)<<endl;
	// 			fwrite(t_a, sizeof(unsigned int), rem, f_d);
	// 			// fwrite(t_f1, sizeof(float), rem, f_s1);
	// 			// fwrite(t_f2, sizeof(float), rem, f_s2);
	// 	}
	// 	fclose(f_d);
	// 	// fclose(f_s1);
	// 	// fclose(f_s2);

 //    }


}
/*To generate pair lists in Binary Format*/