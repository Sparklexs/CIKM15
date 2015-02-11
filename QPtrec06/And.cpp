/*
 * And.cpp
 * Generate index for the pairs, only two-term query (or one) will come in
 *  Created on: Oct 16th, 2014
 *      Author: Qi
 */

#include "And.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include "CluewebReader.h"
#include <fstream>
using namespace std;

struct pinfo{
	unsigned int did;
	unsigned int freq;
	float score;
};

bool myfc (pinfo a, pinfo b){
	return (a.score > b.score);
}

void And::operator() (CluewebReader* Reader, lptrArray& lps) {

	assert(("Only expect 2-term queries now", lps.size() == 2));

	cout<<lps[0]->term<<" "<<lps[0]->unpadded_list_length<<endl;
	cout<<lps[1]->term<<" "<<lps[1]->unpadded_list_length<<endl;
	const int topK = min(lps[0]->unpadded_list_length, lps[1]->unpadded_list_length);
	// const int topK = 10;
	
	QpResult res[topK];

	// initial sorting by did
	lps.sort();

	int smallest_did = lps[0]->did;
	int af = 1;
	float final_score = 0;

	// initialize results heap
	for ( int i = 0; i < topK; ++i)  {
		res[i].did = -1;
		res[i].score = -1.0;
	}

	while(true) {
		// cout<<"inhere?"<<endl;
		// set round's threshold
		const float threshold = res[topK-1].score;
		
		// string terms[2];
		// float scores[2];
		// int freq[2];

		// scores[0] = 0;
		// scores[1] = 0;
		// freq[0] = 0;
		// freq[1] = 0;

		float frequency = 0;
		float score = 0;

		// check for termination condition
		if (smallest_did >= CONSTS::MAXD)
			break;

		// initialize final score
		final_score = 0.0f;

		//check if all lists are in the intersection
		for (int i=0; i<lps.size(); ++i){
			if (lps[i]->did != smallest_did) {
				af = 0;
				break;	
			}
		}
		// cout<<"af: "<<af<<endl;
		// if not, advance all the lists, enter another round
		if(af == 0){
			for (int i=0; i<lps.size(); ++i){
			 	lps[i]->did = lps[i]->nextGEQ( smallest_did + 1 );
		 }
		 	lps.sort();
		 	smallest_did = lps[0]->did;
		 	// cout<<"smallest_did: "<<smallest_did<<endl;
		 	// int test;
		 	// cin>>test;
		 	af = 1;
		 	continue;
		}

		// evaluate all dids with did == smallest_did
		for (int i=0; i<lps.size(); ++i){

			// terms[i] = lps[i]->term;

			if (lps[i]->did == smallest_did) {
				//PROFILER(CONSTS::EVAL);
				//PROFILER(CONSTS::GETFREQ);
				const float frequency = lps[i]->getFreq();
				const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
				// scores[i] = score;
				// freq[i] = frequency;
				final_score += score;

				// move safely to the next did
				lps[i]->did = lps[i]->nextGEQ( smallest_did + 1 );
//				cout<<"next did: "<<lps[i]->did<<endl;
				//PROFILER(CONSTS::NEXTGEQ);
			}
		}

		// if calculated score more than threshold, heapify
		if (Fcompare(final_score, threshold)==1) {
			//PROFILER(CONSTS::HEAPIFY);
			int j;
			for (j = topK-2; (j >= 0) && (Fcompare(final_score, res[j].score)==1); j--)
				res[j+1]=res[j];
			// res[j+1].setRQi(smallest_did,final_score,terms,scores,freq);
			res[j+1].setR(smallest_did,final_score);
		}

		// sort list by did
		lps.sort();

		// set new smallest did for the next round
		smallest_did = lps[0]->did;
	} //end while


	/*Just to print out all the lists with docids for test only*/
	for(int i=0; i<2; i++){
		RawIndexList Rlist = Reader->load_raw_list(lps[i]->term, Reader->term_map[lps[i]->term]);
	    cout<<lps[i]->term<<": "<< Reader->term_map[lps[i]->term]<<" "<<Reader->listLen_[Reader->term_map[lps[i]->term]]<<endl;
	    vector<pinfo> t_list;

	   	for(int h = 0; h < Rlist.unpadded_list_length; h++){
			pinfo t_p;
			t_p.did  = Rlist.doc_ids.at(h);
			t_p.freq = Rlist.freq_s.at(h);
			t_p.score = Rlist.scores.at(h);
			t_list.push_back(t_p);
			// cout<<t_p.did<<" "<<t_p.freq<<" "<<t_p.score<<endl;
			// int ts;
			// cin>>ts;
		}
		sort(t_list.begin(), t_list.end(), myfc);

		ofstream index;
		index.open("/data/qw376/11.10/"+ Rlist.term, ofstream::app);
		for(int j = 0; j < Rlist.unpadded_list_length; j++){
			// index<<t_list.at(j).did<<endl;
			index<<t_list.at(j).did<<" "<<t_list.at(j).score<<endl;
			// index<<t_list.at(j).did<<" "<<t_list.at(j).freq<<" "<<t_list.at(j).score<<endl;
		}
		index.close();
	}
	/*Just to print out all the lists with docids for test only*/

	// cout<<lps[0]->term<<" "<<lps[1]->term<<endl;
	// for(int i=0; i<topK; i++){
	// 	cout<<res[i].did<<": "<<res[i].score<<endl;
	// }

	int position_in_topk;
	for (position_in_topk = topK-1; position_in_topk >=0; --position_in_topk)  {
		//std::cout << position_in_topk << "\t"<< res[position_in_topk].score << "\t" << res[position_in_topk].did << std::endl;
		if (( res[position_in_topk].score > -1.0) && (res[position_in_topk].did < CONSTS::MAXD+1))
			break;
	}

	ofstream index_depth;
	const string dir = "/data/qw376/intersection_index/" + lps[0]->term + "++" + lps[1]->term;
	index_depth.open(dir.c_str(), ofstream::app);
	for(int i=0; i<position_in_topk; i++){
		index_depth<<res[i].did<<": "<<res[i].score<<endl;
		// index_depth<<res[i].did<<endl;
	}
	index_depth.close();

}