/*
 * algo_init.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: Qi
 */

#include "algo_init.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <string>
#include <bitset>

 using namespace std;

 bool myfunc (const sinfo& a, const sinfo& b){
 	return (a.score > b.score);
 }

 bool sortcdt (const fullinfo& a, const fullinfo& b){
 	return (a.score > b.score);
 }	

 bool cmp_by_value(const PAIR& lhs, const PAIR& rhs) {  
 	return lhs.second < rhs.second;  
 }  

 void algo_init::operator() (CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res, profilerC& p) {

	//number of pairs we have
    number_of_pairs = pls.lengths.size();
	//number of single terms we have
    number_of_singles = lps.size();
	//total structures we have
    total_number_of_structures = pls.lengths.size() + lps.size();

	// cout<<"origin singleterms: "<<lps.size()<<endl;
	// cout<<"origin pairs: "<<pls.lengths.size()<<endl;

 	slice_docidspace();

 	for(int i=0; i<=slices; i++){
 		cout<<"slice # "<<i<<": "<<slice_offsets[i]<<endl;
 	}


 	load_singlelists_tomap(Reader, lps);

	decide_hashtablesize();

	decide_termbits();

	load_pairlists_tomap(pls);

	cout<<"singleterms after purge: "<<singleinfos.size()<<endl;
	cout<<"pairs after purge: "<<pairinfos.size()<<endl;

	p.start(CONSTS::ALLQS);

	taat_merging();

	p.end(CONSTS::ALLQS);

	writeout_results();

}


algo_init::algo_init(unsigned int* pgs){

	pages = pgs;

	/*To determine the hashtable size*/
	totaldids = 0;

	//For a 4 term query, the filter is 00001111 
 	short filter = (1<<number_of_singles) - 1;
	// bitset<8> y(filter);
	// cout<<y<<endl;

	/*we estimate the final number of results won't be larger than 15000*/
	fresults.reserve(15000);

	/*to devide the pairs*/
	dem = "+";

	elem = 0;

}

void algo_init::slice_docidspace(){

	slice_offsets[0] = 0;
 	slice_offsets[slices] = CONSTS::MAXD;
 	int range = CONSTS::MAXD/slices;

 	for(int i=1; i<slices; i++){
 		slice_offsets[i] = slice_offsets[i-1] + range;
 	}

}

void algo_init::load_singlelists_tomap(CluewebReader* Reader, lptrArray& lps){

	vector<sinfo> temps_v;
 	vector<size_t> ssize;

 	for (int i=0; i<lps.size(); ++i){

 		const string term = lps[i]->term;
 		const int length = lps[i]->unpadded_list_length;

		// termbits[term] = 0; //for termbits, just insert the term in the map, the value 0 doesn't matter
		term_orders.push_back(make_pair(term, length)); //for termbits, pair value is the listlength


		const int l_sqrt = sqrt(lps[i]->unpadded_list_length);
		// cout<<term<<": "<<length<<endl;

		if( lps.size()>2 && length >= 10000000){
			cout<<"Query size greater than 3, throw away singles larger than 10m: "<< term <<endl;
			continue;
		}

		int term_id = Reader->term_map[term];
		// cout<<term<<" "<<term_id<<endl;
		RawIndexList Rlist = Reader->load_raw_list(term,term_id);

		vector<sinfo> t_list;

		for(int h = 0; h < length; h++){
			sinfo t_p;
			t_p.did  = Rlist.doc_ids.at(h);
	 		// t_p.freq = Rlist.freq_s.at(h);
			t_p.score = Rlist.scores.at(h);
			t_list.push_back(t_p);
		}

		sort(t_list.begin(), t_list.end(), myfunc);

		size_t cut = 1000;
		if(length > cut)
	 		t_list.resize(cut);// resize according to Configuration file
		totaldids += cut;

	 	// cout<<t_list.size()<<endl;

		ssize.push_back(0); //the first should be 0;

		for(int h=0; h<slices; h++){
			int ct = 0;
			for(int j = 0; j < t_list.size(); j++){
				if(t_list.at(j).did > slice_offsets[h] && t_list.at(j).did <= slice_offsets[h + 1]){
					temps_v.push_back(t_list.at(j));
					ct ++;
				}
			}
	 	  	// cout<<"ct: "<<ct<<endl;
			ssize.push_back(ct+ssize.back());
		}



		singleinfos[term] = temps_v;
		singleslicesizes[term] = ssize;
		cout<<term<<": "<<temps_v.size()<<endl;
		cout<<term<<": "<<ssize.size() - 1<<endl;

		temps_v.clear();
		ssize.clear();

	}

}

void algo_init::load_pairlists_tomap(pairlists& pls){

	/*to parse term1+term2 in a pair*/
	position = 0;

	vector<pinfo> tempp_vu;
 	vector<pinfo> tempp_vo;
 	vector<pinfo> tempp_v;
 	vector<size_t> psize;

 	string queryline;
	string did_s;
	int did;
	string sc1_s;
	float sc1;
	string sc2_s;
	float sc2;
	int count = 0;
  	// string ts_s;
  	// float ts;

  	if(pls.lengths.size()>0){/*if pair list size is not 0*/

		for(int k = 0; k < pls.lengths.size(); k++){

			position = pls.pairnames.at(k).find(dem);
			term1 = pls.pairnames.at(k).substr(0, position);
			term2 = pls.pairnames.at(k).substr(position+1,pls.pairnames.at(k).size());
			// bitset<8> x(termbits[term1]);
			// bitset<8> y(termbits[term2]);
			// bitset<8> z(((~termbits[term1])>>3) | ((~termbits[term2])>>3));
			// cout<<"in loading pairs stage: "<<term1<<": "<<x<<", "<<term2<<": "<<y<<": "<<z<<endl;

			if( ( ((~termbits[term1])>>3) | ((~termbits[term2])>>3) ) > 0)
				continue;

			ifstream pair_stream;
			string pair_dir = CONSTS::pair_index + pls.pairnames.at(k);
	  		// string pair_dir = "/data/qw376/pair_index/" + pls.pairnames.at(k);
			pair_stream.open(pair_dir.c_str());

			count = 0;

	  	// cout<<pls.pairnames.at(k)<<endl;

			while(getline(pair_stream, queryline)){

				pinfo temp;

	  		// cout<<queryline<<endl;
				string::iterator itr = queryline.begin();
				string::iterator start = itr;

				while(itr != queryline.end() && !isspace(*itr)){
					++itr;
				}

				did_s = string(start, itr);
				did = atoi(did_s.c_str());

			//take the total score
				start = itr+1;
				itr++;
				while(itr != queryline.end() && !isspace(*itr)){
					++itr;
				}
			// ts_s = string(start, itr);
			// ts = atof(ts_s.c_str());

			//ignore the first freq
				start = itr+1;
				itr++;
				while(itr != queryline.end() && !isspace(*itr)){
					++itr;
				}

			//take the first score
				start = itr+1;
				itr++;
				while(itr != queryline.end() && !isspace(*itr)){
					++itr;
				}

				sc1_s = string(start, itr);
				sc1 = atof(sc1_s.c_str());

			//ignore the second freq
				start = itr+1;
				itr++;
				while(itr != queryline.end() && !isspace(*itr)){
					++itr;
				}

			//take the second score
				start = itr+1;
				itr++;
				while(itr != queryline.end() && !isspace(*itr)){
					++itr;
				}

				sc2_s = string(start, itr);
				sc2 = atof(sc2_s.c_str());


				temp.did = did;
			// temp.s = ts;
				temp.s1 = sc1;
				temp.s2 = sc2;

				tempp_vo.push_back(temp);

				count ++;		
			// t_map[didl] = count;
				if (count == pls.lengths.at(k)){
					totaldids += count;
					break;
				}

			}

		  	psize.push_back(0); //the first should be 0;

		  		// cout<<pls.pairnames.at(k)<<": "<<tempp_vo.size()<<endl;

		  	for(int j=0; j<slices; j++){
		  		int ct = 0;
		  		for(int i=0; i<tempp_vo.size(); i++){	
		  			if(tempp_vo.at(i).did > slice_offsets[j] && tempp_vo.at(i).did <= slice_offsets[j + 1]){
		  				tempp_v.push_back(tempp_vo.at(i));
		  				ct ++;
		  			}
		  		}
			 	  	// cout<<"ct: "<<ct<<endl;
		  		psize.push_back(ct+psize.back());
		  	}


		  	pairinfos[pls.pairnames.at(k)] = tempp_v;
		  	pairslicesizes[pls.pairnames.at(k)] = psize;
		  	cout<<pls.pairnames.at(k)<<": "<<tempp_v.size()<<endl;
		  	cout<<pls.pairnames.at(k)<<": "<<psize.size() - 1<<endl;

		  	tempp_vo.clear();
		  	tempp_v.clear();
		  	psize.clear();

		  	pair_stream.close();
	  }
	  position = 0;//later will use position again to parse the pairs
	}//if pair list size is not 0
}

void algo_init::decide_hashtablesize(){
	if( 5*totaldids/slices < 8191)
		ht = initHash(totaldids/slices, 0); 
	else
		ht = initHash(8191, 1); //table size 8191
}


void algo_init::decide_termbits(){

	int orders = 0;
	sort(term_orders.begin(), term_orders.end(), cmp_by_value);

	for(int i=0; i<term_orders.size(); i++){
		cout<<term_orders.at(i).first<<": "<<term_orders.at(i).second<<endl;
		termbits[term_orders.at(i).first] = ~ (1<<orders++);
	}

	for(map<string, int>::iterator it = termbits.begin(); it!=termbits.end(); ++it){
		bitset<8> x(it->second);
		cout<<it->first<<": "<<x<<endl;

	}
}

void algo_init::taat_merging(){

		for(int k=0; k<slices; k++){

			/*for the pair lists hashing*/
			/*the reason put pairs before singles is that pairs is loaded later than singles in loading process, kind of cache warm-up*/
			for(map<string, vector<pinfo>>::iterator it = pairinfos.begin(); it!=pairinfos.end(); ++it){
				position = it->first.find(dem);
				term1 = it->first.substr(0, position);
				term2 = it->first.substr(position+1,it->first.size());
				for(int i=pairslicesizes[it->first].at(k); i<pairslicesizes[it->first].at(k+1); i++){

					int pos = insertHash(ht, it->second.at(i).did, elem, 0, fresults);
	    			if (pos){ //key already existed

	    				fresults.at(ht->table[pos-1]-1).score = fresults.at(ht->table[pos-1]-1).score + it->second.at(i).s1 * ( fresults.at(ht->table[pos-1]-1).kbits & (termbits[term1]) != filter) + it->second.at(i).s2 * ( fresults.at(ht->table[pos-1]-1).kbits & (termbits[term2]) != filter);
	    				fresults.at(ht->table[pos-1]-1).kbits = (fresults.at(ht->table[pos-1]-1).kbits) & (termbits[term1]) & (termbits[term2]);

	   				}else{ //successfully inserted

					 	fullinfo ftemp;
					 	ftemp.did = it->second.at(i).did;
					 	ftemp.score = it->second.at(i).s1 + it->second.at(i).s2;
					 	ftemp.kbits = (termbits[term1]) & (termbits[term2]);
					 	fresults.push_back(ftemp);

					 	elem++;
	   				}
	   			}
			}/*pairs finished*/

			/*for the single lists hashing*/
			for(map<string, vector<sinfo>>::iterator it = singleinfos.begin(); it!=singleinfos.end(); ++it){
				for(int i=singleslicesizes[it->first].at(k); i<singleslicesizes[it->first].at(k+1); i++){

					int pos = insertHash(ht, it->second.at(i).did, elem, 0, fresults);
	    			if (pos){ //key already existed

	    				fresults.at(ht->table[pos-1]-1).score = fresults.at(ht->table[pos-1]-1).score + it->second.at(i).score;
	    				fresults.at(ht->table[pos-1]-1).kbits = fresults.at(ht->table[pos-1]-1).kbits & (termbits[it->first]);

					}else{ //successfully inserted

					 	fullinfo ftemp;
					 	ftemp.did = it->second.at(i).did;
					 	ftemp.score = it->second.at(i).score;
					 	ftemp.kbits = termbits[it->first];
					 	fresults.push_back(ftemp);

					 	elem++;
					}
				}
			}/*singles finished*/

	   		clearHash(ht);
   	}

	sort(fresults.begin(), fresults.end(), sortcdt);
	fresults.resize(200);
}

void algo_init::writeout_results(){

	// cout<<"final did #: "<<didresults.size()<<endl;
	cout<<"final did #: "<<fresults.size()<<endl;
	cout<<"final did #: "<<elem<<endl;
	// cout<<"attempts: "<<hit<<endl;
	// cout<<"rate: "<<(float)hit/(float)offset<<endl;
	// cout<<"final did #: "<<fresults.size()<<endl;

	// for(int i=0; i<didresults.size(); i++){
	// 	bitset<8> x(kbitsresults.at(i));
	// 	cout<<didresults.at(i)<<", "<<scoreresults.at(i)<<", "<<x<<endl;
	// }

	int count = 0;
	int lookups = 0;
	// int hit = 0;

	ofstream out_stream;
	// out_stream.open(CONSTS::Candidates_200.c_str(), ofstream::app);


	for(int i=0; i<fresults.size(); i++){

		int n = ~(fresults.at(i).kbits);
		count = 0;
		while (n>0) { 
			count = count + (n&1);
			n=n>>1; //Right shift by 1 
		}
		lookups += number_of_singles - count; 
    		bitset<8> x(fresults.at(i).kbits);
			cout<<fresults.at(i).did<<", "<<fresults.at(i).score<<", "<<x<<", "<<number_of_singles - count<<endl;
			// out_stream<<fresults.at(i).did<<" ";
	}
	// out_stream<<endl;
	// out_stream.close();

	cout<<"Total lookups needed: "<<lookups<<endl;


}
