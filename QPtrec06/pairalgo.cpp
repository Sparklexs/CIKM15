/*
 * pairalgo.cpp
 *
 *  Created on: July 17, 2014
 *      Author: Qi
 */

 /*To generate the score lists*/

#include "pairalgo.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>
#include <unistd.h>
 
using namespace std;

struct scores{
	float s1;
	float s2;
	int f1;
	int f2;
};

struct pinfo{
	unsigned int did;
	unsigned int freq;
	float score;
};

bool myfunc (pinfo a, pinfo b){
	return (a.score > b.score);
}

void pairalgo::operator() (CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res) {

	vector<map<ulong, scores>> pimv;
	vector<vector<pinfo>> layer_lists;
	const int buffer_size = 100*1024;
	cout<<"buffer_size: "<<buffer_size<<endl;

	string queryline;
	ulong results[10];
	int known_l[10];
	for(int i = 0; i < 10; i++){
		known_l[i] = 0;
	}

    char * docid_s;
	// int docid;
	int count = 0;
	int test;

	/*---Get the top 10 results from Complicated ranking func----*/

	// ifstream input_stream;
	// // input_stream.open("../../../Dropbox/WSDM_Index_Script/result_log");  //for or
	// input_stream.open("../../../Dropbox/WSDM_Index_Script/500_cr");  //for cr

	// for(int i=0; i < qn-1; ++i){
 //        input_stream.ignore(numeric_limits<streamsize>::max(),'\n');
 //    }

	// getline(input_stream, queryline);

	// char * ql = new char [queryline.length()+1];
 //    std::strcpy (ql, queryline.c_str());

	// // cout<<queryline<<endl;


	// docid_s = strtok (ql," ");

	// while (docid_s != NULL)
 //  	{
 //    	// printf ("%s\n",docid_s);
 // 		results[count] = atol(docid_s);   	
 // 		count++;
 //    	docid_s = strtok (NULL, " ,.-");
 //    	if(count == 10)
 //    		break;
 //  	}

 //  	input_stream.close();

  	// for(int b=0; b<count; b++){
  	// 	cout<<results[b]<<endl;
  	// }
  	/*-----------------------*/

  	// string did_s;
  	// ulong didl;
  	// string s1_s;
  	// float s1;
  	// string s2_s;
  	// float s2;
  	// string freq1_s;
  	// int freq1;
  	// string freq2_s;
  	// int freq2;


  	/*---load the pair lists into map and pimv vector from files, and cal known lvl for pair lists----*/
 //  	for(int k = 0; k < pls.lengths.size(); k++){

 //  	map<ulong, scores> t_map;
 //  	scores t_score;

 //  	ifstream pair_stream;
 //  	string pair_dir = "../../../Dropbox/pair_ind/" + pls.pairnames.at(k);
 //  	pair_stream.open(pair_dir.c_str());

 //  	count = 0;
 //  	cout<<pls.pairnames.at(k)<<endl;

 //  	while(getline(pair_stream, queryline)){

 //  		// cout<<queryline<<endl;
 //  		string::iterator itr = queryline.begin();
	// 	string::iterator start = itr;

 //   		 while(itr != queryline.end() && !isspace(*itr)){
	// 		++itr;
	// 	}

	// 	did_s = string(start, itr);
	// 	didl = atol(did_s.c_str());


	// 	start = itr+1;
 //  	  	itr++;
 //    	while(itr != queryline.end() && !isspace(*itr)){
	// 		++itr;
	// 	}

	// 	s1_s = string(start, itr);
	// 	s1 = atof(s1_s.c_str());
	// 	t_score.s1 = s1;

	// 	start = itr+1;
 //  	  	itr++;
	// 	while(itr != queryline.end() && !isspace(*itr)){
	// 		++itr;
	// 	}

	// 	freq1_s = string(start, itr);
	// 	freq1 = atoi(freq1_s.c_str());
	// 	t_score.f1 = freq1;

	// 	start = itr+1;
 //  	  	itr++;
	// 	while(itr != queryline.end() && !isspace(*itr)){
	// 		++itr;
	// 	}

	// 	s2_s = string(start, itr);
	// 	s2 = atof(s2_s.c_str());
	// 	t_score.s2 = s2;

	// 	start = itr+1;
 //  	  	itr++;
 //  	  	while(itr != queryline.end() && !isspace(*itr)){
	// 		++itr;
	// 	}

	// 	freq2_s = string(start, itr);
	// 	freq2 = atoi(freq1_s.c_str());
	// 	t_score.f2 = freq2;

	// 	t_map[didl] = t_score;

 //  		count ++;
 //  		if (count == pls.lengths.at(k))
 //  			break;

 //  	}
 //  		pimv.push_back(t_map);
	// 	// for (map<ulong,scores>::iterator it=t_map.begin(); it!=t_map.end(); ++it)
 //  //  			 cout << it->first << " => " << it->second.s1 << endl;

 //  		cout<<t_map.size()<<endl;

 //  		for(int i = 0; i < 10; i++){

 //  			map<ulong, scores>:: iterator it;
	// 		it = t_map.find(results[i]);
	// 		if(it != t_map.end()){
	// 			// cout<<results[i]<<" is known in "<<pls.pairnames.at(k)<<endl;
	// 			known_l[i] ++;
	// 		}

 //  		}


 //  	  	pair_stream.close();
	// }
	/*-----------------------*/

	/*cal known lvl for single terms*/
	for (int i=0; i<lps.size(); ++i){

		// if(lps[i]->term!="bahamas")
		// 	continue;

		const int l_sqrt = sqrt(lps[i]->unpadded_list_length);
		cout<<lps[i]->term<<": "<<lps[i]->unpadded_list_length<<endl;
		// cout<<lps[i]->term<<" "<<lps[i]->unpadded_list_length<<endl;
		const string term = lps[i]->term;
		int term_id = Reader->term_map[term];
		cout<<term<<" "<<term_id<<endl;
	 	RawIndexList Rlist = Reader->load_raw_list(term,term_id);

	 	vector<pinfo> t_list;

	 	for(int h = 0; h < lps[i]->unpadded_list_length; h++){
	 		pinfo t_p;
	 		t_p.did  = Rlist.doc_ids.at(h);
	 		t_p.freq = Rlist.freq_s.at(h);
	 		t_p.score = Rlist.scores.at(h);
	 		t_list.push_back(t_p);
	 	}


	 	sort(t_list.begin(), t_list.end(), myfunc);

	 	/*list less than a block*/
	 // 	if(lps[i]->unpadded_list_length<=buffer_size){

	 // 	const int size = lps[i]->unpadded_list_length;
	 // 	cout<<"p1"<<endl;
	 // 	unsigned int t_a [size];
	 // 	cout<<"p2"<<endl;

	 // 	for(int j = 0; j < lps[i]->unpadded_list_length; j++){
	 // 		t_a[j] = t_list.at(j).did;
	 // 	}

	 // 	const string path = "/home/qi/Dropbox/score_index_binary/"+lps[i]->term;
	 // 	FILE *f_b = fopen(path.c_str(),"wb");
	 // 	cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), size, f_b)<<endl;
	 // 	fflush(f_b);
	 // 	fclose(f_b);


	 // 	ofstream index;
	 // 	index.open("/home/qi/Dropbox/score_index/"+lps[i]->term, ofstream::app);
	 // 	for(int j = 0; j < lps[i]->unpadded_list_length; j++){
		// 		// index<<t_a[j]<<endl;
	 // 			index<<t_list.at(j).did<<endl;
		// }
		// index.close();
	 // }

	 /*list more than a block*/
	 // if(lps[i]->unpadded_list_length>buffer_size){

	 // 	 int block = lps[i]->unpadded_list_length / buffer_size;
	 // 	 const int rem = lps[i]->unpadded_list_length % buffer_size;
	 // 	 cout<<"block: "<<block<<" rem: "<<rem<<endl;

		// const string path = "/home/qi/Dropbox/score_index_binary/"+lps[i]->term;
	 // 	FILE *f_b = fopen(path.c_str(),"wb");

	 // 	 for(int k = 0; k < block; k++){
	 // 	 	unsigned int t_a [buffer_size];
	 // 	 	for(int j = 0; j < buffer_size; j++){
	 // 			t_a[j] = t_list.at(j+k*buffer_size).did;
	 // 		}
	 // 	 	cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), buffer_size, f_b)<<endl;
	 // 	 }	
	 	 
	 // 	 if(rem != 0){
	 // 	 unsigned int t_a [rem];
	 // 	 for(int j = 0; j < rem; j++){
	 // 		t_a[j] = t_list.at(j+block*buffer_size).did;
	 // 	 }
	 // 	 cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), rem, f_b)<<endl;
	 // 	}
	 // }

	 /*---------*/
	 	ofstream index;
	 	index.open("/data/qw376/singletermlists/"+lps[i]->term, ofstream::app);
	 	for(int j = 0; j < lps[i]->unpadded_list_length; j++){
				// index<<t_list.at(j).did<<endl;
	 			index<<t_list.at(j).did<<" "<<t_list.at(j).freq<<" "<<t_list.at(j).score<<endl;
		}
		index.close();

	}
	/*-----------------------*/
}
