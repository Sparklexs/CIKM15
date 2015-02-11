/*
 * pairrank.cpp
 * generate depth for bigrams
 *  Created on: Oct 21, 2014
 *      Author: Qi
 */

#include "pairrank.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>
using namespace std;

void pairrank::operator() (int qn, pairlists& pls) {

	// vector<map<ulong, int>> pimv;
	// vector<map<ulong, int>> timv;
	const int p_l = pls.cutoffs.size();
	// cout<<qn<<" "<<"pairs: "<<pls.cutoffs.size()<<endl;

	ofstream pair_num;
	pair_num.open(CONSTS::num_of_pairs.c_str(), ofstream::app);
	pair_num<<pls.cutoffs.size()<<endl;
	pair_num.close();

	int p_e[10][p_l]; 

	const int buffer_size = 1024*100;
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

	ifstream input_stream;
	// input_stream.open("../../../Dropbox/WSDM_Index_Script/result_log");  //for or
	// input_stream.open("/data/qw376/WSDM_Index_Script/500_cr");  //for cr
	input_stream.open(CONSTS::cr_results.c_str());  //for cr


	for(int i=0; i < qn-1; ++i){
        input_stream.ignore(numeric_limits<streamsize>::max(),'\n');
    }

	getline(input_stream, queryline);

	char * ql = new char [queryline.length()+1];
    std::strcpy (ql, queryline.c_str());

	// cout<<queryline<<endl;


	docid_s = strtok (ql," ");

	while (docid_s != NULL)
  	{
    	// printf ("%s\n",docid_s);
 		results[count] = atol(docid_s);   	
 		count++;
    	docid_s = strtok (NULL, " ,.-");
    	if(count == 10)
    		break;
  	}

  	input_stream.close();

  	// for(int b=0; b<count; b++){
  	// 	cout<<results[b]<<endl;
  	// }
  	/*-----------------------*/


	/*---load the pair term lists into map and outputs the depth (From Binary lists)----*/
	if(pls.cutoffs.size()>0){
		for(int k = 0; k<pls.cutoffs.size(); k++){

			ifstream index_score;
			string pair_dir = CONSTS::intersection_index_did + pls.pairnames[k];

			const int size = pls.cutoffs[k];

			// cout<<pls.pairnames[k]<<": "<<size<<endl;

			/*if list length is less than a block*/
			if (size <= buffer_size){
				index_score.open(pair_dir.c_str(), ios::binary);
				unsigned int t_a[size];
				index_score.read((char*)t_a, sizeof(unsigned int)*size);

	  		for(int i = 0; i < 10; i++){
	  				int j;
	  				for(j = 0; j<size; j++){
	  					if (results[i] == t_a[j])
	  						break;
	  				}

					if(j!=size){
						// cout<<results[i]<<" is at depth: "<<j+1<<" of "<<lps[k]->term<<endl;
						p_e[i][k] = j+1;
					
					}else{
						// cout<<results[i]<<" is at depth: -"<<" of "<<lps[k]->term<<endl;
						p_e[i][k] = 0;
					}

	  			}
	  			index_score.close();
	  		}
	  		/*----*/

	  		/*if the list length is more than a block*/
	  		if (size >= buffer_size){

	  			const int block =  size / buffer_size;
				const int rem = size % buffer_size;
				// cout<<"block: "<<block<<endl;
				// cout<<"rem: "<<rem<<endl;
				int depths[10];
				for(int m = 0; m < 10; m++){
					depths[m] = -1;
				}

				//interate the 10 results
	  			for(int i = 0; i < 10; i++){

	  				index_score.open(pair_dir.c_str(), ios::binary);

	  				int depth = -1;

	  				//iterate all the blocks
	  				for (int h = 0; h < block; h++){

	  					// cout<<"p1"<<endl;
	  					unsigned int t_a[buffer_size];
	  					// cout<<"p2"<<endl;
	  					index_score.read((char*)t_a, sizeof(unsigned int)*buffer_size);

	  					int j;
	  					for (j = 0; j < buffer_size; j++){
	  						// cout<<t_a[j]<<": "<<depth<<endl;
	  						if(results[i] == t_a[j]){
	  							depth = j + buffer_size * h + 1;
	  							break;
	  						}
	  					}

	  					if(depth > 0) {
	  						depths[i] = depth;
	  						break;
	  					}
	  					// int test;
	  					// cin>>test;
	  				}

	  				//iterate remaining parts
	  				if(depth == -1){

	  				   unsigned int t_a[rem];

	  				   index_score.read((char*)t_a, sizeof(unsigned int)*rem);

	  				   int j;
	  					for (j = 0; j < rem; j++){
	  						if(results[i] == t_a[j]){
	  							depth = j + buffer_size * block;
	  							break;
	  						}
	  					}

	  					if(depth > 0) {
	  						depths[i] = depth;
	  					}
	  				}

	  				index_score.close();
	  			}

	  			for(int m = 0; m < 10; m++){
					if(depths[m] == -1){
						// cout<<results[m]<<" is at depth: -"<<" of "<<lps[k]->term<<endl;
						p_e[m][k] = 0;
					}else{
						// cout<<results[m]<<" is at depth: "<<depths[m]<<" of "<<lps[k]->term<<endl;
						p_e[m][k] = depths[m] + 1;
					}
				}
	  		

	  		}

	  		/*----*/

		}

	// ofstream index_depth;
	// index_depth.open(CONSTS::pair_depth.c_str(), ofstream::app);
	// for(int i = 0; i < 10; i++){
	// 	for(int n = 0; n < p_l; n++){
	// 		// index_depth<<pls.pairnames[n]<<" "<<p_e[i][n]<<" ";
	// 		index_depth<<pls.cutoffs[n]<<" "<<p_e[i][n]<<" ";
	// 	}
	// 	index_depth<<endl;
	// }
	// index_depth.close();

	/*used for cutoff learning*/
	ofstream index_depth;
	index_depth.open(CONSTS::pair_depth.c_str(), ofstream::app);
	for(int n = 0; n < p_l; n++){
		// cout<<pls.cutoffs[n]<<endl;
		index_depth<<pls.cutoffs[n]<<" ";
		for(int i = 0; i < 10; i++){
			// index_depth<<pls.pairnames[n]<<" "<<p_e[i][n]<<" ";
			// index_depth<<pls.cutoffs[n]<<" "<<p_e[i][n]<<" ";
			index_depth<<p_e[i][n]<<" ";
		}
	}
	index_depth<<endl;
	index_depth.close();

	/*used for pair depth index*/
	// ofstream index_depth;
	// index_depth.open(CONSTS::pair_depth.c_str(), ofstream::app);
	// for(int i = 0; i < 10; i++){
	// 	for(int n = 0; n < p_l; n++){
	// 		index_depth<<pls.pairnames.at(n)<<" "<<p_e[i][n]<<" ";
	// 	}
	// 	index_depth<<endl;
	// }
	// index_depth.close();
  }

}
