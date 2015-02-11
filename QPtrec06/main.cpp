#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <math.h>

#include "TrecReader.h"
#include "CluewebReader.h"
#include "qp.h"
#include "optionparser.h"
using namespace std;

 enum  optionIndex { UNKNOWN, HELP, BUILD, ONDEMAND, OUTPATH, QP, QPATH, OFFSET, LIMIT, VERB, TOPK, LAYER, INDEXROOT, DATAANLY, SCOREIND };
 const option::Descriptor usage[] =
 {
  {UNKNOWN, 	0,"" , ""    ,		option::Arg::None, "USAGE: example [options]\n\nOptions:" },
  {HELP,    	0,"h", "help",		option::Arg::None, "  --help  \tPrint usage and exit." },
  {BUILD,   	0,"b", "buildindex",option::Arg::Optional, "  --buildindex[=<raw index path>], -b [<raw index path>]  \tBuild index with optional path to raw index" },
  {ONDEMAND,   	0,"d", "demandBld", option::Arg::Optional, "  --demandBld[=<raw index path>], -d [<raw index path>]  \tBuild term lists on demand with optional path to raw index" },
  {OUTPATH,   	0,"o", "outputpath",option::Arg::Optional, "  --outputpath[=<index path>], -o [<index path>]  \tWhere to store the index" },
  {QPATH,   	0,"p", "qpath",     option::Arg::Optional, "  --qpath[=<query log path>], -p [<query log path>]  \tWhere to get the queries from" },
  {QP,   		0,"q", "qp",		option::Arg::Optional, "  -q [<buckets>], \t--qp[=<buckets>] \tRun query proc. with optional buckets" },
  {LIMIT,  		0,"l", "limit",		option::Arg::Optional, "  -l [<limit>], \t--limit[=<limit>] \tRun at most limit queries" },
  {OFFSET, 		0,"s", "start",		option::Arg::Optional, "  -s [start], \t--start[=<start>] \tThe first index to start at -> for build, but when run selects the exp. block size" },
  {VERB, 		0,"v", "verbosity",	option::Arg::Optional, "  -v [int], \t--verbosity[=int] \tLogging level -- lower level - more messages" },
  {TOPK, 		0,"k", "topk",	    option::Arg::Optional, "  -k [int], \t-- number of top-k documents to retrieve" },
  {LAYER, 		0,"r", "layer",		option::Arg::Optional, "  -r [int], \t-- Layer off: 0, Layer on: 1" },
  {INDEXROOT, 	0,"i", "indexroot",	option::Arg::Optional, "  -i [str], \t--indexroot[=str] \tThe path to trecroot" },
  {DATAANLY,    0,"a", "dataanly",  option::Arg::None, "-a"},
  {SCOREIND,    0,"c", "scoreind",  option::Arg::None, "-c"},
  {0,0,0,0,0,0}
 };

struct pinfo{
	unsigned int did;
	unsigned int freq;
	float score;
};

struct dinfo{
	int kl;
	int listlength;
	int depth;
};

bool myf (pinfo a, pinfo b){
	return (a.score > b.score);
}

int main(int argc, char * argv[] ){
	   argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	   option::Stats  stats(usage, argc, argv);
	   option::Option options[stats.options_max], buffer[stats.buffer_max];
	   option::Parser parse(usage, argc, argv, options, buffer);

	   if (parse.error())
	     return 1;

	   if (options[HELP]) {
	     option::printUsage(std::cout, usage);
	     return 0;
	   }

	   if (options[VERB] && options[VERB].arg)
		   Log::setGlobalVerbosityForAllLoggers(atoi(options[VERB].arg));

	   int buckets = (options[QP] && options[QP].arg) ? atoi(options[QP].arg) : 0;
	   int limit = (options[LIMIT] && options[LIMIT].arg) ? atoi(options[LIMIT].arg) : CONSTS::MAXD;
	   int offset = (options[OFFSET] && options[OFFSET].arg) ? atoi(options[OFFSET].arg) : 0;
	   int topk = (options[TOPK] && options[TOPK].arg) ? atoi(options[TOPK].arg) : 10;
	   int layer = (options[LAYER] && options[LAYER].arg) ? atoi(options[LAYER].arg) : 0;

	   if(options[BUILD]) {
		   // const std::string&  trecRawRoot = options[BUILD].arg ? std::string(options[BUILD].arg) : CONSTS::trecRawRoot;
		   std::string outputPath = (options[OUTPATH] && options[OUTPATH].arg) ? std::string(options[OUTPATH].arg): CONSTS::trec_output;
		   stringIntVectorsPair tmap; //empty one
//		   TrecReader* Reader =  TrecFactory(trecRawRoot);
//		   TrecFactory(*Reader,outputPath, offset,limit,tmap);
//		   RawIndexList first_R = Reader->load_raw_list("planet",1);
//		   CompressedList first_C(first_R);
//
////		   first_C.serializeToFS(trec_output);
//		   SqlProxy sql(CONSTS::trec_output+"indx.sqlite");
//		   first_C.serializeToDb(CONSTS::trec_output, sql);
//		   delete Reader;
//		   return 0;


		   CluewebReader* Reader  = CluewebFactory();
		   Reader->load_baby_index();
		   delete Reader;

		   return 0;
	   }

	   // if(options[DATAANLY]) {

	   // 		 cout<<"Known Level Calculation.."<<endl;

	   // 		 int pair_th = 1000;    //pair list threshold

	   // 		 string infoline;
	   // 		 int linenum = 0;

	   // 		 CluewebReader* Reader  = CluewebFactory();

	   // 		 ifstream index_depth;
	   // 		 index_depth.open(CONSTS::index_depth);

	   // 		 ofstream index_kl;	
	   // 		 index_kl.open(CONSTS::index_kl, ofstream::app);   		 

	   // 		 while(getline(index_depth, infoline)){
	   // 		 	map<string, dinfo> kmap;
	   // 		 	// map<string, int>pmap;
	   // 		 	vector<string> tvector;
	   // 		 	istringstream iss(infoline);
	   // 		 	int known_lvl = 0;

	   // 		 	do
    // 			{
    //     			string sub;
    //     			iss >> sub;
    //     			if(sub!=""){
    //     			// cout << "Substring: " << sub << endl;
    //     			tvector.push_back(sub);
    //     		    }
    // 			} while (iss);

    // 			for(int i = 0; i<tvector.size(); i++){
    // 				// cout<<tvector.at(i)<<endl;
    // 				if(i % 2 == 0){// if it's term or pair

    // 					if(tvector.at(i).find('+')==string::npos){ //if it's not pair
    // 					dinfo t_d;
    // 					int term_id = Reader->term_map[tvector.at(i)];
    // 					t_d.listlength = Reader->listLen_[term_id];
    // 					t_d.kl = 0;
    // 					t_d.depth = 0;
    // 					kmap[tvector.at(i)] = t_d;
    // 					}else{//if it's a pair

    // 						//do things when depth come in

    // 					}
    // 				}else{ //if it's a depth
    // 					if(tvector.at(i-1).find('+')==string::npos){ //if it's not for a pair
    // 					kmap[tvector.at(i-1)].depth = atoi(tvector.at(i).c_str());
    
    // 					// if(kmap[tvector.at(i-1)].depth <= min(pow(pow(kmap[tvector.at(i-1)].listlength, 1.0/3), 2), (double)10000) && kmap[tvector.at(i-1)].depth != 0){  //single term threshold
    // 					// if(kmap[tvector.at(i-1)].depth <= kmap[tvector.at(i-1)].listlength && kmap[tvector.at(i-1)].depth != 0){
    // 					// if(kmap[tvector.at(i-1)].depth <= (min(sqrt(kmap[tvector.at(i-1)].listlength),(double)10000) + 1) && kmap[tvector.at(i-1)].depth != 0){
    // 					// if(kmap[tvector.at(i-1)].depth <= (min(kmap[tvector.at(i-1)].listlength,7087) + 1) && kmap[tvector.at(i-1)].depth != 0){
    // 					// if(kmap[tvector.at(i-1)].depth <= (min(0.2*kmap[tvector.at(i-1)].listlength, (double)10000) + 1) && kmap[tvector.at(i-1)].depth != 0){
    // 					if(kmap[tvector.at(i-1)].depth == -10){ //no singlelists
    // 						kmap[tvector.at(i-1)].kl = 1;
    // 					}

    // 				  } else {//if it's for a pair
    // 				  		if(atoi(tvector.at(i).c_str())<= pair_th && atoi(tvector.at(i).c_str())!=0){
    // 				  		// cout<<"parsing a pair:"<<endl;
    // 				  		string delimiter = "+";
    // 				  		string s = tvector.at(i-1);
				// 			size_t pos = 0;
				// 			string token;
				// 			while ((pos = s.find(delimiter)) != std::string::npos) {
   	// 							 token = s.substr(0, pos);
    // 							 // std::cout << token << std::endl;
    // 							 if(kmap[token].kl==0)
				// 				 kmap[token].kl = 1;
    // 							 s.erase(0, pos + delimiter.length());
				// 			}
				// 			// cout << s << endl;
				// 			if(kmap[s].kl==0)
				// 			kmap[s].kl = 1;

				// 		}
    // 				  }
    // 				}

    // 			}

    // 			for(map<string, dinfo>::iterator it = kmap.begin(); it!=kmap.end(); ++it){
    // 				// cout<<it->first<<": "<<it->second.kl<<" "<<it->second.listlength<<" "<<sqrt(it->second.listlength)<<" "<<it->second.depth<<endl;
    // 				known_lvl = known_lvl + it->second.kl;
    // 			}

    // 			linenum ++;

    // 			if(linenum % 10 == 1 && linenum != 1)
    // 			index_kl << endl;

    // 			if(linenum % 10 == 1)
    // 			index_kl << kmap.size()<< " ";

    // 			index_kl << known_lvl << " ";

    // 			// if(linenum == 20)
    // 			// 	break;


	   // 		 }

	   // 		 index_kl.close();
	   // 		 delete Reader;
	   // 		 return 0;
	   	
	   // }

	   /*building top layer impact sorted index*/

	   // if(options[SCOREIND]) {

		  //  // cout<<"Building Scoreind Here"<<endl;
		  //  // string fd("/home/qi/top_layer_index_did/");
		  //  // string fd("/data/qw376/top_layer_index_did/");  //on pangolin
		  //  // ofstream index;
		  //  const int buffer_size = 100*1024;
		  //  CluewebReader* Reader  = CluewebFactory();
		  //  for (size_t i=0; i<0; ++i){
		  //  // for (size_t i=0; i<Reader->num_terms; i++){
		  //  // for (size_t i=7500; i<7791; i++){
	 		// 	   const string term = Reader->term_[i];
	 		// 	   // if(term!="009"){
	 		// 	   // 	continue;
	 		// 	   // }
	 		// 	   RawIndexList Rlist = Reader->load_raw_list(Reader->term_[i],i);
	 		// 	   cout<<term<<": "<< i<<endl;
	 		// 	   vector<pinfo> t_list;

	 		// 	   	for(int h = 0; h < Rlist.unpadded_list_length; h++){
	 		// 			pinfo t_p;
	 		// 			t_p.did  = Rlist.doc_ids[h];
	 		// 			t_p.freq = Rlist.freq_s[h];
	 		// 			t_p.score = Rlist.scores[h];
	 		// 			t_list.push_back(t_p);
	 		// 			// cout<<t_p.did<<" "<<t_p.freq<<" "<<t_p.score<<endl;
	 		// 			// int ts;
	 		// 			// cin>>ts;
	 		// 		}
	 		// 		sort(t_list.begin(), t_list.end(), myf);


				// 	/*list less than a block*/
	 				// if(Rlist.unpadded_list_length<=buffer_size){

		 			// 	const int size = Rlist.unpadded_list_length;
		 			// 	unsigned int t_a [size];
		 			// 	float t_f [size];

		 			// 	for(int j = 0; j < Rlist.unpadded_list_length; j++){
		 			// 		t_a[j] = t_list[j].did;
		 			// 		t_f[j] = t_list[j].score;
		 			// 	}

		 			// 	// const string path = fd + Rlist.term;
		 			// 	const string path_d = CONSTS::top_layer_index_did + Rlist.term;
		 			// 	const string path_s = CONSTS::top_layer_index_score + Rlist.term;
		 			// 	FILE *f_d = fopen(path_d.c_str(),"wb");
		 			// 	FILE *f_s = fopen(path_s.c_str(),"wb");
		 			// 	// FILE *f_b = fopen(CONSTS::top_layer_index_did.c_str(),"wb");
		 			// 	// cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), size, f_b)<<endl;
		 			// 	fwrite(t_a, sizeof(unsigned int), size, f_d);
		 			// 	fflush(f_d);
		 			// 	fclose(f_d);
		 			// 	fwrite(t_f, sizeof(float), size, f_s);
		 			// 	fflush(f_s);
		 			// 	fclose(f_s);

	 				// }

				// 	 /*list more than a block*/
	 		// 		if(Rlist.unpadded_list_length>buffer_size){

	 		// 			int block = Rlist.unpadded_list_length / buffer_size;
	 		// 			const int rem = Rlist.unpadded_list_length % buffer_size;
	 		// 			// cout<<"block: "<<block<<" rem: "<<rem<<endl;
	 		// 			// cout<<"buffer_size: "<<buffer_size<<endl;

				// 		const string path_d = CONSTS::top_layer_index_did + Rlist.term;
				// 		const string path_s = CONSTS::top_layer_index_score + Rlist.term;
	 		// 			FILE *f_d = fopen(path_d.c_str(),"wb");
	 		// 			FILE *f_s = fopen(path_s.c_str(),"wb");

	 		// 			for(int k = 0; k < block; k++){
	 		// 			  // cout<<"p0"<<endl;
	 	 // 				  unsigned int t_a [buffer_size];
	 	 // 				  float t_f [buffer_size];
	 	 // 				  // cout<<"p1"<<endl;
	 	 // 				  for(int j = 0; j < buffer_size; j++){
	 	 // 				  	// cout<<"j: "<<j<<endl;
	 		// 				t_a[j] = t_list[j+k*buffer_size].did;
	 		// 				t_f[j] = t_list[j+k*buffer_size].score;
	 		// 			  }
	 	 // 					// cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), buffer_size, f_b)<<endl;
	 		// 			  fwrite(t_a, sizeof(unsigned int), buffer_size, f_d);
	 		// 			  fwrite(t_f, sizeof(float), buffer_size, f_s);
	 		// 			}	

	 		// 			if(rem != 0){
		 	// 				// cout<<"p0"<<endl;
		 	// 				unsigned int t_a [rem];
		 	// 				float t_f [rem];
		 	// 				// cout<<"p1"<<endl;
		 	//  				for(int j = 0; j < rem; j++){
			 // 					t_a[j] = t_list[j+block*buffer_size].did;
			 // 					t_f[j] = t_list[j+block*buffer_size].score;
		 	//  				}
		 	// 				// cout<<"writen bytes: "<<fwrite(t_a, sizeof(unsigned int), rem, f_b)<<endl;
		 	//  				fwrite(t_a, sizeof(unsigned int), rem, f_d);
		 	//  				fwrite(t_f, sizeof(float), rem, f_s);
		 	// 			}
	 		// 			fclose(f_d);
	 		// 			fclose(f_s);

	 		// 	    }
	 		// 		/*---------*/
	 		// 	//   	 ofstream index;
	 		// 	// 	index.open("/home/qi/Dropbox/score_index_m/"+ Rlist.term, ofstream::app);
	 		// 	// 	for(int j = 0; j < Rlist.unpadded_list_length; j++){
				// 	// 	index<<t_list.at(j).did<<endl;
	 		// 	// 			// index<<t_list.at(j).did<<" "<<t_list.at(j).freq<<" "<<t_list.at(j).score<<endl;
				// 	// }
				// 	// index.close();
				//  }

				// delete Reader;

		  // 		return 0;
	  	// }

	   std::string queryFile = (options[QPATH] && options[QPATH].arg) ? std::string(options[QPATH].arg): CONSTS::testingQuery;

	   SingleHashBlocker::expectedBlockSize = (options[OFFSET] && options[OFFSET].arg) ? atoi(options[OFFSET].arg) : 8;
	   //COUT3 << "Expected block size: " << SingleHashBlocker::expectedBlockSize << Log::endl; // togo

	   //togo
	   //std::cout << "########### expected block size: 2^" << SingleHashBlocker::expectedBlockSize << std::endl;

//	   if(nice(-3)) {}
//	   cout<<"dropbox"<<endl;

	   // termsCache Cache;
	   // if(options[ONDEMAND] )
		  //  onDemandCpool::initPool(CONSTS::STORAGE_CAPACITY); //prepare the cache but load none
	   // else {
		  //  // load different terms_mapping
		  //  if (layer == 0)
			 //   Cache.fill_cache(CONSTS::termsMapping.c_str()); //load terms, this one here loading all the terms in the pool directory, ListIterator.cpp 105 load() function
		  //  else
			 //   Cache.fill_cache(CONSTS::termsMapping_Layer.c_str()); //load terms Layering
	   // }

	   termsCache Cache;
	   Cache.fill_cache(); //load terms, this one here loading all the terms in the pool directory, ListIterator.cpp 105 load() function

	   QueryProcessing qp(Cache);
	   topk = 500;
	   CluewebReader* Reader  = CluewebFactory();
	   qp(Reader, queryFile.c_str(), buckets, limit, topk, layer); //process the queries
	   qp.printReport();
	   COUT1 << "start cleanup" << Log::endl;
	   delete Reader;
	   return 0;
}
