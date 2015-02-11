/*
 * qp.cpp
 *
 *  Created on: Nov 24th, 2014
 *      Author: Qi
 */

#include <string.h>
#include <sstream>
#include <math.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <bitset>
#include <dirent.h>

#include "qp.h"
#include "profiling.h" 	// profiling
#include "pfor.h" 		// compression
#include "utils.h"		// helper functions

// (A) Algorithms without Block-Max structures
#include "exhaustiveOR.h"
#include "Wand.h"
#include "Maxscore.h"

// (B) Algorithms with Block-Max structures
// 		(B)-1 Postings-Oriented BM algorithms
#include "PostingOriented_BMW.h"
#include "PostingOriented_BMM.h"
#include "PostingOriented_BMM_NLB.h"

// 		(B)-2 DocID-Oriented BM algorithms
#include "DocidOriented_BMW.h"
#include "DocidOriented_BMM.h" // both BMM and BMM-NLB

// Quantization
#include "DocidOriented_BMM_BMQ.h"
#include "DocidOriented_BMW_BMQ.h"

// Layering
#include "DocidOriented_BMM_Layering.h"

//my algo
#include "intersection.h"
// #include "union_or.h"
#include "algo_toplayer.h"
#include "CluewebReader.h"
#include "And.h"
#include "pairrank.h"
#include "rankinfo.h"


using namespace std;

/* return -1 means finish, 1 means valid */
int readline(FILE* fd, char* line) {
	if(fd==NULL) {
		CERR << "input file is wrong" << Log::endl;
		exit(1);
	}

	char* tt;
	if(!feof(fd)) {
		memset(line,0,5120);
		if(fgets(line,5120,fd) == NULL) {
			//			CERR<<"end of file"<<endl;
			return -1;
		}

		if((tt=strchr(line,'\n')) == NULL) {
			CERR<<"No enter....... \n"<<Log::endl;
			CERR<<"line is "<<line<<Log::endl;
			exit(1);
		}
		*tt = ' ';
	}
	return 1;
}

/* output of query processing: (i) result_log (docID score), (ii) kth score of result set or (iii) create ground truth for results */
class resultsDump {
	FILE *fresult_log;
	FILE *fqscore;
public:
	resultsDump(const char* result_= CONSTS::fResultLog.c_str(), const char* fresult_=CONSTS::fQueryScore.c_str()) {
		fresult_log = fopen(result_,"w");
		fqscore = fopen(fresult_,"w");
	}
	~resultsDump() { fclose(fresult_log); fclose(fqscore); }
	// dump the result for checking
	float operator()(int qn, const std::vector<std::string>& word_l, QpResult* res, int topk) {
		// where the last position that has a valid score and docid is in the topk
		int position_in_topk;
		for (position_in_topk = topk-1; position_in_topk >=0; --position_in_topk)  {
			//std::cout << position_in_topk << "\t"<< res[position_in_topk].score << "\t" << res[position_in_topk].did << std::endl;
			if (( res[position_in_topk].score > -1.0) && (res[position_in_topk].did < CONSTS::MAXD+1))
				break;
		}

//		fprintf(fresult_log, "Query:%d\t%s\t%4.30f\n", qn, word_l[0].c_str(), res[position_in_topk].score);
//		fprintf(fresult_log, "Query: %d\t%4.30f ", qn, res[position_in_topk].score);

//		for(int i = 0; i<word_l.size(); i++){
//			fprintf(fresult_log, "%s ", word_l[i].c_str());
//		}
//		fprintf(fresult_log, "\n");

		/*Not using it for now*/	
		for (int i = 0 ; i <= position_in_topk; i++)
			// fprintf(fresult_log,"%d %4.30f\n", res[i].did, res[i].score);
			fprintf(fresult_log,"%d ", res[i].did );

		fprintf(fresult_log,"\n");
		/*Not using it for now*/	

/*			// create golden standard
			FILE *Golden_threshold_handler = fopen((CONSTS::Golden_Threshold).c_str(),"a");
			for (int i = 0; i<word_l.size(); i++)
				fprintf(Golden_threshold_handler, "%s ", word_l[i].c_str());
			fprintf(Golden_threshold_handler, "%4.30f\n", res[position_in_topk].score);
			fclose(Golden_threshold_handler);
			// end of golden standard
*/
//		fprintf(fqscore,"%4.30f\n", res[position_in_topk].score);

		/*Not using it for now*/	
		// fprintf(fqscore,"Query: %d\t%4.30f ", qn, res[position_in_topk].score);
		// 		for(int i = 0; i<word_l.size(); i++){
		// 			fprintf(fqscore, "%s ", word_l[i].c_str());
		// 		}
		// 		fprintf(fqscore, "\n");
		// return res[position_in_topk].score;
		/*Not using it for now*/
	}

	void skip() {fprintf(fqscore,"0.0\n");}

};

/* profiling and loading lexicon and document lengths */
QueryProcessing::QueryProcessing(termsCache& cache) : qn(0), Cache(cache), trecReader(0) {
	profilerC& p = profilerC::getInstance();
	p.add(" Total and average time overhead of OTF BMG (ms): ",CONSTS::OTFBMG);
	p.add(" Total and average preparation time (ms): ", CONSTS::ESSENTIAL);
	p.add(" Total and average TAAT merging time (ms): ", CONSTS::SKIPS);
	p.add(" Total and average sorting time (ms): ", CONSTS::SORT);
	p.add(" Total and average lookup time (ms): ", CONSTS::SKIPSB);
	p.add(" Total and average Query Processing time (ms): ",CONSTS::ALLQS);
	p.add("candidateS",CONSTS::CAND); // add more values here and the CONSTS ENUM (globals.h) for counters and profiling

	for(size_t i=0; i<CONSTS::NOP; ++i)
		p.initNewCounter();

	for(size_t i=0; i<Cache.size(); ++i)
		termsInCache[std::string(Cache[i].term)]=i;

//	mappingForOnDemand = getTrecWordsMappingMap(CONSTS::zeroMapping);
//	loaddoclen();
	loaddoclen_clueweb();
}

/* Query Processing */
void QueryProcessing::operator()(CluewebReader* Reader, const char* queryLog, const int buckets, const int limit, const int topk, const int layer)	{
	resultsDump resLogger;
	profilerC& p = profilerC::getInstance();

	// load terms -> termID into map for Standard Index
	termsMap lex;
//	fillTermsMap(lex, CONSTS::basic_table);
	fillTermsMap_clueweb(lex, CONSTS::basic_table);

	// termsMap pair_lex;
	// fillTermsMap_pairs(pair_lex);  //pair indexes needed when do qp with pairs

	/*11.8*/
//	ulong available_budget = ComputeSpaceBudget(CONSTS::percent_of_index);
	// termsMap pair_lex = load_termpair_cutoffs(CONSTS::cutoffs_path + CONSTS::term_pair_cutoffs+CONSTS::percent_of_index_str+"/Prob_Learned");
//	termsMap pair_lex = load_termpair_cutoffs("/home/qw376/pair_lex/trec_06_testing_lex_sorted");
	//termsMap pair_lex = load_termpair_cutoffs(CONSTS::term_pair_cutoffs+CONSTS::percent_of_index_str+"/Prob_Learned_Small");
	termsMap pair_lex = load_termpair_cutoffs("/home/qw376/pair_lex/mix_pair_cutoffs");
	termsMap opair_lex = load_termpair_cutoffs("/home/qw376/pair_lex/trec_06_testing_lex_sorted");
//	termsMap pair_lex;
	//termsMap toplayer_lex = load_toplayer_cutoffs(CONSTS::cutoffs_path + CONSTS::top_layer_cutoffs+CONSTS::percent_of_index_str+"/original_ind"); 
	//termsMap toplayer_lex = load_toplayer_cutoffs(CONSTS::top_layer_cutoffs+CONSTS::percent_of_index_str+"/Prob_Learned_Small");
//	termsMap toplayer_lex = load_toplayer_cutoffs("/home/qw376/pair_lex/mix_term_cutoffs_replace20k");
	termsMap toplayer_lex = load_toplayer_cutoffs("/home/qw376/pair_lex/mix_term_cutoffs_noreplace20k");
	termsMap otoplayer_lex = load_toplayer_cutoffs("/home/qw376/cutoffs/unigram/trec06/100/original_ind");
	bool is_top_layer = true;

  	//Load top-layer model
 	vector<vector<float> > topLayerModel = LoadModel(is_top_layer);

 	//Load term-pair layer model
 	is_top_layer = false;
 	vector<vector<float> > termPairLayerModel = LoadModel(is_top_layer);

 	//Load the blocks
	vector<uint> posting_block_boundary = LoadBlock(0);  
	vector<uint> posting_block_sizes = LoadBlock(1);  
	vector<uint> listlen_block_boundary = LoadBlock(2);  
	vector<uint> intesection_block_boundary = LoadBlock(3); 

	cout<<"percent_of_index: "<<CONSTS::percent_of_index_str<<endl;
	cout<<"Query_Budget_TopLayer: "<<CONSTS::Query_Budget_TopLayer<<endl;
	cout<<"Query_Budget_TermPair: "<<CONSTS::Query_Budget_TermPair<<endl;
	cout<<"Num_Doc_for_Lookups: "<<CONSTS::Num_Doc_for_Lookups<<endl;
	cout<<"lookup_budget: "<<CONSTS::lookup_budget<<endl;
	cout<<"num_of_candidate: "<<CONSTS::num_of_candidate<<endl;


	// Auxiliary functions
	//List_Length_Distribution(List_Lengths);
	//Bucketized_Maxscore_Stats(); // Usage: load file with <maxscore, unpadded list length> print out distribution
	//Query_log_Maxscore_Length_Correlation(lex);
	//Get_Average_Document_Length(); // auxiliary testing function

	// Standard Index
	// QueryLogManager logManager(queryLog,&lex);
	int version = 1;
	QueryLogManager logManager(queryLog,&lex,version);

	COUT1 << logManager.size() << " queries are available" << Log::endl;

	QueryLogManager::queriesFileIterator queriesFIterator(logManager);
	queriesFIterator.changeSkipPolicy(limit,buckets);

	cout << "################################################################" << endl;
	cout << "QP configuration\ntopk:" << topk << "\tbuckets: " << buckets << "\tlimit:" << limit << "\tlayer: " << layer << endl;
	cout << "################################################################" << endl;

	/* perform query processing for each query */
	while( qn < 27) { 

		++queriesFIterator;
		if(queriesFIterator == logManager.end()) {
//			COUT << "finished benchmark" << Log::endl;
			break;
		}

		std::vector<std::string> word_l = (*queriesFIterator);
		int qid = queriesFIterator.Getqid();

		while(word_l.size()>=10) {	// no queries more than 10 - in our query log of 1k there is just 1 query of this size
			++queriesFIterator;
			CERR << "skip " << word_l << Log::endl;
			word_l = (*queriesFIterator);
			qid = queriesFIterator.Getqid();
		}

		/*Control query size*/
//		if(word_l.size()<5)
//		 	continue;
		qn++;
 		//if(qn <= 0)
                //  continue;
		// ########################################################################################
		// Initializing (opening) list (lps=list pointer) structures based on standard index
		lptrArray lps = openLists(word_l, lex); // Standard Index, need to be opened anytime

		pairlists pls = openPLists(word_l, pair_lex);//pair indexes needed when do qp with pairs

		pairlists opls = openPLists(word_l, opair_lex);

		toplayers tls = openTopLayers(word_l, toplayer_lex); //toplayers

		toplayers otls = openTopLayers(word_l, otoplayer_lex);



		// NOTE: Uncomment one of the following algorithms and make sure you call it with the appropriate arguments and returns the appropriate value
		// ########################################################################################
		// (1) Algorithms without Block-Max Indexes
		// ExhaustiveOR wand(pages);                  // Exhaustive OR

		// Wand wand(pages);                    		 // WAND
		// Maxscore wand(pages);              	 	 // Maxscore
		// pairalgo wand(pages);                	//to generate single term binary files(debug purpose)
		// rankinfo wand(pages);				//Cal depth info for layers 
		// pairrank wand(pages);				//Cal depth info for pairs
		// intersection wand(pages);					//Generate pair index
		// union_or wand(pages);					//get the union
		// And wand(pages);					//Generate pair index
		algo_toplayer wand(pages);				//test the speed


		// ########################################################################################
		// (D) Run algorithms (except OPT algorithms)
		// QpResult res[topk]; 			// result set - comment line if PriorityArray is used to maintain the results !!!!!!****uncomment if want to use res output

		// cout<<"the query: ";
		// for(int h = 0; h < word_l.size(); h++){
		// 	cout<<word_l[h]<<" ";
		// }
		// cout<<endl;

		// CluewebReader* Reader  = CluewebFactory(); //for new algo
		// p.start(CONSTS::ALLQS); 		// Start measuring qp time - NOTE: that OTF BMG is already measured if DocID-Oriented Block-Max structures are used (DOCIDBLOCKMAX is defined)

		// various default parameters for running algorithms
		// wand(lps, topk, res);   //or, pair
		// wand(qn, lps);   //rankinfo
		// wand(qn, pls);   //pairrank
		// wand(Reader, lps);   //intersection, And, union
		// wand(Reader, word_l);   //intersection, And, union
		wand(topLayerModel, termPairLayerModel, posting_block_boundary, posting_block_sizes, listlen_block_boundary, intesection_block_boundary, Reader, qn, tls, otls, pls, opls,lps, topk, p, qid);   //algo_toplayer
		// wand(lps, topk, 0.0f);
		// wand(lps, topk, res, 0.0f);  //Maxscore
		// wand(lps, topk, res, 0.0f, 0);
		// PriorityArray<QpResult> resultsHeap = wand(lps, topk);
		// p.end(CONSTS::ALLQS); 			// Stop measuring qp time - NOTE: that OTF BMG is already measured if DocID-Oriented Block-Max structures are used (DOCIDBLOCKMAX is defined)
		// delete Reader;	//for new algo
		// ########################################################################################
		// obtain result set after query processing and perform sanity check with ground truth
		// Note: depending on the algorithm the results are either returned to (a) QpResult or (b) PriorityArray<QpResult>, so comment and uncomment the next lines accordingly
		// float score = resLogger(qn, word_l, res, topk);		// (a) QpResult Qi's commenting needed!!!!   !!!!!!****uncomment if want to use res output
		// resultsHeap.sortData();							// (b) PriorityArray<QpResult>
		// float score = resultsHeap.getV()[topk-1].score;	//     PriorityArray<QpResult>

		// ########################################################################################
	} // end of Query Processing for the current Query
}

/* pair cutoffs 11.8*/
/*
termsMap QueryProcessing::load_termpair_cutoffs(const string input){

	termsMap lex;
	FILE* flex = fopen(input.c_str(),"r");
	if(flex==NULL)
		cout << input << " could not be opened" << endl;

	char term1[1000];
	char term2[1000];
	string pair;
	int length;
	while( fscanf(flex,"%s %s %u\n", term1, term2, &length)!= EOF ){
		string str1(term1);
		string str2(term2);
		pair = str1 + "+" +str2;
		lex[pair] = length;
	}
	cout<<"loading term pair cutoffs of size: "<<lex.size()<<" from file: "<<input<<endl;
	fclose(flex);
	assert(lex.size());
	return lex;
}
*/

termsMap QueryProcessing::load_termpair_cutoffs(const string input){

        termsMap lex;
        FILE* flex = fopen(input.c_str(),"r");
        if(flex==NULL)
                cout << input << " could not be opened" << endl;

        char pair[1000];
        int length;
        while( fscanf(flex,"%s %u\n", pair, &length)!= EOF ){
                string str(pair);
                lex[str] = length;
//		lex[str] = 300;
        }
        cout<<"loading term pair cutoffs of size: "<<lex.size()<<" from file: "<<input<<endl;
        fclose(flex);
        assert(lex.size());
        return lex;
}
/* toplayer cutoffs 11.8*/
termsMap QueryProcessing::load_toplayer_cutoffs(const string input){

	termsMap lex;
	FILE* flex = fopen(input.c_str(),"r");
	if(flex==NULL)
		cout << input << " could not be opened" << endl;

	char term[1000];
	int length;
	while( fscanf(flex,"%s %u\n", term, &length)!= EOF )
		lex[string(term)] = length;
	cout<<"loading toplayer cutoffs of size: "<<lex.size()<<" from file: "<<input<<endl;
	fclose(flex);
	assert(lex.size());
	return lex;
}

/* load pairs -> listlen into map from file */
// void QueryProcessing::fillTermsMap_pairs(termsMap& lex){
// 	DIR *dir;
// 	struct dirent *ent;

// 	// if( (dir = opendir("/data/qw376/pair_index/") )!= NULL) {
// 	if( (dir = opendir(CONSTS::pair_index.c_str()))!= NULL) {
// 		 /* print all the files and directories within directory */
//   		while ((ent = readdir (dir)) != NULL) {
//    	 	  // printf ("%s\n", ent->d_name);
//   		  string filename = string(ent->d_name);
//   		  if((filename.compare(".")!=0)&&(filename.compare("..")!=0)){
//    	 	  		lex[filename] = 500; //pair depth
//    	 	}
//   		}
//   		cout<<"pairlex size: "<<lex.size()<<endl;
//   		closedir (dir);
// 		} else {
//  		 /* could not open directory */
//   		cout<<"open dir failed"<<endl;
// 	}
// }

/* load pairs -> listlen into map from file */
void QueryProcessing::fillTermsMap_pairs(termsMap& lex){

	FILE* flex = fopen(CONSTS::pair_lex.c_str(),"r");
	if(flex==NULL)
		CERR << CONSTS::pair_lex << " could not be opened" << EFATAL;

	char term[2000];
	int length;
	while( fscanf(flex,"%s %d\n",term, &length)!= EOF )
	//	lex[string(term)] = length;
	//	lex[string(term)] = 300;
	cout<<"pairlex size: "<<lex.size()<<endl;
	fclose(flex);
	assert(lex.size());

}

/* load terms -> listlen into map from file */
void QueryProcessing::fillTermsMap_clueweb(termsMap& lex, const std::string& path) {

	ifstream doc_lexikon_stream;
	doc_lexikon_stream.open(path.c_str());

	string lexicon_line;
	string term;
	string termid_s;
	string listLen_s;
	int listLen;
	int termid;

	while (getline(doc_lexikon_stream, lexicon_line)){

	//			cout<<"linenum: "<<linenum<<endl;
	//			if(linenum == 10)
	//				break;

				//term, first field
				string::iterator itr1 = lexicon_line.begin();
				string::iterator start1 = itr1;
				while(itr1 != lexicon_line.end() && !isspace(*itr1)){
					++itr1;
				}
				term = string(start1, itr1);
				// cout<<"term: "<<term<<" ";
				//term

				 //termID, come after the 1st space
			  	 start1 = itr1+1;
			  	 itr1++;
			  	 while (!isspace(*itr1)) {
			  	 ++itr1;
			  	 }

			  	 termid_s = string(start1,itr1);
			  	 termid = atoi(termid_s.c_str());
			  	 // cout<<"termid: "<<termid<<" ";
			  	 //termID

			  	 //listLen, come after the 2nd space
			  	 start1 = itr1+1;
			  	 itr1++;
			  	 while (!isspace(*itr1)) {
			  	 ++itr1;
			  	 }
			  	 listLen_s = string(start1, itr1);
			  	 listLen = atoi(listLen_s.c_str());
				 // cout<<"listLen: "<<listLen<<endl;

			  	 //listLen

			  	lex[string(term)] = listLen;
	}

	cout<<"termlex size: "<<lex.size()<<endl;

	doc_lexikon_stream.close();

}

/* load terms -> listlen into map from file */
void QueryProcessing::fillTermsMap(termsMap& lex, const std::string& path) {
	FILE* flex = fopen(path.c_str(),"r");
	if(flex==NULL)
		CERR << path << " could not be opened" << EFATAL;

	char term[1024];
	int length;
	while( fscanf(flex,"%s\t%d\n",term, &length)!= EOF )
		lex[string(term)] = length;

	fclose(flex);
	assert(lex.size());
}

/* load document length from file */
unsigned int* QueryProcessing::loaddoclen(const char* fname)	{
	int docn = CONSTS::MAXD;
	pages = new unsigned int[ docn + 128];

	FILE *fdoclength = fopen64(fname,"r");
	if( fdoclength == NULL)
		CERR <<" doc length file is missing "<< EFATAL;

	if( fread(pages,sizeof(int), docn, fdoclength) != docn )
		CERR<<"wrong doc len "<< EFATAL;

	fclose(fdoclength);

	for(int i =0;i<128 ; i++)
		pages[docn + i] = docn/2;

	return pages;
}

/* load document length from file for clueweb*/
unsigned int* QueryProcessing::loaddoclen_clueweb(const char* fname)	{

	int docn = CONSTS::MAXD;
	pages = new unsigned int[ docn + 128];

	ifstream inputstream;
	inputstream.open(fname);
	string curr_line_str;
	string doc_id_s;
	string doc_length_s;
	int doc_id;
	int doc_length;

	while (getline(inputstream, curr_line_str)) {
		string::iterator itr = curr_line_str.begin();
		string::iterator start = itr;

		//the first is doc id
		while(itr != curr_line_str.end() && !isspace(*itr)){
			++itr;
		}

		doc_id_s = string(start,itr);
//			  cout<<"doc_id: "<<doc_id_s<<endl;
		doc_id = atoi(doc_id_s.c_str());

		start = itr + 1;
		itr++;

		//the second is doc length
	    while (itr != curr_line_str.end()) {
			++itr;
	    }
	    doc_length_s = string(start,itr);
//			  cout<<"doc_length: "<<doc_length_s<<endl;
	    doc_length = atoi(doc_length_s.c_str());

	    pages[doc_id] = doc_length;
	}

	inputstream.close();

	for(int i =0;i<128 ; i++)
		pages[docn + i] = docn/2;

	return pages;
}

/* Print reports
Usage: 1. add to globals.h, CONSTS namespace in the counter enumerator for more counters
       2. add code for counters in the appropriate place and
	   3. use COUT3 and the name of the counter to print the results */
void QueryProcessing::printReport() const {
	profilerC& p = profilerC::getInstance();
	p.printReport();

	std::locale loc("");
	std::cout.imbue(loc);

	//COUT3 << "TOTAL EVALUATION: " << p.getCounter(CONSTS::EVAL) << Log::endl;
	//COUT3 << "Total nextGEQ: " << p.getCounter(CONSTS::NEXTGEQ) << Log::endl;
}

//shouldn't use this now
void QueryProcessing::onDemandCall(const std::string& term) {
	// COUT1 << "on demand: " << term << Log::endl;
	// stringIntVectorsPair tmap;
	// tmap.first.push_back(term);
	// tmap.second.push_back(mappingForOnDemand[term]);
	// if (! trecReader)
	// 	trecReader = TrecFactory(CONSTS::trecRawRoot);
	// TrecFactory(*trecReader,CONSTS::trecRoot,0,10,tmap);
	// Cache.addSingleToCache(term);
	// termsInCache[term]=Cache.size()-1;
}

/* core openList function */
lptrArray QueryProcessing::openLists(const std::vector<std::string>& word_l, termsMap& lex){
	// create structure and reserve resources
	lptrArray lps;
	lps.reserve(word_l.size());
	//prepare peers
	std::vector<std::string> peers( word_l.size());

	for(size_t i=0; i<word_l.size(); i++) {
		if(termsInCache.find(word_l[i]) == termsInCache.end()) //no term! try on-demand processing from rawindex; The load() func in the fill_term() func put everything in the cahce(centry) ListIterator.cpp 130
			// onDemandCall(word_l[i]);
			cout<<"onDemandCall problem, you probably didn't build the baby index or gave the wrong dir.."<<endl;
		peers.push_back(Cache[ termsInCache[word_l[i]]].cpool.getName() );
	}

	// note: this one breaks when query has duplicate terms
	for(size_t i=0; i<word_l.size(); i++)	{
		lptr* lpsi = &(Cache[termsInCache[word_l[i]]]);
		lpsi->open(peers,lex.find((word_l[i]))->second);
		lps.push_back(lpsi);
	}
	return lps;
}


/* core openList function, just for the intersection computation 12.3*/
// lptrArray QueryProcessing::openLists(const std::vector<std::string>& word_l, termsMap& lex){
// 	// create structure and reserve resources
// 	lptrArray lps;
// 	lps.reserve(word_l.size());
// 	for(size_t i=0; i<word_l.size(); i++) {	
// 		lptr c;
// 		c.term = word_l[i];
// 		cout<<c.term<<endl;
// 		lptr* lpsi = &c;
// 		lps.push_back(lpsi);
// 	}
// 	return lps;
// }

/*open pairs list, if pairs in the dir is not sorted alphabatically*/
// pairlists QueryProcessing::openPLists(const std::vector<std::string>& word_l, termsMap& lex){

// 	// cout<<"-------"<<endl;
// 	pairlists pls;
// 	string pair;

// 	for(int i = 0; i<word_l.size(); i++){
// 	  		for(int j = i+1; j< word_l.size(); j++){
// 	  			// cout<<word_l.at(i)<<" "<<word_l.at(j)<<endl;
// 	  	// 		if(word_l.at(i).compare(word_l.at(j))<0)
// 				//     pair = word_l[i] + "+" + word_l[j];
// 				// else
// 	  	// 		    pair = word_l[j] + "+" + word_l[i];
// 	  			// cout<<"generated: "<<pair<<endl;
// 	  			pair = word_l[i] + "+" + word_l[j];
// 	  			termsMap::iterator itr = lex.find(pair);
// 	  			if( itr!= lex.end() ) {			
// 					pls.pairnames.push_back(pair);
// 					pls.lengths.push_back(itr->second);
// 	  		 	}
// 	  		 	else{
// 	  		 		pair = word_l[j] + "+" + word_l[i];
// 	  		 		if( itr!= lex.end() ) {			
// 					pls.pairnames.push_back(pair);
// 					pls.lengths.push_back(itr->second);
// 	  		 		}
// 	  		 	}
// 	  		 	// cout<<"selected: "<<endl;
// 	  		 	// for(int i = 0; i<pls.lengths.size(); i++){
// 	  		 	// 	cout<<pls.pairnames[i]<<" "<<pls.lengths[i]<<endl;
// 	  		 	// }
// 	  		}
// 	 }

// 	return pls;
// }

/*open pairs list, if pairs in the dir is sorted alphabatically*/
pairlists QueryProcessing::openPLists(const std::vector<std::string>& word_l, termsMap& lex){

	// cout<<"-------"<<endl;
	pairlists pls;
	string pair;

	for(int i = 0; i<word_l.size(); i++){
	  		for(int j = i+1; j< word_l.size(); j++){
	  			// cout<<word_l.at(i)<<" "<<word_l.at(j)<<endl;
	  			if(word_l[i].compare(word_l[j])<0)
				    pair = word_l[i] + "+" + word_l[j];
				else
	  			    pair = word_l[j] + "+" + word_l[i];

	  			// cout<<"generated: "<<pair<<endl;
	  			termsMap::iterator itr = lex.find(pair);
	  			if( itr!= lex.end() ) {			
					pls.pairnames.push_back(pair);
					pls.cutoffs.push_back(itr->second);
	  		 	}
	  		 	// cout<<"selected: "<<endl;
	  		 	// for(int i = 0; i<pls.lengths.size(); i++){
	  		 	// 	cout<<pls.pairnames[i]<<" "<<pls.lengths[i]<<endl;
	  		 	// }
	  		}
	 }

	return pls;
}

/*open top layers 11.8*/
toplayers QueryProcessing::openTopLayers(const std::vector<std::string>& word_l, termsMap& lex){

	toplayers tls;
	string term;

	for(int i = 0; i<word_l.size(); i++){
  		termsMap::iterator itr = lex.find(word_l[i]);
  			if( itr!= lex.end() ) {			
				tls.terms.push_back(word_l[i]);
				tls.cutoffs.push_back(itr->second);
  		 	}
	 }
	return tls;
}

/* Load Queries - Version for Standard Index (Layering) */
QueryLogManager::QueryLogManager(const char* fname, termsMap *l) :lex(l){
	assert(l->size());

	FILE *fq = fopen(fname,"r");
	queriesD.reserve(1000);
	std::vector<std::string> terms;
	while(1)  {
		terms.clear();
		char line[5120];

		if( readline(fq, line) == -1 ) //this is the end
			break;

		char *word;
		word = strtok(line," ");

		if( lex->find(std::string(word))!= lex->end() ) //only if term is in cache...
			terms.push_back(word);
		else
			CERR << word << " is not in lex" << Log::endl;

		while((word = strtok(NULL," ")) != NULL){
			if( lex->find(std::string(word))!= lex->end() ) //only if term is in cache...
				terms.push_back(word);
			else
				CERR << word << " is not in lex" << Log::endl;
		}

		// if(terms.size()==2) //for pairs only added by Qi, for pairinfo only, for others need to comment this line
		queriesD.push_back(terms);
		setScoreForQuery(terms,0.0);
	}
	//COUT << queriesD.size() << " queries loaded" << Log::endl;
	fclose(fq);

}

/* Load Queries - Version for Standard Index (Layering) */
QueryLogManager::QueryLogManager(const char* fname, termsMap *l, int& version) :lex(l){
	assert(l->size());

	FILE *fq = fopen(fname,"r");
	queriesD.reserve(1000);
	std::vector<std::string> terms;
	while(1)  {
		terms.clear();
		char line[5120];

		if( readline(fq, line) == -1 ) //this is the end
			break;

		//the first one is qid
		char *word;
		word = strtok(line," ");
		qidD.push_back(atoi(word));

		//the remaining parts are the query terms
		while((word = strtok(NULL," ")) != NULL){
			if( lex->find(std::string(word))!= lex->end() ) //only if term is in cache...
				terms.push_back(word);
			else
				CERR << word << " is not in lex" << Log::endl;
		}

		// if(terms.size()==2) //for pairs only added by Qi, for pairinfo only, for others need to comment this line
		queriesD.push_back(terms);
		setScoreForQuery(terms,0.0);
	}
	//COUT << queriesD.size() << " queries loaded" << Log::endl;
	fclose(fq);

}

/* Load Queries - Version for Impact-Sorted Index (Layering) */
QueryLogManager::QueryLogManager(const char* fname, termsMap *l, bool& nothing) :lex(l){
	assert(l->size());

	FILE *fq = fopen(fname,"r");
	queriesD.reserve(1000);
	std::vector<std::string> terms;

	while(1)  {
		terms.clear();
		char line[5120];

		if( readline(fq, line) == -1 )
			break;

		char *word;
		word = strtok(line," ");

		// First query term in line
		if ( lex->find(std::string(word))!= lex->end() ) // if term exists in lex (with no layer), add
			terms.push_back(std::string(word));
		else {
			std::string check_term (std::string(word)+CONSTS::GOOD_TERM_SUFFIX);
			if ( lex->find(check_term)!= lex->end() ) // if term exists in lex, but with layers, add
				terms.push_back(std::string(word));
			else	// if term does not exist, print msg
				CERR << std::string(word) << " was not found in lexicon" << Log::endl;
		}

		// Remaining query terms in the same line
		while((word = strtok(NULL," ")) != NULL){
			if ( lex->find(std::string(word))!= lex->end() ) // if term exists in lex (with no layer), add
				terms.push_back(std::string(word));
			else {
				std::string check_term (std::string(word)+CONSTS::GOOD_TERM_SUFFIX);
				if ( lex->find(check_term)!= lex->end() ) // if term exists in lex, but with layers, add
					terms.push_back(std::string(word));
				else	// if term does not exist, print msg
					CERR << std::string(word) << " was not found in lexicon" << Log::endl;
			}
		}

		queriesD.push_back(terms);
		setScoreForQuery(terms,0.0);
	}

	COUT2 << queriesD.size() << " queries loaded" << Log::endl;
	fclose(fq);

}

// Given a vector<uint> print out for sanity checking.
void PrintBlockVec(const vector<uint> vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << vec[i] << " ";
  }
  cout << endl;
}

// Given the type of block to load, load from the appropriate file the block.
vector<uint> LoadBlock(const uint type) {
  string input;
  if (type == 0) {  // posting block boundary
    input = BLOCKS::posting_block_boundary;
  } else if (type == 1) {  // posting block sizes
    input = BLOCKS::posting_block_sizes;
  } else if (type == 2) {  // listlen block boundary
    input = BLOCKS::listlen_block_boundary;
  } else if (type == 3) {  // intersection block boundary
    input = BLOCKS::intersection_block_boundary;
  } else {
    cout << "No support for this type of blocks. Exiting ..." << endl;
  } 

  vector<uint> vec; // changed from ulong (64x machines)
  ifstream fH(input.c_str());

  if (fH.good()) {
    while (fH.good()) {
      string line;
      getline(fH, line);
      vec = TokenizeStringAndConvertToUint(line);
      break;  // since we only have one line.
    }
  }

  // Reporting.
  cout << "Block was loaded from file: " << input << endl; 

  return vec;
}

// Given a boolean declaring whether this is the top_layer or the term_pair,
// load from the corresponding file (from MODEL) to a vec of vec of floats.
vector<vector<float> > LoadModel(bool is_top_layer) {
  vector<vector<float> > model;
  string input;
  if (is_top_layer) {
    input = MODEL::top_layer_model_scores;
  } else {
    input = MODEL::term_pair_model_scores;
  }
  ifstream fH(input.c_str());

  if (fH.good()) {
    while (fH.good()) {
      string line; 
      getline(fH, line);
      vector<float> scoreVec = TokenizeStringAndConvertToFloat(line);
      model.push_back(scoreVec);
    }
  }
  fH.close();

  // Reporting.
  cout << "Loading model from file: " << input << endl;

  return model;
}

// Similar to LoadModel, but we return a single vector and the number of columns.
vector<float> LoadModelVec(bool is_top_layer, uint &num_cols) {
  vector<float> model;
  uint cnt = 0;
  string input;
  if (is_top_layer) {
    input = MODEL::top_layer_model_scores;
  } else {
    input = MODEL::term_pair_model_scores;
  }
  ifstream fH(input.c_str());

  if (fH.good()) {
    while (fH.good()) {
      string line; 
      getline(fH, line);
      vector<float> scoreVec = TokenizeStringAndConvertToFloat(line);
      for (size_t i = 0; i < scoreVec.size(); ++i) {
        model.push_back(scoreVec[i]);
      }
      // Set number of columns assuming all lines have the same size.
      if (cnt == 0) {
        num_cols = scoreVec.size(); 
      }
      ++cnt;
    }
  }
  fH.close();

  // Reporting.
  cout << "Loading model from file: " << input << endl;

  return model;
}

// Given a string, tokenize it and convert values to float. Return vector of floats.
vector<float> TokenizeStringAndConvertToFloat(string line) {
  vector<float> tokenizedFloats;
  istringstream iss(line);
  do {
    string token;
    iss >> token;
    if (token=="") continue;
    tokenizedFloats.push_back(convertStrToFloat(token));
  } while (iss);
  return tokenizedFloats;
}

// Given a string, tokenize it and convert values to uint. Return vector of uints.
vector<uint> TokenizeStringAndConvertToUint(string line) {
  vector<uint> tokenizedUints;
  istringstream iss(line);
  do {
    string token;
    iss >> token;
    if (token=="") continue;
    tokenizedUints.push_back(convertStrToUint(token));
  } while (iss);
  return tokenizedUints;
}

// Given a string, convert it to float and return it.
// TODO(dimopoulos): templatize it ?
float convertStrToFloat(const string str) {
  istringstream buffer(str);
  float tmp;
  buffer >> tmp;
  return tmp;
}

// Given a string, convert it to uint and return it.
uint convertStrToUint(const string str) {
  istringstream buffer(str);
  uint tmp;
  buffer >> tmp;
  return tmp;
}

/*11.8*/
ulong QueryProcessing::ComputeSpaceBudget(const double percent_of_index) {
  assert(percent_of_index > 0 && percent_of_index < 1);  // range of percent (0,1)
  ulong space_budget = percent_of_index * (double) CONSTS::kClueTotalPostings;
  cout << "The " << percent_of_index << " percent of index corresponds to "
   	<< space_budget << " postings." << endl;
  return space_budget;
}

/* load terms -> termID into map from file for LAYERING */
void QueryProcessing::Term_To_list_length_Map(termsMap& lex, const std::string& path, const bool& unpadded_list_length) {
	FILE* flex = fopen(path.c_str(),"r");
	if (flex==NULL)
		CERR << path << " could not be opened" << EFATAL;

	char term[1024];
	int length;
	int unpadded_length;
	while( fscanf(flex,"%s\t%d\t%d\n",term, &length, &unpadded_length)!= EOF ) {
		if (unpadded_list_length)
			lex[string(term)] = unpadded_length;
		else
			lex[string(term)] = length;
	}

	fclose(flex);
	assert(lex.size());
}

/* Write lists to file for quick experiments */
void QueryProcessing::Write_Inverted_Lists_in_File(std::vector<std::string>& word_l, lptrArray& lps) {
	// check for missing lists and do not store anything
	if (lps.size()==word_l.size()) {
		int list = 0;

		// folder to store inverted lists info
		std::string folder = "list/";
		std::string list_length_suffix = "_len";

		// file name to store query log
		std::string query_log = folder + "Query_log";
		// create file handler for query log, open it and write Query
		std::ofstream query_log_handler;
		query_log_handler.open (query_log.c_str(), ios::in | ios::out | ios::app);
		query_log_handler << word_l.size() << " ";

		// write query
		for (int i=0; i<word_l.size(); i++)
			query_log_handler << word_l[i] << " ";
		query_log_handler << "\n";
		query_log_handler.close();

		// Obtain scores, docids, freq for each list and dump them to file
		for(auto it = lps.cbegin(); it!=lps.cend(); ++it) {
			RawIndexList rilist = trecReader->getRawList((*it)->term,mappingForOnDemand[(*it)->term]);

			// file name based on term
			std::string filename = word_l[list].c_str();

			// construct full path
			std::string filename_path = folder + filename;

			// create file handler and open it
			std::ofstream file_handler;
			file_handler.open (filename_path.c_str(), ios::in | ios::out | ios::trunc);

			int list_length = 0;
			// for all docids print out to file the docid and its BM25 score (or frequency optional)
			for (int i=0; i<rilist.doc_ids.size(); i++) {
				// print only the ones inside the docid space
				if (rilist.doc_ids.at(i)<CONSTS::MAXD) {
					list_length++;
					// human readable variables
					int docid = rilist.doc_ids.at(i);
					float score = rilist.scores.at(i);
					//float frequency = rilist.freq_s.at(i);
					//float final_score = lps[list]->calcScore(frequency, pages[rilist.doc_ids.at(i)]);

					// Print out to file either of these versions
					file_handler << docid << " " << score << "\n";
					// full version with docid, score, freq, computed score
					//file_handler << docid << " " << score << " " << frequency << " " << final_score << "\n";
				}
			}
			// close handler
			file_handler.close();

			// write the list length in a separate file
			std::ofstream file_list_length_handler;
			std::string filename_list_length_path = folder + filename + list_length_suffix;
			file_list_length_handler.open (filename_list_length_path.c_str(), ios::in | ios::out | ios::trunc);
			file_list_length_handler << list_length << "\n";
			file_list_length_handler.close();

			// increase list counter
			++list;
		}
	}
}

/* Compare the Space requirements of Postings-Oriented Block-Max Indexes and DocID-Oriented ones */
void QueryProcessing::Space_Comparison(std::vector<int>& List_Lengths) {
	long long BMW_space = 0;
	unsigned long long Hash_space = 0;
	unsigned long long Hash_oracle_space = 0;
	long BMW_block_max = 64; // 1 float value for block max per 64 postings
	int Bucket_size = 19;
	long long GB_convertion = 1024*1024*1024;
	std::vector<unsigned long long> Bucket_Hash_Values (Bucket_size, 0);
	std::vector<int> Bucket (Bucket_size, 0);

	// parameters
	Bucket.at(0) = 10; // 2^7   -- 14
	Bucket.at(1) = 10; // 2^8   -- 14
	Bucket.at(2) = 10; // 2^9   -- 13
	Bucket.at(3) = 10; // 2^10  -- 12
	Bucket.at(4) = 6; // 2^11 --- 12
	Bucket.at(5) = 6; // 2^12 -- 10
	Bucket.at(6) = 7; // 2^13 -
	Bucket.at(7) = 7; // 2^14 --- 9
	Bucket.at(8) = 7; // 2^15 -- 8
	Bucket.at(9) = 8; // 2^16 - -----8
	Bucket.at(10) = 8; // 2^17 --- 8
	Bucket.at(11) = 7; // 2^18
	Bucket.at(12) = 6; // 2^19 -- was 7
	Bucket.at(13) = 6; // 2^20 -- was 7
	Bucket.at(14) = 6; // 2^21 -- was 7
	Bucket.at(15) = 6; // 2^22 -- was 7
	Bucket.at(16) = 6; // 2^23 -- was 7
	Bucket.at(17) = 6; // 2^24 -- was 7

	std::cout << "########## Parameter List #############" << std::endl;
	// print out parameter list
	for (int i=0; i<Bucket.size(); i++)
		std::cout << Bucket.at(i) << std::endl;
	std::cout << "#################" << std::endl;

	// for each list loaded
	for (int i=0; i<List_Lengths.size(); i++) {
		// compute BMW space (both BMW_space are working but the second is more elegant)
		//BMW_space += (List_Lengths.at(i)%BMW_block_max == 0) ? (int) List_Lengths.at(i)/BMW_block_max : (int) List_Lengths.at(i)/BMW_block_max + 1;
		BMW_space += 1 + (( List_Lengths.at(i) - 1 ) / BMW_block_max );

		// compute Hash space
		unsigned int lenBits = intlog2(List_Lengths.at(i));
		int Hash_block_max = lenBits < 8 ? 10 : 6;
		float Hash_blocksize = pow(2, (float) Hash_block_max);
		float Total = pow(2, (float) 25);
		Hash_space += 1 +  (( Total - 1 ) /  Hash_blocksize );

/*		// Expected postings per block space computation
		unsigned int expectedBits = 6;
		unsigned int maxdBits = CONSTS::MAXDBITS;
		unsigned int bbits = (maxdBits-lenBits+expectedBits);
		float blocksize = pow(2, (float)bbits);
		Hash_oracle_space += 1 + (( CONSTS::MAXD - 1 ) / blocksize );
*/
		unsigned int bbits;
		int B_counter = 0;
		for (int i=7; i<=25; i++) {
			if (lenBits<i) {
				bbits = Bucket.at(B_counter);
				float blocksize = pow(2, (float) bbits);
				Hash_oracle_space += 1 + (( CONSTS::MAXD - 1 ) / blocksize );
				Bucket_Hash_Values.at(B_counter) += 1 + (( CONSTS::MAXD - 1 ) / blocksize ); // to change maxd ?
				break;
			} else
				++B_counter;
		}
	}

	long long BMW_size = BMW_space*sizeof(float);
	long long Hash_size = Hash_space*sizeof(float);
	unsigned long long Hash_size_oracle = Hash_oracle_space*sizeof(float);

	double BMW_GB = (double) BMW_size/GB_convertion;
	double Hash_GB = (double) Hash_size/GB_convertion;
	double Hash_oracle_GB = (double) Hash_size_oracle/GB_convertion;

	// Print Results
	std::cout << "-----------------------------------------------------" << std::endl;
	std::cout << "Total number of distinct terms: " << List_Lengths.size() << std::endl;
	std::cout << "BMW space: " << BMW_space << " block maxes values" << std::endl;
	std::cout << BMW_size << " bytes total space and " << BMW_GB << " GB of space" << std::endl;
	std::cout << "-----------------------------------------------------" << std::endl;
	std::cout << "Hash space (2^6): " << Hash_space << " block maxes values" << std::endl;
	std::cout << Hash_size << " bytes total space and " << Hash_GB << " GB of space" << std::endl;
	std::cout << "-----------------------------------------------------" << std::endl;
	std::cout << "(if > billions not reliable) Hash space (oracle): " << Hash_oracle_space << " block maxes values" << std::endl;
	std::cout << Hash_size_oracle << " bytes total space and " << Hash_oracle_GB << " GB of space" << std::endl;
	std::cout << "-----------------------------------------------------" << std::endl;

	// variable to maintain the # bytes per bucket
	std::vector<unsigned long long> Bucket_Hash_Bytes(Bucket_size, 0);
	std::vector<double> Bucket_Hash_GBs(Bucket_size, 0);
	double total;

	// print bucketized results
	for (int i=0; i<Bucket_Hash_Values.size(); i++) {
		Bucket_Hash_Bytes.at(i) = Bucket_Hash_Values.at(i)*sizeof(float);
		Bucket_Hash_GBs.at(i) = (double) Bucket_Hash_Bytes.at(i)/GB_convertion;
		total+= Bucket_Hash_GBs.at(i);
		std::cout << "Bucket (" << i << ") \t (2^" << (i+7) << ")\t# block-max values: " << Bucket_Hash_Values.at(i) << "\t # bytes: " << Bucket_Hash_Bytes.at(i) << "\tspace: " << Bucket_Hash_GBs.at(i) << " GBs " << std::endl;
	}
	std::cout << "-----------------------------------------------------" << std::endl;
	std::cout << "Total space in GBs: " << total << std::endl;
}

/* Rule out queries not matching the given policies */
bool QueryLogManager::queriesFileIterator::isSkip() const {
	if(limit && count>=limit)
		return true;
	const size_t sz = queriesP[curr].size();
	if(sz == 0) {
		CERR << "empty query in log: " << curr << Log::endl;
		return true;
	}

	return bucket > 0 ? sz != bucket : false;
}

QueryLogManager::queriesFileIterator& QueryLogManager::queriesFileIterator::operator++() {
	++curr;
	while(curr < queriesP.size() && isSkip()) {
		++curr;
	}
	++count;
	return *this;
}

std::string joinStringsVector(const std::vector<std::string>& terms) {
	std::stringstream buffer;
	std::copy(terms.begin(), terms.end(), std::ostream_iterator<std::string>( buffer ) );
	return buffer.str();
}

float QueryLogManager::score(const std::vector<std::string>& terms) const {
	strToFloatMap::const_iterator it = mapQueryToScore.find(joinStringsVector(terms));
	return ( it != mapQueryToScore.end()) ? ((*it).second) : (-1.0);
}

void QueryLogManager::setScoreForQuery(const std::vector<std::string>& terms, float v) {
	mapQueryToScore[joinStringsVector(terms)] = v;
}

/* Load scores (Ground Truth) from file */
void QueryLogManager::loadScoresFromFile(const char* fname) {
	FILE *fq = fopen(fname,"r");
	std::vector<std::string> terms;
	while(1)  {
		terms.clear();
		char line[5120];

		if( readline(fq, line) == -1 ) //this is the end
			break;

		std::string token, text(line);
		std::istringstream iss(text);
		while ( getline(iss, token, ' ') )
			terms.push_back(token);

		//well, last one is not a term, but a score
		float score = atof(terms[terms.size()-1].c_str());
		terms.pop_back();
		setScoreForQuery(terms,score);
	}
	fclose(fq);
}


/* Usage: Given the list length output the block bits used for DocID-Oriented Block-Max Indexes */
void QueryProcessing::bitOracle(unsigned int& length, int& lenBits) {
	lenBits = intlog2(length);
	std::vector<int> Bucket (18, 0);
	Bucket.at(0) = 10; // 2^7
	Bucket.at(1) = 10; // 2^8
	Bucket.at(2) = 10; // 2^9
	Bucket.at(3) = 10; // 2^10
	Bucket.at(4) = 6; // 2^11 --
	Bucket.at(5) = 6; // 2^12 --
	Bucket.at(6) = 7; // 2^13 -
	Bucket.at(7) = 7; // 2^14 --
	Bucket.at(8) = 7; // 2^15 --
	Bucket.at(9) = 8; // 2^16 -
	Bucket.at(10) = 8; // 2^17 -
	Bucket.at(11) = 7; // 2^18
	Bucket.at(12) = 6; // 2^19 --
	Bucket.at(13) = 6; // 2^20 --
	Bucket.at(14) = 6; // 2^21 --
	Bucket.at(15) = 6; // 2^22 --
	Bucket.at(16) = 6; // 2^23 --
	Bucket.at(17) = 6; // 2^24 --

	// pick the right bucket
	int B_counter = 0;
	for (int i=7; i<25; i++) {
		if (lenBits<i) {
			lenBits = Bucket.at(B_counter);
			break;
		} else
			++B_counter;
	}
}

/* Set the quantiles needed for the Quantized Index
   Input: lptr pointers
   Output: it sets the quantile value per term */
void QueryProcessing::Set_Quantiles(lptrArray& lps) {
	for (int i=0; i<lps.size(); i++)
		lps[i]->quantile = (float) lps.getListMaxScore(i)/CONSTS::Quantization;
	//std::cout << lps.getListMaxScore(i) << " Quantiz: " << CONSTS::Quantization << " computed quantile: " << lps[i]->quantile << std::endl;
}

/* Get the average Document Length by reading the appropriate file */
float QueryProcessing::Get_Average_Document_Length() {
	FILE *handler = fopen(CONSTS::doclenFileName.c_str(), "r");
	if (handler==NULL)
		CERR << CONSTS::doclenFileName.c_str() << " could not be opened" << EFATAL;

	int *doc = new int[CONSTS::MAXD];
	if( fread( doc, sizeof(int), CONSTS::MAXD, handler) != CONSTS::MAXD)
		CERR << "Document File: " << CONSTS::doclenFileName << " does not contain " << CONSTS::MAXD << " values " << EFATAL;

	long long sum;
	for (int i=0; i<CONSTS::MAXD; i++)
		sum += doc[i];

	fclose(handler);
	free(doc);
	float avg = (float) sum/CONSTS::MAXD;
	COUT1 << "Total sum of document lengths: " << sum << " and the average document length: " << avg << Log::endl;
	return avg;
}

/* Print Maxscore List Length correlation stats */
void QueryProcessing::Query_log_Maxscore_Length_Correlation(termsMap& lex) {
	std::vector<int> List_Lengths;
	termsMap queries_loaded;
	int total_query_terms = 0;
	int distinct_query_terms = 0;
	int same_term = 0;
	FILE *query_log_handler = fopen(CONSTS::ikQuery.c_str(),"r");
	while(1)  {
		char line[5120];
		if( readline(query_log_handler, line) == -1 )
			break;

		char *term;
		term = strtok(line," ");

		// First query term in line
		std::unordered_map<std::string, int>::const_iterator term_exists = lex.find(std::string(term));
		if ( term_exists != lex.end() ) { // if term exists
			std::unordered_map<std::string, int>::const_iterator query_exists = queries_loaded.find(std::string(term));
			if ( query_exists == queries_loaded.end() ) { // if term does not exists in our loaded queries
				std::pair<std::string, int> new_entry (term_exists->first, term_exists->second);
				queries_loaded.insert(new_entry);
				//List_Lengths.push_back(term_exists->second);
				++distinct_query_terms;
			} else
				++same_term;
		} else // term does not exist in our lexicon
			CERR << std::string(term) << " was not found in lexicon" << Log::endl;
		++total_query_terms;

		// Remaining query terms in the same line
		while ((term = strtok(NULL," ")) != NULL){
			++total_query_terms;
			std::unordered_map<std::string, int>::const_iterator term1_exists = lex.find(std::string(term));
			if ( term1_exists != lex.end() ) { // if term exists
				std::unordered_map<std::string, int>::const_iterator query_exists = queries_loaded.find(std::string(term));
				if ( query_exists == queries_loaded.end() ) { // if term does not exists in our loaded queries
					std::pair<std::string, int> new_entry (term1_exists->first, term1_exists->second);
					queries_loaded.insert(new_entry);
					//List_Lengths.push_back(term1_exists->second);
					++distinct_query_terms;
				} else
					++same_term;
			} else // term does not exist in our lexicon
				CERR << std::string(term) << " was not found in lexicon" << Log::endl;
		}
	}
	// close handler
	fclose(query_log_handler);

	std::vector<long> Buckets;
	// finer granularity observation of list length distribution
	for (int i=0; i<26; i++)
		Buckets.push_back(pow(2, i));

	std::vector<long long> List_Length_Buckets (Buckets.size(), 0);
	std::vector<float> Maxscore_Sum_Buckets (Buckets.size(), 0.0f);

	int terms = 0;
	long long list_length_sum = 0;
	for(std::unordered_map<std::string, int>::const_iterator it = queries_loaded.begin(); it != queries_loaded.end(); ++it) {
		//std::cout << it->first << " and " << it->second << std::endl;
		std::vector<std::string> term (1, std::string(it->first));
		lptrArray lps = openLists(term, lex);

		++terms;
		list_length_sum += lps[0]->unpadded_list_length;
/*		for (int j=0; j<Buckets.size()-1; j++) {
			if ((it->second >= Buckets.at(j))&&(it->second < Buckets.at(j+1))) {
				List_Length_Buckets.at(j) += 1;
				Maxscore_Sum_Buckets.at(j) += lps[0]->maxScoreOfList;
				break;
			}
		}
*/
	}

/*
	std::cout << "# total terms: " << terms << std::endl;
	std::cout << "# list_length_sum: " << list_length_sum << std::endl;
	std::cout << "avg list_length_sum: " << (float) list_length_sum/terms << std::endl;

	std::cout << "# same terms: " << same_term << std::endl;
	std::cout << "# distinct query terms: " << distinct_query_terms << std::endl;
	std::cout << "# query terms: " << total_query_terms << std::endl;

	for (int j=0; j<Buckets.size(); j++)
		std::cout << "Bucket(" << j << "): #query terms: " << List_Length_Buckets.at(j) << " and sum of maxscore: " << Maxscore_Sum_Buckets.at(j) << std::endl;
*/
}

/* Get distinct terms' list length distribution from the Query log (1000 queries)
   Input: map containing terms and their unpadded list lengths
   Output: vector of int that contains all distinct's terms list length from the query log */
std::vector<int> QueryProcessing::Load_Query_log_List_lengths(termsMap& lex) {
	// load query file
	std::vector<int> List_Lengths;
	termsMap queries_loaded;
	int total_query_terms = 0;
	int distinct_query_terms = 0;
	FILE *query_log_handler = fopen(CONSTS::ikQuery.c_str(),"r");
	while(1)  {
		char line[5120];
		if( readline(query_log_handler, line) == -1 )
			break;

		char *term;
		term = strtok(line," ");

		// First query term in line
		std::unordered_map<std::string, int>::const_iterator term_exists = lex.find(std::string(term));
		if ( term_exists != lex.end() ) { // if term exists
			//std::unordered_map<std::string, int>::const_iterator query_exists = queries_loaded.find(std::string(term));
			//if ( query_exists == queries_loaded.end() ) { // if term does not exists in our loaded queries
			//	std::pair<std::string, int> new_entry (term_exists->first, term_exists->second);
			//	queries_loaded.insert(new_entry);
			List_Lengths.push_back(term_exists->second);
			//	++distinct_query_terms;
			//}
		} else // term does not exist in our lexicon
			CERR << std::string(term) << " was not found in lexicon" << Log::endl;

		// Remaining query terms in the same line
		while ((term = strtok(NULL," ")) != NULL){
			std::unordered_map<std::string, int>::const_iterator term1_exists = lex.find(std::string(term));
			if ( term1_exists != lex.end() ) { // if term exists
				//std::unordered_map<std::string, int>::const_iterator query_exists = queries_loaded.find(std::string(term));
				//if ( query_exists == queries_loaded.end() ) { // if term does not exists in our loaded queries
				//	std::pair<std::string, int> new_entry (term1_exists->first, term1_exists->second);
				//	queries_loaded.insert(new_entry);
				List_Lengths.push_back(term1_exists->second);
				//	++distinct_query_terms;
				//}
			} else // term does not exist in our lexicon
				CERR << std::string(term) << " was not found in lexicon" << Log::endl;
		}
	}
	// close handler
	fclose(query_log_handler);
	COUT2 << "Total distinct terms in query trace that exist in our lexicon: " << distinct_query_terms << Log::endl;
	return List_Lengths;
}

/* Generate a block max array given a lptr
   Input: lptr of a specific term and max array to store the block maxes
   Output: vector of block maxes
   Note: Assuming we initialized pages[] that that we have the document length of all documents */
void QueryProcessing::on_the_fly_max_array_generation(lptr*& lps, std::vector<float>& max_array, int& bits) {
	int total_blocks = (lps->lengthOfList/CONSTS::BS); // max blocks based on pfd compression
	int cur_block_number = 0; // hash block number
	float score = 0.0f;

	// for all docids of lps less than CONSTS::MAXD, so as not to count dids out of range and junk scores
	// increase by one so that in the next round, the method skipToDidBlockAndDecode to decode the next block
	for (int blocks=0; blocks<total_blocks; blocks++) {
		// for all blocks except the first one, decode block
		if ( blocks != 0) {
			// decode block
			lps->skipToDidBlockAndDecode(lps->currMax+1);
		}

		// if docid in docid range
		if (lps->did < CONSTS::MAXD) {
			// get block number
			cur_block_number = lps->did >> bits;
			// compute score
			score = lps->calcScore(lps->getFreq(), pages[lps->did]);
			// check if we need to update current value of block
			if ( Fcompare(max_array[cur_block_number], score) == -1 )
				max_array[cur_block_number] = score;
		} else
			break;

		// for all dids of block, compute did (prefix sum), and push to vector
		for (lps->elem+=1; lps->did<lps->currMax; lps->elem++) {
			lps->did +=lps->dbuf[lps->elem];

			// if docid in docid range
			if (lps->did < CONSTS::MAXD) {
				// get block number
				cur_block_number = lps->did >> bits;
				// compute score
				score = lps->calcScore(lps->getFreq(), pages[lps->did]);
				// check if we need to update current value of block
				if ( Fcompare(max_array[cur_block_number], score) == -1 )
					max_array[cur_block_number] = score;
			} else
				break;
		}
	}

	// reset pointers, so that we can use the lptr structure again in query processing (assuming that was a preprocessing step)
	lps->reset_list();
}

/* Usage: it generates a RawIndexList from lps structure (faster)
   Input: lps structure for the specific term
   Output: RawIndexList structure */
RawIndexList QueryProcessing::lps_to_RawIndexList(lptr*& lps) {
	BasicList B_Term(lps->term, lps->termId);
	RawIndexList Raw_list(B_Term);
	int total_blocks = (lps->lengthOfList/CONSTS::BS); // lengtoflist is the padded one, so if we divide by bs the result is always an integer value

	//debug
	//std::cout << lps->maxScoreOfList <<" " << lps->did << " \t freq: " << lps->getFreq() << "\tscore: " << lps->calcScore(lps->getFreq(), pages[lps->did]) <<  " and element: " << lps->elem << " and curMax " << lps->currMax << std::endl;

	// for all docids of lps (lengthofList has plus the padded ones, but it's ok)
	// increase by one so that in the next round, the method skipToDidBlockAndDecode to decode the next block
	for (int blocks=0; blocks<total_blocks; blocks++) {
		// for all blocks except the first one, decode block
		if ( blocks != 0) {
			// decode block
			lps->skipToDidBlockAndDecode(lps->currMax+1);
		}

		Raw_list.doc_ids.push_back(lps->did);
		Raw_list.freq_s.push_back(lps->getFreq());
		Raw_list.scores.push_back(lps->calcScore(lps->getFreq(), pages[lps->did]));

		// for all dids of block, compute did (prefix sum), and push to vector
		for (lps->elem+=1; lps->did<lps->currMax; lps->elem++) {
			lps->did +=lps->dbuf[lps->elem];

			// fill vectors
			Raw_list.doc_ids.push_back(lps->did);
			Raw_list.freq_s.push_back(lps->getFreq());
			Raw_list.scores.push_back(lps->calcScore(lps->getFreq(), pages[lps->did]));
		}
	}

	// optional: fix padded scores to 0.0f
	for (int i=lps->unpadded_list_length; i<lps->lengthOfList; i++)
		Raw_list.scores.at(i) = (0.0f);

	// reset pointers, so that we can use the lptr structure
	lps->reset_list();

	return Raw_list;
}


// Usage: it generates a RawIndexList from lps structure (plus the padded info) (if you want unpadded, change in for loop the length_list with CONSTS::MAXD
// Input: lps structure for the specific term
// Output: RawIndexList structure
RawIndexList QueryProcessing::naive_lps_to_RawIndexList(lptr*& lps) {
	BasicList B_Term(lps->term, lps->termId);
	RawIndexList Raw_list(B_Term);
	int dids_retrieved = 0;
	int max_did = (CONSTS::MAXD+CONSTS::BS);
	for (int did=0; (dids_retrieved<lps->lengthOfList)&&(did<max_did); ++did) {
		lps->did = lps->nextGEQ(did);
		// if did exists push
		if (lps->did == did) {
			Raw_list.doc_ids.push_back(did);
			Raw_list.freq_s.push_back(lps->getFreq());
			Raw_list.scores.push_back(lps->calcScore(lps->getFreq(), did));
			++dids_retrieved;
		}
	}
	return Raw_list;
}

// Usage: it sets metainfo about layering to lps structure, such as has_layers, is_essential boolean, list_id and computes paddded list length
// Input: lists, list_ids vector, mapping of terms to unpadded_length
void QueryProcessing::set_layered_metainfo_to_lps(lptrArray& lps, std::vector<int>& list_ids, termsMap& Unpadded_Length_Map, std::vector<std::string>& layered_terms) {
	int list_counter = 0;
	for (int i=0; i<list_ids.size(); i+=2) {
		// check if bad term exists and if so it means that we have layer
		if (list_ids.at(i+1) != -1) {
			// good term
			lps[list_counter]->has_layers = true;
			lps[list_counter]->list_id = list_ids.at(i);
			lps[list_counter]->is_essential = true;
			lps[list_counter]->layer_status = 0;

			// find good term in our terms to unpadded length map and set correct unpadded list length and padded list length for the good term
			std::unordered_map<std::string, int>::const_iterator good_term_exists = Unpadded_Length_Map.find(layered_terms.at(list_counter));
			if ( good_term_exists != Unpadded_Length_Map.end() ) {
				lps[list_counter]->unpadded_list_length = good_term_exists->second;
				lps[list_counter]->lengthOfList = lps[list_counter]->unpadded_list_length + (CONSTS::BS - (lps[list_counter]->unpadded_list_length%CONSTS::BS));
			}

			// increase counter
			++list_counter;

			// bad term
			lps[list_counter]->has_layers = true;
			lps[list_counter]->list_id = list_ids.at(i+1);
			lps[list_counter]->is_essential = true;
			lps[list_counter]->layer_status = 0;

			// find bad term in our terms to unpadded length map and set correct unpadded list length and padded list length for the bad term
			std::unordered_map<std::string, int>::const_iterator bad_term_exists = Unpadded_Length_Map.find(layered_terms.at(list_counter));
			if ( bad_term_exists != Unpadded_Length_Map.end() ) {
				lps[list_counter]->unpadded_list_length = bad_term_exists->second;
				lps[list_counter]->lengthOfList = lps[list_counter]->unpadded_list_length + (CONSTS::BS - (lps[list_counter]->unpadded_list_length%CONSTS::BS));
			}
		} else { // only original term exists
			lps[list_counter]->has_layers = false;
			lps[list_counter]->list_id = list_ids.at(i);
			lps[list_counter]->is_essential = true;
			lps[list_counter]->layer_status = 5;
		}

		// increase counter
		++list_counter;
	}
}

/* Translation of the query terms to layered terms and assigns ids
   Input: vector of string containing query terms, lexicon, vector to store translated terms, vector to store ids */
void QueryProcessing::translate_and_assign_ids_to_layered_terms(std::vector<std::string>& word_l, termsMap& lex, std::vector<string>& translated_terms, std::vector<int>& list_ids) {
	int list_id_counter = 0;
	const int no_bad_layer = -1;

	for (int i=0; i<word_l.size(); i++) {
		std::unordered_map<std::string, int>::const_iterator term_exists = lex.find(word_l.at(i));

		// if term does not exist, it means we have layers so push both good and bad term to our query vector
		if ( term_exists == lex.end() ) {
			// add translated terms into vector
			translated_terms.push_back(word_l.at(i)+CONSTS::GOOD_TERM_SUFFIX);
			translated_terms.push_back(word_l.at(i)+CONSTS::BAD_TERM_SUFFIX);

			// assign ids to list (even = good lists and odd bad lists)
			list_ids.push_back(list_id_counter);
			list_ids.push_back(++list_id_counter);
		} else { // else push the original query term, since no layers were created for this term
			// add terms as is, since term has no layers
			translated_terms.push_back(word_l.at(i));

			// assign accordingly the ids
			list_ids.push_back(list_id_counter);
			list_ids.push_back(no_bad_layer);

			// increase counter
			++list_id_counter;
		}
		++list_id_counter;
	}
}

/* Generates a block max array
   Input: vector of dids, scores and vector for output
   Output: vector of block maxes */
void QueryProcessing::on_the_fly_max_array_generation(std::vector<float>& max_array, RawIndexList& rilist, int& bits) {
	int block_number = 0;
	// for all docids
	for (int i=0; i<rilist.doc_ids.size(); i++) {
		if (rilist.doc_ids.at(i)<CONSTS::MAXD) {
			// get block number
			block_number = rilist.doc_ids.at(i) >> bits;
			// check if we need to update current value of block
			if ( Fcompare(max_array.at(block_number), rilist.scores.at(i)) == -1 )
				max_array.at(block_number) = rilist.scores.at(i);
		}
	}
}

/* Load the unpadded postings of the loaded index
   Input: path of file to load the real postings
   Output: map of <term, # real postings> */
std::unordered_map<std::string, int> QueryProcessing::Load_unpadded_postings(const std::string& path) {
	std::unordered_map<std::string, int> real_postings;
	FILE* real_postings_handler = fopen(path.c_str(),"r");
	if(real_postings_handler==NULL)
		CERR << path << " could not be opened" << EFATAL;

	char term[1024];
	int length;
	while( fscanf(real_postings_handler,"%s\t%d\n",term, &length)!= EOF )
		real_postings[string(term)] = length;

	fclose(real_postings_handler);
	assert(real_postings.size());
	return real_postings;
}

/* Load from the raw index the list lengths into a vector
   Input: nothing
   Output: vector of list lengths of all terms in our index */
// std::vector<int> QueryProcessing::Load_Entire_Index_List_Lengths() {
// 	std::vector<int> List_Lengths;
// 	std::string path = CONSTS::trecRawRoot + CONSTS::MERGED_BOOL_PATH + CONSTS::INFO_INDEX;
// 	FILE* List_Length_Handler = fopen(path.c_str(), "r");
// 	if (List_Length_Handler==NULL)
// 		CERR << path << " could not be opened" << EFATAL;

// 	int number_of_terms_in_collection;
// 	fread( &number_of_terms_in_collection, sizeof(int), 1, List_Length_Handler); // first int contains the # of terms in the entire collection

// 	unsigned int *inf_buffer = new unsigned int[ 4 * number_of_terms_in_collection];
// 	fread( inf_buffer, sizeof(int), 4 * number_of_terms_in_collection, List_Length_Handler);  // read all inf file

// 	for (int i=0; i<number_of_terms_in_collection; i++)
// 		List_Lengths.push_back(inf_buffer[i*4+1]); // the unpadded size is the second integer from the 4-tuple we load		fclose(List_Length_Handler);

// 	// print report
// 	COUT2<<"There are "<<number_of_terms_in_collection<<" terms in our index"<<Log::endl;
// 	COUT2 << List_Lengths.size() << " distinct terms loaded in vector";
// 	return List_Lengths;
// }

/* List length Distribution for the entire index */
void QueryProcessing::List_Length_Distribution(std::vector<int>& List_Lengths) {
	std::vector<long> Buckets;
	// finer granularity observation of list length distribution
	for (int i=0; i<26; i++)
		Buckets.push_back(pow(2, i));

	std::vector<long long> List_Length_Buckets (Buckets.size(), 0);

	long total_distinct_terms = 0;
	for (long i=0; i<List_Lengths.size(); i++) {
		for (int j=0; j<Buckets.size()-1; j++) {
			if ((List_Lengths.at(i) >= Buckets.at(j))&&(List_Lengths.at(i) < Buckets.at(j+1))) {
				List_Length_Buckets.at(j) += 1;
				total_distinct_terms++;
				break;
			}
		}
	}

	// Print Results
	std::cout << "-----------------------------------------------------" << std::endl;
	std::cout << "Total number of distinct terms in the index: " << total_distinct_terms << std::endl;
	for (int i=0; i<Buckets.size()-1; i++)
		std::cout << "# Lists with size [" << Buckets.at(i) << ", " << Buckets.at(i+1) << ") : " << List_Length_Buckets.at(i) << std::endl;
	std::cout << "-----------------------------------------------------" << std::endl;
}

/* Get List Lengths */
// Input: unordered map <string, int> lex
// Output: vector of list lengths of all terms in our loaded index
std::vector<int> QueryProcessing::Get_List_Lengths(termsMap& lex) {
	std::vector<int> List_Lengths;
	for (termsMap::iterator it = lex.begin(); it != lex.end(); ++it)
		List_Lengths.push_back(it->second);

	return List_Lengths;
}
