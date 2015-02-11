/*
 * globals.h
 *
 *  Created on: Nov 19th, 2014
 *      Author: qi
 */

/*
Whatever you do, create baby lex first, use ./qp -b; Modify two places, one is the lex file dir, the other is the number of terms in the lex
*/

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <ostream>
#include <vector>
#include <assert.h>
#include <stdint.h>

#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

// typedefs
typedef unsigned int uint;
typedef unsigned long ulong;

//#define PROFILING
#ifdef GPROFILING
#define MYNOINLINE __attribute__((noinline))
#define PROFILER(a) //
#else
#define MYNOINLINE
#define PROFILER(a) profilerC::getInstance().stepCounter(a);
#endif

#ifdef CPP0X
	#include <unordered_map>
	#define hashMapStd std::unordered_map
	typedef std::unordered_map<std::string, int> termsMap;
#else
	#include <map>
	typedef std::map<std::string, int> termsMap;
	#define hashMapStd std::map
#endif




namespace Log {
	enum verbLevels { VALL, VDEBUG, VPROFILING, VOUTPUT };
#ifdef GPROFILING
#define COUT1 //
#define COUT2 //
#else
	#define COUT Log::logger()
	#define COUT1 Log::logger() << Log::verb<Log::VDEBUG>
	#define COUT2 Log::logger() << Log::verb<Log::VPROFILING>
	#define COUT3 Log::logger() << Log::verb<Log::VOUTPUT>
	#define COUT4 Log::logger(std::cout,4)
	#define COUT5 Log::logger(std::cout,5)
	#define CERR Log::logger(std::cerr,100)
	#define EFATAL Log::endl << " " << __FILE__ << " " << __LINE__ << Log::fatal
#endif
	void setGlobalVerbosityForAllLoggers(int v); //{ logger::GLOBAL_VERBOSITY = v; }

	typedef std::ostream outstrType; //just in case

	class logger {
		outstrType& out; //the stream of this logger
		int verbosity; //current verbosity
		static int GLOBAL_VERBOSITY;
	public:

		logger(outstrType& s=std::cout, int v=0) : out(s), verbosity(v) {}
		static void setGlobalV(int v) { GLOBAL_VERBOSITY = v; }
		template <typename T>
		logger& operator<<(const T& rhs) {
			if(verbosity >= GLOBAL_VERBOSITY)
			  out <<  rhs;
			return *this;
		 }

		//manipulators of the stream types
		typedef logger& (*endlType)(logger&);
		typedef logger& (*setVerbType)(logger&, int& v);
		logger& operator<<(endlType manip){ return manip(*this); }
		logger& operator<<(setVerbType manip){ return manip(*this,verbosity); }

		void setV(int v) { verbosity = v;}
		void flush() { if(verbosity >= GLOBAL_VERBOSITY) out << std::endl; }
	};

	logger& endl(logger& stream);  //stream flush
	logger& fatal(logger& stream);  //print and throw

	template<int val> //verbosity setter
	logger& verb(logger& stream, int& v) {
		v=val;
		return stream;
	}


/*

  usage
	COUT << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
	Log::setGlobalVerbosityForAllLoggers(1);
	CERR << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
	log::logger() << log::verb<5> << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
	COUT5 << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
*/

};

namespace CONSTS {

	enum counters { CNT_EVICTIONS,
		WAND,GETBLOCKSCORE,GETBLOCKSCORE1, BLOCKFAIL,ALIGN,MISALIGNED, GETFREQ, EVAL,HEAPIFY,HEAPIFYOK,NEXTGEQ,NEXTGEQ1,ALLQS,OTFBMG,
		SKIPS, SKIPSB, DECODE, CAND, MISTRUST, ESSENTIAL, SHALLOWPOINTER, SHALLOWPOINTER1, SORT, DOCIDFILTER1, DOCIDFILTER2,
		EARLYTERMINATION0, EARLYTERMINATION01, EARLYTERMINATION1, EARLYTERMINATION2, EARLYTERMINATION3, EARLYTERMINATION4, EARLYTERMINATION5,
		NONESSENTIAL,STEP1, STEP2, STEP3, NOP };

	const std::string QLogPath("../QueryLog/");
	const std::string fResultLog("../500_or_5");
	const std::string fQueryScore("../query_score");
	const std::string Golden_Threshold("../QueryLog/Golden_Threshold_trec06"); // Golden_Threshold_trec06 // Golden_Threshold_trec05
	const std::string zeroMapping("../QueryLog/10000q_terms.mapping"); //10000q_terms.mapping CHANGED !!
	const std::string termsMapping("../QueryLog/10000q_terms.mapping");
	// 10000q_terms.mapping // normal and reordered
	// 10000q_terms.mapping_layered // layered
	// 10000q_terms.mapping_layered_sorted // reordered + layered
	// version2: 10000q_terms_mapping_layered_version2 // layered (version2)
	// TREC05
	// trec05_mapping
	const std::string termsMapping_Layer("../QueryLog/10000q_terms.mapping_layered");

	// not used
	const std::string Index_term_mapping("word_file");
	const std::string entireIndexListLengths("../QueryLog/Entire_index_list_lengths");
	const std::string Unpadded_10k_postings("../QueryLog/10k_pair_term_unpadded_postings");
	const std::string Entire_Index_Maxscore_Unpadded_List_length("../QueryLog/Entire_Index_Maxscore_Unpadded_List_length"); //../QueryLog/
	const std::string Starting_Threshold("../QueryLog/10th_Starting_Threshold");
	const std::string Top1_Threshold("../QueryLog/Top-1_Threshold");
	const std::string Top10_Threshold("../QueryLog/Top-10_Threshold");

	const std::string ikQueryScores("../QueryLog/1k_scores_float");  // 1k_scores
	const std::string ikQuery("../QueryLog/1000query"); //1k_trec05; //1000query
//	const std::string trecRoot("/research/index/trec06/"); //for trec on dodo


	/*dodo*/
	const std::string trec_output("/home/qw376/index_trec_06/"); //to generate babyindex
	const std::string trecRoot("/home/qw376/index_trec_06/"); //main index 
	// const std::string trec_output("/home/qw376/index_for_LM_11.29/"); //to generate babyindex
	// const std::string trecRoot("/home/qw376/index_for_LM_11.29/"); //main index 
	// const std::string testingQuery("/home/vgc/qi/11.10/PairQueriesTrainingSet_10.11_Sorted_Uniq_Nodup_Nosym");
	// const std::string testingQuery("/home/vgc/qi/11.10/bigram_FinalPair_withoutdups_10.31");
	// const std::string testingQuery("/home/vgc/qi/11.10/unique_testing_pairs");
	// const std::string testingQuery("/home/vgc/qi/11.10/query_id_mapping_10.22_10k_testing");
	// const std::string testingQuery("/home/qw376/1.28/trec06_training_lm_unique");
	// const std::string testingQuery("/home/qw376/pair_lex/bigram_drop_prob_12.1_nodup");
	// const std::string testingQuery("/home/qw376/1.28/trec06_testing_lm_unique");
	// const std::string testingQuery("/home/qw376/special_query_trace/09query_qid"); 
	const std::string testingQuery("/home/qw376/1.28/trec_06_testing_12.12");
	// const std::string testingQuery("/home/qw376/1.28/check_query");
	// const std::string testingQuery("/home/qw376/1.28/trec_06_testing_actualqueries_11.22_no_colon"); //main query
	// const std::string testingQuery("/home/qw376/1.28/trec_06_training_actualqueries_11.22_no_colon"); //main query
	const std::string doclenFileName("/home/qw376/Info_Clueweb/documentLengths"); //doc lengths
	// const std::string basic_table("/home/vgc/qi/experiments/lexicon_baby_712"); //lexicon
	// const std::string basic_table("/home/vgc/qi/11.10/lexicon_baby_testing_11.3"); //lexicon
	const std::string basic_table("/home/qw376/1.28/baby_lexicon_all_trec06_11.19"); //lexicon 
	// const std::string basic_table("/home/qw376/1.28/baby_lexicon_LM_11.29"); //lexicon 
	const uint number_of_terms = 37200; //at CluewebReader.cpp line 36
	// const uint number_of_terms = 1780695; //at CluewebReader.cpp line 36
	const std::string top_layer_index_did("/home/qw376/tl_trec06_all_did/");//impact sorted binary dir
	const std::string top_layer_index_score("/home/qw376/tl_trec06_all_score/");//impact sorted binary dir
	// const std::string intersection_index_did("/home/qw376/pair_and_trec06_training_did/");//impact sorted binary dir
	const std::string intersection_index_did("/home/qw376/pair_and_trec06_testing_did/");//impact sorted binary dir and
	// const std::string union_index_did("/home/qw376/or_pair_and_trec06_training_did/");//impact sorted binary dir or
	// const std::string intersection_index_score_1("/home/qw376/pair_and_trec06_training_score_first/");//impact sorted binary dir term1
	const std::string intersection_index_score_first("/home/qw376/pair_and_trec06_testing_score_first/");//impact sorted binary dir term1
	// const std::string intersection_index_score_2("/home/qw376/pair_and_trec06_training_score_second/");//impact sorted binary dir term2
	const std::string intersection_index_score_second("/home/qw376/pair_and_trec06_testing_score_second/");//impact sorted binary dir term2
	// const std::string pair_lex("/home/qw376/pair_lex/trec_06_testing_lex");//pair lex
	// const std::string pair_lex("/home/qw376/pair_lex/trec_06_training_lex");//pair lex
	const std::string pair_lex("/home/qw376/pair_lex/LM_intersection_size_12.1_4");//pair lex
	const std::string cr_results("/home/qw376/complex_ranker_scores/cr_trec06_testing_actual_queries_11.22");//cr results for testing
	// const std::string cr_results("/home/qw376/complex_ranker_scores/cr_trec06_training_actual_queries_11.22");//cr results for training
	// const std::string index_depth("/home/qw376/depth_index/layer_trec_06_testing_11.24");//depth index
	const std::string pair_depth("/home/qw376/depth_index/pair_trec_06_testing_11.24");//depth for testing pairs
	const std::string num_of_pairs("/home/qw376/depth_index/testing_num_of_pairs_11.24");//depth for testing pairs
	// const std::string pair_depth("/home/qw376/depth_index/pair_trec_06_training_11.24_10");//depth for training pairs
	// const std::string num_of_pairs("/home/qw376/depth_index/training_num_of_pairs_11.24_10");//depth for training pairs
	// const std::string index_kl("/home/qw376/experiments/index_kl_0_1000");//known index
	const std::string clueweb_index("/home/qw376/Info_Clueweb/InvertedIndex");//clueweb whole index
	//const std::string Candidates_Pool("/home/qw376/1.28/Candidates_R_2");//Candidates index
	const std::string top_layer_cutoffs = "/home/qw376/trec06_tl_cutoffs/";
	const std::string term_pair_cutoffs = "/home/qw376/trec06_pairs_cutoffs/";
	//const std::string cutoffs_path = "/home/qw376/cutoffs/";
	const std::string percent_of_index_str = "50"; 


	/*moa*/
	// // const std::string trec_output("/home/qw376/index_trec_06/"); //to generate babyindex
	// // const std::string trecRoot("/home/qw376/index_trec_06/"); //main index 
	// const std::string trec_output("/home/qi/index_for_LM_11.29/"); //to generate babyindex
	// const std::string trecRoot("/home/qi/index_for_LM_11.29/"); //main index 
	// // const std::string testingQuery("/home/vgc/qi/11.10/PairQueriesTrainingSet_10.11_Sorted_Uniq_Nodup_Nosym");
	// // const std::string testingQuery("/home/vgc/qi/11.10/bigram_FinalPair_withoutdups_10.31");
	// // const std::string testingQuery("/home/vgc/qi/11.10/unique_testing_pairs");
	// // const std::string testingQuery("/home/vgc/qi/11.10/query_id_mapping_10.22_10k_testing");
	// // const std::string testingQuery("/home/qw376/1.28/trec06_training_lm_unique");
	// const std::string testingQuery("/home/qi/pair_lex/bigram_drop_prob_12.1_nodup");
	// // const std::string testingQuery("/home/qw376/1.28/trec06_testing_lm_unique");
	// // const std::string testingQuery("/home/qw376/1.28/trec_06_testing_actualqueries_11.22_no_colon"); //main query
	// // const std::string testingQuery("/home/qw376/1.28/trec_06_training_actualqueries_11.22_no_colon"); //main query
	// const std::string doclenFileName("/home/qi/Info_Clueweb/documentLengths"); //doc lengths
	// // const std::string basic_table("/home/vgc/qi/experiments/lexicon_baby_712"); //lexicon
	// // const std::string basic_table("/home/vgc/qi/11.10/lexicon_baby_testing_11.3"); //lexicon
	// // const std::string basic_table("/home/qw376/1.28/baby_lexicon_all_trec06_11.19"); //lexicon 
	// const std::string basic_table("/home/qi/1.28/baby_lexicon_LM_11.29"); //lexicon 
	// // const uint number_of_terms = 37200; //at CluewebReader.cpp line 36
	// const uint number_of_terms = 1780695; //at CluewebReader.cpp line 36
	// const std::string top_layer_index_did("/home/qi/tl_trec06_all_did/");//impact sorted binary dir
	// // const std::string top_layer_index_score("/home/qw376/top_layer_index_score/");//impact sorted binary dir
	// // const std::string intersection_index_did("/home/qw376/pair_and_trec06_training_did/");//impact sorted binary dir
	// // const std::string intersection_index_did("/home/qw376/pair_and_trec06_testing_did/");//impact sorted binary dir and
	// const std::string union_index_did("/home/qi/or_pair_and_trec06_training_did/");//impact sorted binary dir or
	// // const std::string intersection_index_score_1("/home/qw376/pair_and_trec06_training_score_first/");//impact sorted binary dir term1
	// // const std::string intersection_index_score_1("/home/qw376/pair_and_trec06_testing_score_first/");//impact sorted binary dir term1
	// // const std::string intersection_index_score_2("/home/qw376/pair_and_trec06_training_score_second/");//impact sorted binary dir term2
	// // const std::string intersection_index_score_2("/home/qw376/pair_and_trec06_testing_score_second/");//impact sorted binary dir term2
	// // const std::string pair_lex("/home/qw376/pair_lex/trec_06_testing_lex");//pair lex
	// // const std::string pair_lex("/home/qw376/pair_lex/trec_06_training_lex");//pair lex
	// const std::string pair_lex("/home/qi/pair_lex/LM_intersection_size_12.1_6");//pair lex
	// const std::string cr_results("/home/qi/complex_ranker_scores/cr_trec06_testing_actual_queries_11.22");//cr results for testing
	// // const std::string cr_results("/home/qw376/complex_ranker_scores/cr_trec06_training_actual_queries_11.22");//cr results for training
	// // const std::string index_depth("/home/qw376/depth_index/layer_trec_06_testing_11.24");//depth index
	// const std::string pair_depth("/home/qi/depth_index/pair_trec_06_testing_11.24");//depth for testing pairs
	// const std::string num_of_pairs("/home/qi/depth_index/testing_num_of_pairs_11.24");//depth for testing pairs
	// // const std::string pair_depth("/home/qw376/depth_index/pair_trec_06_training_11.24_10");//depth for training pairs
	// // const std::string num_of_pairs("/home/qw376/depth_index/training_num_of_pairs_11.24_10");//depth for training pairs
	// // const std::string index_kl("/home/qw376/experiments/index_kl_0_1000");//known index
	// const std::string clueweb_index("/home/qi/Info_Clueweb/InvertedIndex");//clueweb whole index
	// // const std::string Candidates_Pool("/home/qw376/11.10/Candidates_test_g");//Candidates index

	/*pangolin*/
	// // const std::string trec_output("/home/qw376/index_trec_06/"); //to generate babyindex
	// // const std::string trecRoot("/home/qw376/index_trec_06/"); //main index 
	// const std::string trec_output("/data/qw376/index_for_LM_11.29/"); //to generate babyindex
	// const std::string trecRoot("/data/qw376/index_for_LM_11.29/"); //main index 
	// // const std::string testingQuery("/home/vgc/qi/11.10/PairQueriesTrainingSet_10.11_Sorted_Uniq_Nodup_Nosym");
	// // const std::string testingQuery("/home/vgc/qi/11.10/bigram_FinalPair_withoutdups_10.31");
	// // const std::string testingQuery("/home/vgc/qi/11.10/unique_testing_pairs");
	// // const std::string testingQuery("/home/vgc/qi/11.10/query_id_mapping_10.22_10k_testing");
	// // const std::string testingQuery("/home/qw376/1.28/trec06_training_lm_unique");
	// const std::string testingQuery("/data/qw376/pair_lex/bigram_drop_prob_12.1_nodup");
	// // const std::string testingQuery("/home/qw376/1.28/trec06_testing_lm_unique");
	// // const std::string testingQuery("/home/qw376/1.28/trec_06_testing_actualqueries_11.22_no_colon"); //main query
	// // const std::string testingQuery("/home/qw376/1.28/trec_06_training_actualqueries_11.22_no_colon"); //main query
	// const std::string doclenFileName("/data/constantinos/Index/documentLengths"); //doc lengths
	// // const std::string basic_table("/home/vgc/qi/experiments/lexicon_baby_712"); //lexicon
	// // const std::string basic_table("/home/vgc/qi/11.10/lexicon_baby_testing_11.3"); //lexicon
	// // const std::string basic_table("/home/qw376/1.28/baby_lexicon_all_trec06_11.19"); //lexicon 
	// const std::string basic_table("/data/qw376/1.28/baby_lexicon_LM_11.29"); //lexicon 
	// // const uint number_of_terms = 37200; //at CluewebReader.cpp line 36
	// const uint number_of_terms = 1780695; //at CluewebReader.cpp line 36
	// const std::string top_layer_index_did("/data/qw376/tl_trec06_all_did/");//impact sorted binary dir
	// // const std::string top_layer_index_score("/home/qw376/top_layer_index_score/");//impact sorted binary dir
	// // const std::string intersection_index_did("/home/qw376/pair_and_trec06_training_did/");//impact sorted binary dir
	// // const std::string intersection_index_did("/home/qw376/pair_and_trec06_testing_did/");//impact sorted binary dir and
	// const std::string union_index_did("/data/qw376/or_pair_and_trec06_training_did/");//impact sorted binary dir or
	// // const std::string intersection_index_score_1("/home/qw376/pair_and_trec06_training_score_first/");//impact sorted binary dir term1
	// // const std::string intersection_index_score_1("/home/qw376/pair_and_trec06_testing_score_first/");//impact sorted binary dir term1
	// // const std::string intersection_index_score_2("/home/qw376/pair_and_trec06_training_score_second/");//impact sorted binary dir term2
	// // const std::string intersection_index_score_2("/home/qw376/pair_and_trec06_testing_score_second/");//impact sorted binary dir term2
	// // const std::string pair_lex("/home/qw376/pair_lex/trec_06_testing_lex");//pair lex
	// // const std::string pair_lex("/home/qw376/pair_lex/trec_06_training_lex");//pair lex
	// const std::string pair_lex("/data/qw376/pair_lex/LM_intersection_size_12.1_10");//pair lex
	// const std::string cr_results("/data/qw376/complex_ranker_scores/cr_trec06_testing_actual_queries_11.22");//cr results for testing
	// // const std::string cr_results("/home/qw376/complex_ranker_scores/cr_trec06_training_actual_queries_11.22");//cr results for training
	// // const std::string index_depth("/home/qw376/depth_index/layer_trec_06_testing_11.24");//depth index
	// const std::string pair_depth("/data/qw376/depth_index/pair_trec_06_testing_11.24");//depth for testing pairs
	// const std::string num_of_pairs("/data/qw376/depth_index/testing_num_of_pairs_11.24");//depth for testing pairs
	// // const std::string pair_depth("/home/qw376/depth_index/pair_trec_06_training_11.24_10");//depth for training pairs
	// // const std::string num_of_pairs("/home/qw376/depth_index/training_num_of_pairs_11.24_10");//depth for training pairs
	// // const std::string index_kl("/home/qw376/experiments/index_kl_0_1000");//known index
	// const std::string clueweb_index("/data/constantinos/InvertedIndex");//clueweb whole index
	// // const std::string Candidates_Pool("/home/qw376/11.10/Candidates_test_g");//Candidates index


	/*vida*/
	// // const std::string trec_output("/home/vgc/qi/trec_output/"); //to generate babyindex
	// // const std::string trecRoot("/home/vgc/qi/trec_output/"); //main index
	// const std::string trec_output("/home/vgc/qi/index_trec_06/"); //to generate babyindex
	// const std::string trecRoot("/home/vgc/qi/index_trec_06/"); //main index
	// // const std::string testingQuery("/home/vgc/qi/11.10/PairQueriesTrainingSet_10.11_Sorted_Uniq_Nodup_Nosym");
	// // const std::string testingQuery("/home/vgc/qi/11.10/bigram_FinalPair_withoutdups_10.31");
	// // const std::string testingQuery("/home/vgc/qi/11.10/unique_testing_pairs");
	// // const std::string testingQuery("/home/vgc/qi/11.10/query_id_mapping_10.22_10k_testing");
	// // const std::string testingQuery("/home/vgc/qi/1.28/trec06_training_lm_unique");
	// // const std::string testingQuery("/home/vgc/qi/1.28/trec06_testing_lm_unique");
	// const std::string testingQuery("/home/vgc/qi/pair_lex/bigram_drop_prob_12.1_nodup");
	// // const std::string testingQuery("/home/vgc/qi/11.10/query_id_mapping_10.9_50k"); //main query
	// // const std::string testingQuery("/data/qw376/experiments/test_query"); //main query
	// const std::string doclenFileName("/home/vgc/qi/Info_Clueweb/documentLengths"); //doc lengths
	// // const std::string basic_table("/home/vgc/qi/experiments/lexicon_baby_712"); //lexicon
	// // const std::string basic_table("/home/vgc/qi/11.10/lexicon_baby_testing_11.3"); //lexicon
	// // const std::string basic_table("/home/vgc/qi/1.28/baby_lexicon_all_trec06_11.19"); //lexicon
	// const std::string basic_table("/home/vgc/qi/1.28/baby_lexicon_LM_11.29"); //lexicon 
	// // const uint number_of_terms = 37200; //at CluewebReader.cpp line 36
	// const uint number_of_terms = 1780695; //at CluewebReader.cpp line 36
	// const std::string top_layer_index_did("/home/vgc/qi/top_layer_index_did/");//impact sorted binary dir
	// // const std::string top_layer_index_score("/home/vgc/qi/top_layer_index_score/");//impact sorted binary dir
	// // const std::string intersection_index_did("/home/vgc/qi/pair_and_trec06_training_did/");//impact sorted binary dir
	// const std::string intersection_index_did("/home/vgc/qi/pair_and_trec06_testing_did/");//impact sorted binary dir
	// // const std::string intersection_index_score_1("/home/vgc/qi/pair_and_trec06_training_score_first/");//impact sorted binary dir term1
	// // const std::string intersection_index_score_1("/home/vgc/qi/pair_and_trec06_testing_score_first/");//impact sorted binary dir term1
	// // const std::string intersection_index_score_2("/home/vgc/qi/pair_and_trec06_training_score_second/");//impact sorted binary dir term2
	// // const std::string intersection_index_score_2("/home/vgc/qi/pair_and_trec06_testing_score_second/");//impact sorted binary dir term2
	// // const std::string pair_index("/home/vgc/qi/pair_index_testq_or/");//binary index pair
	// const std::string pair_lex("/home/vgc/qi/pair_lex/LM_intersection_size_14.3kk-14.44kk");//pair lex
	// // const std::string cr_results("/home/vgc/qi/11.10/cr_500_11.3_testing_10k");//cr results
	// // const std::string index_depth("/home/vgc/qi/11.10/depth_index_testing_10");//depth index
	// // const std::string pair_depth("/home/vgc/qi/11.10/pair_depth_index");//depth for pairs
	// const std::string clueweb_index("/home/vgc/qi/Info_Clueweb/InvertedIndex");//clueweb whole index
	// // const std::string Candidates_Pool("/home/vgc/qi/11.10/Candidates_test_g");//Candidates index



	// indexes
	// /data5/constantinos/trec06/
	// /data5/constantinos/trec06_sorted/
    // /data5/constantinos/trec06_layer/
	// /data5/constantinos/trec06_sorted_layer/
	// /data5/constantinos/trec06_layer2/

	// Other options for index
	// Classical Index
	// /home/sergeyn/BMW/index64trec06/
	// /data2/BMW/BACKUPindex64trec06/ --
	// home/constantinos/Desktop/BMW_next/index/index64trec06/ (external disk)
	// correct one:
	// /data2/BMW/trec06/
	// SORTED
	// /data2/BMW/index64trec06/
	// Layered
	// /data2/BMW/trec06_Layers/

	#define EXPECTED_BLOCKSIZE "64"
	const int EXPECTED_BLOCK_SIZE(atoi(EXPECTED_BLOCKSIZE));
//const std::string trecRoot("/data1/Indexing/index" + std::string(EXPECTED_BLOCKSIZE) + "/"); // index location
// for BMW const std::string trecRoot("/data2/BMW/BMW_Code_Backup/index64trec06/"); or //const std::string trecRoot("/data/sding/BMW_Code_Backup/index64trec06/");
// ("/data1/Indexing/index" + std::string(EXPECTED_BLOCKSIZE) + "/");

	// const std::string trecRawRoot("/home/constantinos/Trec06_RawIndx/");///home/sergeyn/BMW/index_result/");
	// /data/sding/TwoLevelTrec/index_result/
	// media/Book_/Indexing Backup/index_result/


//	const std::string doclenFileName(trecRoot+"doclen_file"); // /data/sding/TwoLevelTrec/index_result/doclen_file or doclen_file_sorted
	const std::string sqlPath(trecRoot+"indx.sqlite");
//	const std::string basic_table(trecRoot+"basic_table"); //trec
//	const std::string basic_table("/home/qi/Dropbox/WSDM_Index_Script/lexicon_baby"); //clueweb
	// const std::string basic_table("../../../Dropbox/WSDM_Index_Script/lexicon_baby"); //clueweb
	const std::string TERM_PATH(trecRoot+"pool");
	const std::string FLAG_PATH(trecRoot+"flag");
	const std::string MAX_PATH(trecRoot+"max");
	const std::string SCORE_PATH(trecRoot+"score");
	const std::string SIZE_PATH(trecRoot+"size");

	const std::string GOOD_TERM_SUFFIX("_good");
	const std::string BAD_TERM_SUFFIX("_bad");

	const std::string INFO_INDEX("8_1.inf");
	const std::string INDEX("8_1.dat");
	const std::string SORTED("_sorted");
	const std::string DOCUMENT_LENGTH("doclen_file");
	const std::string MERGED_BOOL_PATH("merged_bool/");
	const std::string WORD_FILE("/home/constantinos/Trec06_RawIndx/word_file");

//	const int MAXD(25205179); //trec
	const int MAXD(50220423); //clue_web
//	const int MAXTERMS(32818983); //trec
	const int MAXTERMS(86532822); //clue_web
        const float AVGD(860.917); //clue_web
	const int MAXDBITS(25);
	const int TOPK(500);
	const unsigned int STORAGE_CAPACITY(218); //2187 -- the lowest value to have 0 evictions on our 1K queries sample
	const int BS(64);

	const unsigned int bitvector_size = 1<<15;

	// Quantization
	const int Quantization(255);

	// Layering
	const int LAYER_SPLIT_POSTINGS_THRESHOLD(4096); // 65536 : Version1 (BMW): 50000 // Version2: 4096
	const float LAYER_SCORE_PARAMETER(0.25); // Version1 (BMW): 0.02  // Version2: 0.25  //0.5f // 0.25f

	const double FRAC(0.10); //no idea what those two are...
	const int S(16);

	/*11.8*/
	const ulong kClueTotalPostings = 17075485964;
//	const double percent_of_index = 0.50;	
	const uint Query_Budget_TopLayer = 1000;
	const uint Query_Budget_TermPair = 3000;
	const uint Query_Budget = 5000;
	const std::string Candidates_Pool("/home/qw376/1.28/Candidates_R");//Candidates index
	const uint Num_Doc_for_Lookups = 2000;
	const uint lookup_budget = 20000;
	const uint num_of_candidate = 2000;

};

typedef	std::vector<unsigned int> vecUInt;
typedef	std::vector<unsigned int>::iterator vecUIntItr;
typedef	std::vector<unsigned int>::const_iterator vecUIntCItr;

/*vida*/
// namespace MODEL {
//  const string term_pair_model_scores("/home/vgc/qi/QPcode/term_pair_model_scores");
//  const string top_layer_model_scores("/home/vgc/qi/QPcode/top_layer_model_scores");
// }
// namespace BLOCKS {
//   const string posting_block_boundary("/home/vgc/qi/QPcode/posting_block_boundary");  // type 0
//   const string posting_block_sizes("/home/vgc/qi/QPcode/posting_block_sizes");  // type 1
//   const string listlen_block_boundary("/home/vgc/qi/QPcode/listlen_block_boundary");  // type 2
//   const string intersection_block_boundary("/home/vgc/qi/QPcode/intersection_block_boundary");  // type 3
// }

/*dodo*/
namespace MODEL {
 const string term_pair_model_scores("/home/qw376/QPcode/bigram_model_scores");
 const string top_layer_model_scores("/home/qw376/QPcode/unigram_model_scores");
}
namespace BLOCKS {
  const string posting_block_boundary("/home/qw376/QPcode/posting_block_boundary");  // type 0
  const string posting_block_sizes("/home/qw376/QPcode/posting_block_sizes");  // type 1
  const string listlen_block_boundary("/home/qw376/QPcode/listlen_block_boundary");  // type 2
  const string intersection_block_boundary("/home/qw376/QPcode/intersection_block_boundary");  // type 3
}



#endif /* GLOBALS_H_ */
