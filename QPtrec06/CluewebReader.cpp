#include "pfor.h"
#include "CluewebReader.h"
#include "math.h"
#include "globals.h"

using namespace std;

inline FILE* safeFopen(const std::string& path) {
	FILE* res = fopen64(path.c_str(), "r");
	if(!res)
		CERR << "failed to open a file at:" << path  << EFATAL;
	return res;
}
CluewebReader::CluewebReader(int _docn){

//	findex = safeFopen(index_path); // open 8_1.dat file
//	finf = safeFopen( inf_path); // open 8_1.inf file
//	flex = safeFopen(word_file.c_str());  // open word_file
//	fdoclength = safeFopen(doclength_path.c_str()); // open doc_len file
//
//	int infn;
//	fread( &infn, sizeof(int), 1, finf); // first int contains the # of terms in the entire collection
//	COUT1<<"Totally there are "<<infn<<" terms"<<Log::endl;
//	inf_buffer = new unsigned int[ 4 * infn]; // read
//	fread( inf_buffer, sizeof(int), 4 * infn, finf);
//	//now we have all our sizes in an inf_buffer. Let's partial sum them
//	inf_prefix_sum = new size_t[infn+1];
//	inf_prefix_sum[0] = 0;
//	for(size_t i=0; i<4 * infn; i+=4)
//		inf_prefix_sum[i/4 + 1] = (2*inf_buffer[i+1]*sizeof(unsigned int))+inf_prefix_sum[i/4];

	// inf_buffer[0] = term_id and inf_buffer[1] == padded postings
		docn = _docn; // CONSTS::MAXD Total # of documents in collection
		doclen = new unsigned int[docn + 256];

		termid_ = new unsigned int[CONSTS::number_of_terms];
		listLen_ = new unsigned long[CONSTS::number_of_terms];
		offset_ = new unsigned long[CONSTS::number_of_terms];
		listLenBytes_ = new unsigned long[CONSTS::number_of_terms];

		load_doclength();
		cout<<"doclength is loaded"<<endl;

		// string input_source( "/data/qw376/experiments/lexicon_baby_712" );
		ifstream doc_lexikon_stream;
		// doc_lexikon_stream.open(input_source.c_str());
		doc_lexikon_stream.open(CONSTS::basic_table.c_str());

		string lexicon_line;
		string term;
		string termid_s;
		string listLen_s;
		string offset_s;
		string listLenBytes_s;
		int termid;
		long listLen;
		long offset;
		long listLenBytes;
		int linenum = 0;



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
			term_.push_back(term);
			term_map[term] = linenum;
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
		  	 termid_[linenum] = termid;
		  	 //termID

		  	 //listLen, come after the 2nd space
		  	 start1 = itr1+1;
		  	 itr1++;
		  	 while (!isspace(*itr1)) {
		  	 ++itr1;
		  	 }
		  	 listLen_s = string(start1, itr1);
		  	 listLen = atol(listLen_s.c_str());
			 // cout<<"listLen: "<<listLen<<" ";
			 listLen_[linenum] = listLen;
		  	 //listLen

			//offset, come after the 3th space
		  	start1 = itr1+1;
		  	itr1++;
		  	while (!isspace(*itr1)) {
		  	++itr1;
		  	}
		  	offset_s = string(start1, itr1);
		  	offset = atol(offset_s.c_str());
			// cout<<"offset: "<<offset<<" ";
			offset_[linenum] = offset;
		  	//offset

		  	//listLenBytes, come after the 4th space
		  	start1 = itr1+1;
		  	itr1++;
		  	while (!isspace(*itr1)) {
		  	   ++itr1;
		  	}
		  	listLenBytes_s = string(start1, itr1);
		  	listLenBytes = atol(listLenBytes_s.c_str());
		  	// cout<<"listLenBytes: "<<listLenBytes<<endl;
		  	listLenBytes_[linenum] = listLenBytes;
			//listLenBytes
			linenum++;
		}
		num_terms = linenum;
		cout<<"number of terms: "<<num_terms<<endl;
		doc_lexikon_stream.close();

//		for(int i=0; i<10; i++){
//			cout<<term_[i]<<" ";
//			cout<<termid_[i]<<" ";
//			cout<<listLen_[i]<<" ";
//			cout<<offset_[i]<<" ";
//			cout<<listLenBytes_[i]<<endl;
//		}

//		for (map<string, int>::iterator it = term_map.begin() ; it != term_map.end(); ++it){
//		        	cout<<it->first<<" "<<it->second<<endl;
//		}
}

CluewebReader::~CluewebReader(void)
{
//	fclose(findex);
//	fclose(finf);
//	fclose(flex);
//	fclose(fdoclength);

	delete[] doclen;

	delete[] termid_;
	delete[] listLen_;
	delete[] offset_;
	delete[] listLenBytes_;
}

void CluewebReader::load_doclength(){

	// string input_source( "/data/constantinos/Index/documentLengths" );
	ifstream inputstream;
	inputstream.open(CONSTS::doclenFileName.c_str());
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

	    doclen[doc_id] = doc_length;
	}

	inputstream.close();

	for(int i = 0; i< 256; i++)
		doclen[ i + docn ] = 0;
}

void  CluewebReader::loadRawListIn(RawIndexList& tList) {

//	size_t wid = tList.termId;
//	assert(wid == inf_buffer[ 4 * (wid-1) ]);
//	size_t listsize = inf_buffer[4*(wid-1) + 1];
//	fseek(findex,inf_prefix_sum[wid-1],SEEK_SET); //jump to the right position
//	if( (listsize*2)!= fread( hold_buffer, sizeof(int), listsize*2, findex ))
//				CERR << "can not read " << EFATAL; //<<listsize*2
//
//	tList.doc_ids.reserve(listsize+CONSTS::BS);
//	tList.freq_s.reserve(listsize+CONSTS::BS);
//
//	for(int i = 0 ;i < listsize; i++){
//		tList.doc_ids.push_back(hold_buffer[2*i] - 1);
////		cout<<"doc_id: "<<hold_buffer[2*i] - 1<<" ";
//		tList.freq_s.push_back(hold_buffer[2*i + 1]);
////		cout<<"freq "<<hold_buffer[2*i + 1]<<endl;
//	}
//	tList.lengthOfList = listsize;

	size_t tid = tList.termId;
	unsigned long listLen = listLen_[tid];
	unsigned long listLenBytes = listLenBytes_[tid];
	unsigned long offset = offset_[tid];

	compressed_list = new unsigned char[listLenBytes];
	uncompressed_list = new unsigned int[listLen*2];

	// const char* InvertedIndex = "/data/constantinos/InvertedIndex";
	FILE* F_InvertedIndex = fopen(CONSTS::clueweb_index.c_str(),"r");
	if( F_InvertedIndex == NULL ) cout << "Problem! The file: " << CONSTS::clueweb_index << " could not be opened!" << endl;
	//  	  else cout << "Loading from file: " << InvertedIndex << endl;
	 fseek(F_InvertedIndex, offset, SEEK_SET);
	 fread(compressed_list, 1, listLenBytes, F_InvertedIndex);


	 decompression_vbytes(compressed_list, uncompressed_list, listLen*2);

	 tList.doc_ids.reserve(listLen+CONSTS::BS);
	 tList.freq_s.reserve(listLen+CONSTS::BS);

	 for(unsigned int i=0; i<listLen; i++){
//		 cout<<"doc_id: "<<tList.doc_ids.back()+uncompressed_list[2*i]<<" ";
		 if(tList.doc_ids.size()!=0)
		 tList.doc_ids.push_back(tList.doc_ids.back()+uncompressed_list[2*i]);
		 else
	     tList.doc_ids.push_back(uncompressed_list[2*i]);
//		 cout<<"freq: "<<uncompressed_list[2*i+1]<<endl;
		 tList.freq_s.push_back(uncompressed_list[2*i+1]);

	 }

	 tList.lengthOfList = listLen;
	 fclose (F_InvertedIndex);
	 free (compressed_list);
	 free (uncompressed_list);

	// cout<<"tid: "<<tid<<" listLen: "<<listLen<<endl;
}

RawIndexList CluewebReader::getRawList(const std::string& term, size_t wid){
	BasicList blTerm(term,wid);
	RawIndexList riTerm (blTerm);
	loadRawListIn(riTerm); //now it has the docids and freqs
	//debug check what it prints and doublecheck, should be the same as in the unpadded postings file
	//COUT1<<term<<": unpadded length: "<<riTerm.lengthOfList<<Log::endl;
	//unsigned int unpadded_list_length = riTerm.lengthOfList;

	// add dummy entries at the end to next multiple of CHUNK_SIZE
    for (size_t i = 0; (i == 0) || ((riTerm.lengthOfList&(CONSTS::BS-1)) != 0); ++i)     {
    	riTerm.doc_ids.push_back(CONSTS::MAXD + i);
    	riTerm.freq_s.push_back(1);
    	++riTerm.lengthOfList;
    }

    riTerm.rankWithBM25(doclen); // basic version was this but changed with the following one so we compute the correct max score
    //riTerm.rankWithBM25(doclen, unpadded_list_length);
    // Note: the list is ranked and the scores are in riTerm.scores, maxScore for list is set as well
	//COUT1<<term<<":"<<riTerm.maxScoreOfList<<Log::endl;

	return riTerm;
}

// Correct - unpadded score computation of BM25
// Usage: it loads to a RawIndexList structure (vectors of scores, docids, freqs)
// Input: term, term_id
// Output: RawIndexList structure
RawIndexList CluewebReader::load_raw_list(const std::string& term, size_t wid) {
	// arguments: term, term_id
	BasicList basic_list_Term(term, wid);
	RawIndexList Raw_list(basic_list_Term);

	// Get docids, freqs for specific term
	loadRawListIn(Raw_list);

	// Rank all docids
	Raw_list.rankWithBM25(doclen);

	// set unpadded list length
	Raw_list.unpadded_list_length = Raw_list.lengthOfList;

	// add dummy entries at the end to next multiple of CHUNK_SIZE
	// Docids = MAXD + 1, freqs = 1, scores = 0.0f
    for (size_t i = 0; (i == 0) || ((Raw_list.lengthOfList&(CONSTS::BS-1)) != 0); ++i)     {
    	Raw_list.doc_ids.push_back(CONSTS::MAXD + i);
    	Raw_list.freq_s.push_back(1);
    	Raw_list.scores.push_back(0.0f);
    	++Raw_list.lengthOfList;
    }
	return Raw_list;
}
// Usage: Load raw index rank all postings and store maxscore of list and unpadded list length
// Input: map of term termid, path, offset, limit
// Output: file with pairs of maxscore and list lengths
void CluewebReader::dump_To_File_Pairs_of_Maxscore_and_Unpadded_List_Lengths(const stringIntVectorsPair& tmap, const std::string& path, int offset, int limit) {
	// folder + file to store maxscore, list length pairs
	const std::string PATH = path + CONSTS::Entire_Index_Maxscore_Unpadded_List_length;

	// open handler to write basic_table file
	FILE *maxscore_list_length_handler = fopen(PATH.c_str(),"a");

	COUT2 << "Open file for writing: " << PATH << Log::endl;
	float maximum_score_over_index = 0.0f;
	int counter = 0;

	// for all index
	size_t goTo = std::min(tmap.first.size(), size_t(limit));
	for(size_t it = offset; it < goTo; ++it) {
		std::string term =  tmap.first[it]; // term
		size_t wid =  tmap.second[it]; // term_id

		// original term declaration and filling of docids, freqs vectors
		RawIndexList original_Term (term);
		original_Term = load_raw_list(term, wid);
		++counter;

		//std::cout << term << "\t" << original_Term.maxScoreOfList<< "\t" << original_Term.unpadded_list_length << std::endl;

		// find the maximum score in the entire index
		if ( Fcompare(original_Term.maxScoreOfList, maximum_score_over_index) == 1 )
			maximum_score_over_index = original_Term.maxScoreOfList;

		// write to files
		fprintf(maxscore_list_length_handler,"%f\t%u\n", original_Term.maxScoreOfList, original_Term.unpadded_list_length);
	}
	fclose(maxscore_list_length_handler);

	// print reports
	std::cout << "Total terms processed: " << counter << std::endl;
	std::cout << "Maximum score over the entire index: " << maximum_score_over_index << std::endl;
}

/*
 *  This is used to prepare the compressed list
 *
 *  Input:
 *  1. the file which contains the length for each document
 */
int CluewebReader::prepare_list(const stringIntVectorsPair& tmap, const std::string& rootPath, int offset, int limit){
	const std::string BASIC_PATH = rootPath+"basic_table";

	FILE *fbasic = fopen(BASIC_PATH.c_str(),"a");

	//COUT << "Loaded mappings: " << tmap.first.size() << std::endl;

//	sqlLoader sql;
//    const char zSql[] = "CREATE  TABLE  IF NOT EXISTS 'terms'('term' TEXT PRIMARY KEY  NOT NULL  UNIQUE, 'flagB' BLOB, 'maxB' BLOB, 'minB' BLOB, 'scoreB' BLOB, 'sizeB' BLOB)";
//    sqlite3_exec(db, zSql, 0, 0, zerr);

	size_t goTo = std::min(tmap.first.size(), size_t(limit));
	for(size_t it = offset; it < goTo; ++it)	{
		std::string term =  tmap.first[it]; // term
		size_t wid =  tmap.second[it]; // termid
		COUT2 << "compressing term: " << term << Log::endl;

		RawIndexList original_Term (term);
		original_Term = load_raw_list(term, wid);
		// get the compressed entire list for the original term
		CompressedList original_cmp(original_Term);
		//CompressedList cmp(load_raw_list(word,wid)); // changed from getRawList

		fprintf(fbasic,"%s\t%u\n", term.c_str(), original_Term.unpadded_list_length); // change from cmp.lengthOfList
		original_cmp.serializeToFS(rootPath);
		original_cmp.serializeToDb(rootPath, sql);
		//cmp.freeMem();
	}

	fclose(fbasic);
	return 0;
}

// ADDED
int CluewebReader::prepare_list_Layers(const stringIntVectorsPair& tmap, const std::string& rootPath, int offset, int limit){
	// print parameter for layering
	COUT3 << "Layer creation starting with Threshold Parameter = " << CONSTS::LAYER_SCORE_PARAMETER << " and Posting split parameter = " << CONSTS::LAYER_SPLIT_POSTINGS_THRESHOLD << Log::endl;

	// folder + file to store term list length pairs
	const std::string BASIC_PATH = rootPath+"basic_table";

	// open handler to write basic_table file
	FILE *fbasic = fopen(BASIC_PATH.c_str(),"a");

	// open handler to write layered terms mapping
	FILE *term_mapping = fopen((rootPath+"10000q_terms_mapping"+"_layered").c_str(),"a");

	// statistics
	int layered_terms = 0;
	int no_layered_terms = 0;
	int total = 0;
	unsigned long long good_postings = 0;
	unsigned long long bad_postings = 0;
	unsigned long long single_postings = 0;

	// for all lexicon
	size_t goTo = std::min(tmap.first.size(), size_t(limit));
	for(size_t it = offset; it < goTo; ++it) {
		std::string term =  tmap.first[it]; // term
		size_t wid =  tmap.second[it]; // term_id

		// original term declaration and filling of docids, freqs vectors
		RawIndexList original_Term (term);
		original_Term = load_raw_list(term, wid);

		// if list length of original term > layer threshold, create layers, otherwise not
		if ( original_Term.unpadded_list_length > CONSTS::LAYER_SPLIT_POSTINGS_THRESHOLD ) {
			// create good and bad terms
			BasicList gTerm(term+CONSTS::GOOD_TERM_SUFFIX, wid);
			BasicList bTerm(term+CONSTS::BAD_TERM_SUFFIX, wid);

			RawIndexList good_Term (gTerm);
			RawIndexList bad_Term (bTerm);

			// debug
			COUT2 << "Creating layers for term: " << term << Log::endl;

			// given the original Term with vectors of docids, freqs and scores, construct the good and the bad terms
			build_layered_index(original_Term, good_Term, bad_Term);

			// get the compressed entire list the good and bad term
			CompressedList good_cmp(good_Term);
			CompressedList bad_cmp(bad_Term);

			// write in basic_table good and bad term
			//fprintf(fbasic,"%s\t%d\n",(term+CONSTS::GOOD_TERM_SUFFIX).c_str(), good_cmp.lengthOfList );
			//fprintf(fbasic,"%s\t%d\n",(term+CONSTS::BAD_TERM_SUFFIX).c_str(), bad_cmp.lengthOfList );

			// store in basic table: term, original_unpadded list length, unpadded list length of list
			fprintf(fbasic,"%s\t%u\t%u\n",(term+CONSTS::GOOD_TERM_SUFFIX).c_str(), good_Term.original_unpadded_list_length, good_Term.unpadded_list_length);
			fprintf(fbasic,"%s\t%u\t%u\n",(term+CONSTS::BAD_TERM_SUFFIX).c_str(), bad_Term.original_unpadded_list_length, bad_Term.unpadded_list_length);

			// write in terms mapping good and bad term mappings
			fprintf(term_mapping,"%s\t%lu\n",(term+CONSTS::GOOD_TERM_SUFFIX).c_str(), wid );
			fprintf(term_mapping,"%s\t%lu\n",(term+CONSTS::BAD_TERM_SUFFIX).c_str(), wid );

			// write to filesystem the good and bad term
			good_cmp.serializeToFS( rootPath );
			bad_cmp.serializeToFS( rootPath );

			// write to sql
			good_cmp.serializeToDb(rootPath, sql);
			bad_cmp.serializeToDb(rootPath, sql);

			// stats counter
			++layered_terms;
			good_postings += good_Term.unpadded_list_length;
			bad_postings += bad_Term.unpadded_list_length;
		} else {
			// get the compressed entire list for the original term
			CompressedList original_cmp(original_Term);

			// write in basic_table the original term's unpadded list length
            fprintf(fbasic,"%s\t%u\t%u\n",term.c_str(), original_Term.unpadded_list_length, original_Term.unpadded_list_length);

			// write in terms mapping the original term mapping
            fprintf(term_mapping,"%s\t%lu\n",term.c_str(), wid );

			// write to filesystem the original term
			original_cmp.serializeToFS( rootPath );

			// future use for sql
			original_cmp.serializeToDb(rootPath, sql);

			// stats counter
			++no_layered_terms;
			single_postings += original_Term.unpadded_list_length;
		}
		//stats counter
		++total;

		if (layered_terms%1000==0)
			std::cout << "good: " << good_postings << "\t bad: " << bad_postings << std::endl;
		if (no_layered_terms%1000==0)
			std::cout << "single: " << single_postings << std::endl;
	}

	// print stats
	COUT3 << "# layered terms: " << layered_terms << " , # no layered terms: " << no_layered_terms << " , total terms: " << total << Log::endl;
	COUT3 << "# good postings: " << good_postings << " , # bad postings: " << bad_postings << " , single postings: " << single_postings << Log::endl;

	fclose(fbasic);
	fclose(term_mapping);

	return 0;
}


// Usage: Creating Layered index
// Input: Given the original rawindexlist construct the good and the bad term
void CluewebReader::build_layered_index(RawIndexList& original_Term, RawIndexList& good_Term, RawIndexList& bad_Term) {
	// obtain split threshold
	//float split_threshold = get_split_threshold(original_Term.scores); // BMW Version
	float split_threshold = CONSTS::LAYER_SCORE_PARAMETER*original_Term.maxScoreOfList; // Parameters Version

	// for all docids of the original list
	for (int i=0; i<original_Term.unpadded_list_length; i++) {
		// split good and bad terms based on the docid's score
		if ( Fcompare(original_Term.scores.at(i), split_threshold) < 0 ) {
			bad_Term.doc_ids.push_back(original_Term.doc_ids.at(i));
			bad_Term.freq_s.push_back(original_Term.freq_s.at(i));
			bad_Term.scores.push_back(original_Term.scores.at(i));
		} else { // add to good list
			good_Term.doc_ids.push_back(original_Term.doc_ids.at(i));
			good_Term.freq_s.push_back(original_Term.freq_s.at(i));
			good_Term.scores.push_back(original_Term.scores.at(i));
		}
	}

	// fill the max score of the list
	good_Term.maxScoreOfList = original_Term.maxScoreOfList; // it must be in the good layer
	bad_Term.maxScoreOfList = *(std::max_element(bad_Term.scores.begin(), bad_Term.scores.end()));

	// fill unpadded length of lists
	good_Term.unpadded_list_length = good_Term.doc_ids.size();
	bad_Term.unpadded_list_length = bad_Term.doc_ids.size();

	// store original's term unpadded list length in both layers (it's needed for computing the correct pre)
	// we can later store only the deltas TODO
	good_Term.original_unpadded_list_length = original_Term.unpadded_list_length;
	bad_Term.original_unpadded_list_length = original_Term.unpadded_list_length;

	// store the padded list length of lists (we increase it below, when we pad)
	good_Term.lengthOfList = good_Term.unpadded_list_length;
	bad_Term.lengthOfList = bad_Term.unpadded_list_length;

	// add dummy entries at the end to next multiple of CHUNK_SIZE
    for (size_t i = 0; (i == 0) || ((good_Term.lengthOfList&(CONSTS::BS-1)) != 0); ++i) {
    	good_Term.doc_ids.push_back(CONSTS::MAXD + i);
    	good_Term.freq_s.push_back(1);
    	good_Term.scores.push_back(0.0f);
    	++good_Term.lengthOfList;
    }

	// add dummy entries at the end to next multiple of CHUNK_SIZE
    for (size_t i = 0; (i == 0) || ((bad_Term.lengthOfList&(CONSTS::BS-1)) != 0); ++i) {
    	bad_Term.doc_ids.push_back(CONSTS::MAXD + i);
    	bad_Term.freq_s.push_back(1);
    	bad_Term.scores.push_back(0.0f);
    	++bad_Term.lengthOfList;
    }
}

// Usage: When we want to split layers as in BMW (topk scored documents go to the good list)
float CluewebReader::get_split_threshold(std::vector<float>& scores) {
	std::vector<float> local_scores (scores.size(), 0.0f);
	local_scores = scores;
	std::sort(local_scores.begin(), local_scores.end());
	int split_position = (int) scores.size()*CONSTS::LAYER_SCORE_PARAMETER;
	return local_scores[scores.size() - split_position];
}

void CluewebReader::load_baby_index(){
	 SqlProxy sql(CONSTS::trec_output+"indx.sqlite");
	 for (size_t i=0; i<num_terms; i++){
	 			   const string term = term_[i];
	 			   size_t termid = termid_[i];
	 			   RawIndexList Rlist = load_raw_list(term_[i],i);
	 			   CompressedList Clist(Rlist);
	 			   Clist.serializeToDb(CONSTS::trec_output, sql);
	 }
}

int CluewebReader::decompression_vbytes(unsigned char* input, unsigned int* output, int size){
    unsigned char* curr_byte = input;
    unsigned int n;
    for (int i = 0; i < size; ++i) {
	    unsigned char b = *curr_byte;
	    n = b & 0x7F;
//        cout<<"The first byte: "<<n<<endl;
//        print_binary(n);
//        cout<<endl;

	    while((b & 0x80) !=0){
	    	n = n << 7;
	    	++curr_byte;
	        b = *curr_byte;
	        n |= (b & 0x7F);
//	        cout<<"The following byte: "<<n<<endl;
//	        print_binary(n);
//	        cout<<endl;
	    }
    ++curr_byte;
    output[i] = n;
  }

  int num_bytes_consumed = (curr_byte - input);
  return (num_bytes_consumed >> 2) + ((num_bytes_consumed & 3) != 0 ? 1 : 0);

}

// Print unsigned integer in binary format with spaces separating each byte.
void CluewebReader::print_binary(unsigned int num) {
	 int arr[32];
	  int i = 0;
	  while (i++ < 32 || num / 2 != 0) {
	    arr[i - 1] = num % 2;
	    num /= 2;
	  }

	  for (i = 31; i >= 0; i--) {
	    printf("%d", arr[i]);
	    if (i % 8 == 0)
	      printf(" ");
	  }
	  printf("\n");
}


CluewebReader* CluewebFactory() {
//	const std::string mergeBoolPath = index_res + "merged_bool/";
//	const std::string mergeBoolPath = index_res;
//	const std::string dat = mergeBoolPath+"8_1.dat";  // 8_1.dat_sorted // 8_1.dat // TRANSFER TO globals.h
//	const std::string inf = mergeBoolPath+"8_1.inf";
//	const std::string doclen = index_res+"doclen_file";  // doclen_file_sorted // doclen_file
//	const std::string word = index_res+"word_file";

	CluewebReader* reader = new CluewebReader(CONSTS::MAXD);
	return reader;
}

//void CluewebFactory(CluewebReader& reader, const std::string& resultRootPath, int offset, int limit, stringIntVectorsPair& tmap) 	{
//	if(tmap.first.empty()) {
//		FILE *ftrec = fopen(CONSTS::termsMapping.c_str(),"r"); // based on layering
//		if(ftrec == NULL)
//			CERR << "trec file is empty" << EFATAL;
//		tmap = getTrecWordsMappingList(ftrec);
//		fclose(ftrec);
//	}
//
//	//reader.prepare_list(tmap, resultRootPath, offset, limit);
//
//	// ADDON LAYERED version
//	//reader.prepare_list_Layers(tmap, resultRootPath, offset, limit);
//
//	// additional - add this line inside if case in the above code
//	//std::string path = CONSTS::trecRawRoot + CONSTS::Index_term_mapping;
//	//reader.dump_To_File_Pairs_of_Maxscore_and_Unpadded_List_Lengths(tmap, resultRootPath, offset, CONSTS::MAXTERMS);
//	// if you want to use it, change the path that you store the file, because now its in globals.h and has prefix ../QueryLog/
//}
