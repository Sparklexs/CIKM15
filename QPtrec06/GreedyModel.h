//This file is all the parts added by costas, all the code below is in other files

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

entry::entry(uint in_idx, uint in_space, float in_score) {
  idx = in_idx;
  space = in_space;
  score = in_score;
}

entry::entry() {
  idx = 0;
  space = 0;
  score = 0;
}

void entry::setEntry(uint in_idx, uint in_space, float in_score) {
  idx = in_idx;
  space = in_space;
  score = in_score;
}

bool compareScore(const entry& i, const entry& j) { return i.score > j.score; }

// methods
vector<vector<float> > LoadModel(bool is_top_layer);
vector<float> LoadModelVec(bool is_top_layer, uint &num_cols);
vector<uint> LoadBlock(const uint type);
void PrintBlockVec(const vector<uint>);
vector<float> TokenizeStringAndConvertToFloat(string line);
vector<uint> TokenizeStringAndConvertToUint(string line);
float convertStrToFloat(const string str);
uint convertStrToUint(const string str);
uint findIdxInBucket(const vector<uint> vec, const uint input);
void onlineGreedyDepthSelectionAlgorithm(const uint space_budget,
                                         const vector<vector<float> > &model,
                                         const vector<uint> &row_block_boundary,
                                         const vector<uint> &col_block_boundary,
                                         const vector<uint> &posting_block_sizes,
                                         const vector<uint> &depths,
                                         vector<uint> &cutoffs);

// method implementations

inline void onlineGreedyDepthSelectionAlgorithm(const uint space_budget,
                                                const vector<vector<float> > &model,
                                                const vector<uint> &row_block_boundary,
                                                const vector<uint> &col_block_boundary,
                                                const vector<uint> &posting_block_sizes,
                                                const vector<uint> &depths,
                                                vector<uint> &cutoffs) {
  vector<entry> entries (posting_block_sizes.size()*depths.size(), entry());

  // Compute score and add to heap.
  for (size_t i = 0; i < depths.size(); ++i) {
    size_t bucket_idx = findIdxInBucket(row_block_boundary , depths[i]);
    for (size_t j = 0; j < col_block_boundary.size(); ++j) {
      // Note the leftovers, but we avoid ifs.
      entries[i+j].setEntry(i, posting_block_sizes[j], model[i][j]);
    }
  }

  // Sort entries by score.
  sort(entries.begin(), entries.end(), compareScore); 

  // Execute the greedy algorithm. Note fix leftovers ?
  uint current_budget = 0;
  uint idx = 0;
  while (current_budget < space_budget && idx < entries.size() && entries[idx].score > 0) {
    current_budget += entries[idx].space;  // increase budget
    cutoffs[entries[idx].idx] += entries[idx].space;  // increase particular cutoff
    ++idx;
  }
}

// Given an input value, find in which col of the vec belongs to and return
// the index.
// Note that this method must be only used for listlen/intersection vectors.
inline uint findIdxInBucket(const vector<uint> vec, const uint input) {
  uint idx = 0;

  while (input > vec[idx]) {  // loop until out of range
    ++idx;
  }

  return idx;
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