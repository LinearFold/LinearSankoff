/*
 *LinearFold.h*
 header file for LinearFold.cpp.

 author: Kai Zhao, Dezhong Deng, He Zhang, Liang Zhang
 edited by: 03/2021
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>

#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "data_type.h"

using namespace std;

#define MIN_CUBE_PRUNING_SIZE 20

// #ifdef lv
//   typedef int value_type;
// #define VALUE_MIN numeric_limits<int>::lowest()
// #else
//   typedef double value_type;
//   #define VALUE_MIN numeric_limits<double>::lowest()
// #endif

enum LFManner {
  LFMANNER_NONE = 0,              // 0: empty
  LFMANNER_H,                     // 1: hairpin candidate
  LFMANNER_HAIRPIN,               // 2: hairpin
  LFMANNER_SINGLE,                // 3: single
  LFMANNER_HELIX,                 // 4: helix
  LFMANNER_MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
  LFMANNER_MULTI_eq_MULTI_plus_U, // 6: multi = multi + U
  LFMANNER_P_eq_MULTI,            // 7: P = (multi)
  LFMANNER_M2_eq_M_plus_P,        // 8: M2 = M + P
  LFMANNER_M_eq_M2,               // 9: M = M2
  LFMANNER_M_eq_M_plus_U,         // 10: M = M + U
  LFMANNER_M_eq_P,                // 11: M = P
  LFMANNER_C_eq_C_plus_U,     // 12: C = C + U
  LFMANNER_C_eq_C_plus_P,     // 13: C = C + P
};

enum BestTypes {
  TYPE_C = 0,
  TYPE_H,
  TYPE_P,
  TYPE_M,
  TYPE_Multi,
  TYPE_M2,
};

// bool cmp(const tuple<value_type, int, int>& a, const tuple<value_type, int, int>& b) 
// { 
//     if (get<0>(a) != get<0>(b))
//         return (get<0>(a) < get<0>(b));

//     else
//         return  ((get<1>(a) - get<2>(a)) < (get<1>(b) - get<2>(b)));

// } 

// void window_fill(set<pair<int,int> >& window_visited, const int i, const int j, const int seq_length, const int window_size){

//     for (int ii = max(0,i-window_size); ii <= min(seq_length-1, i+window_size); ++ii){
//         for (int jj = max(0,j-window_size); jj <= min(seq_length-1, j+window_size); ++jj){
//             if (ii < jj){
//                 window_visited.insert(make_pair(ii,jj));
//             }
//         }
//     }
//     return;
// }

struct LFState {
    // double score;
    value_type score;
    LFManner manner;

    union TraceInfo {
        int split;
        struct {
            char l1;
            int l2;
        } paddings;
    };

    TraceInfo trace;

    LFState(): manner(LFMANNER_NONE), score(VALUE_MIN) {};
    LFState(value_type s, LFManner m): score(s), manner(m) {};

    void set(value_type score_, LFManner manner_) {
        score = score_; manner = manner_;
    }

    void set(value_type score_, LFManner manner_, int split_) {
        score = score_; manner = manner_; trace.split = split_;
    }

    void set(value_type score_, LFManner manner_, char l1_, int l2_) {
        score = score_; manner = manner_;
        trace.paddings.l1 = l1_; trace.paddings.l2 = l2_;
    }
};


class BeamCKYParser {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;
    bool use_constraints; // lisiz, add constraints
    bool zuker;
    int  window_size; //2 + 1 + 2 = 5 in total, 5*5 window size.
    float zuker_energy_delta;
    bool use_shape = false;
    double m = 1.8;
    double b = -0.6;
    bool is_fasta = false;

    struct DecoderResult {
        string structure;
        value_type score;
        unsigned long num_states;
        double time;
    };

    BeamCKYParser(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false,
                  bool is_constraints=false,
                  bool zuker_subopt=false,
                  float zuker_energy_delta=5.0,
                  string shape_file_path="",
                  bool is_fasta=false); // lisiz, add constraints

    vector<unordered_map<int, LFState>> bestH, bestP, bestM2, bestMulti, bestM;
    //Zuker subopt
    vector<unordered_map<int, LFState>> bestH_beta, bestP_beta, bestM2_beta, bestMulti_beta, bestM_beta;

    int parse(const string& seq, vector<int>* cons);
                
    void outside(vector<int> next_pair[]); //for zuker subopt

private:
    void get_parentheses(char* result, const string& seq);

    // pair<string, string> get_parentheses_outside_real_backtrace(int i, int j, State& state_beta, map<tuple<BestTypes, int, int>, pair<string, string> >& global_visited_outside, map<tuple<BestTypes, int, int>, string>& global_visited_inside, set<pair<int,int> >& window_visited);
    // string get_parentheses_inside_real_backtrace(int i, int j, State& state, map<tuple<BestTypes, int, int>, string>& global_visited_inside, set<pair<int,int> >& window_visited);

    unsigned seq_length;

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    // same as bestM, but ordered
    vector<vector<pair<value_type, int>>> sorted_bestM;

    // hzhang: sort keys in each beam to avoid randomness
    vector<pair<int, LFState>> keys;

    // hzhang: sort keys in each beam to avoid randomness
    void sort_keys(unordered_map<int, LFState> &map, vector<pair<int,LFState>> &sorted_keys);

    void sortM(value_type threshold,
               unordered_map<int, LFState> &beamstep,
               vector<pair<value_type, int>>& sorted_stepM);
    vector<LFState> bestC;
    vector<int> nucs;

    //Zuker subopt
    vector<LFState> bestC_beta;
    
    // SHAPE
    vector<double> SHAPE_data;
    vector<int> pseudo_energy_stack;


    // lisiz: constraints
    vector<int> allow_unpaired_position;
    vector<int> allow_unpaired_range;
    bool allow_paired(int i, int j, vector<int>* cons, char nuci, char nucj);

    void prepare(unsigned len);

    void update_if_better(LFState &state, value_type newscore, LFManner manner) {
      if (state.score < newscore)
            state.set(newscore, manner);
    };

    void update_if_better(LFState &state, value_type newscore, LFManner manner, int split) {
        if (state.score < newscore || state.manner == LFMANNER_NONE)
            state.set(newscore, manner, split);
    };

    void update_if_better(LFState &state, value_type newscore, LFManner manner, char l1, int l2) {
        if (state.score < newscore || state.manner == LFMANNER_NONE)
            state.set(newscore, manner, l1, l2);
    };

    value_type beam_prune(unordered_map<int, LFState>& beamstep);

    // vector to store the scores at each beam temporarily for beam pruning
    vector<pair<value_type, int>> scores;
};

#endif //FASTCKY_BEAMCKYPAR_H
