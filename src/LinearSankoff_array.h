/*
 *LinearSankoff.h*
 header file for LinearSankoff.cpp.

 author: Sizhen Li
 created by: 08/2021
*/

#ifndef SANKOFF_H
#define SANKOFF_H

#include <string>
#include <limits>
#include <vector>
#include <set>
#include <unordered_map>

#include "LinearFold.h"
#include "LinearPartition.h"
#include "LinearSankoff.h"
#include "HMMAlign.h"
#include "check_mem.h"

using namespace std;

struct BoolState {
    bool valid;

    BoolState(): valid(false) {};
};


class SankoffParser {
public:
    float weight;
    int beam;
    float lfbeam;
    int alnbeam;
    float max_energy_diff;
    bool verbose;
    
    BeamAlign hmmalign;
    float similarity;
    float aln_viterbi;

    int num_seqs;
    int seq1_len;
    int seq2_len;
    int sum_len;
    vector<SeqObject> sequences;

    SankoffParser(float aln_weight, int beam_size, int LFbeam, int LAbeam, float energy_diff, bool is_verbose);
    ~SankoffParser(){}

    void parse(const vector<string> &seqs, bool limited, const set<pair<int, int>> &allowed_pairs, vector<pair<int, int>> &out_pairs, int num_pairs);
    void parse(const vector<string> &seqs);
    void outside(bool limited, const set<pair<int, int>> &allowed_pairs);

private:
    // state, cost of folding and alignment, three-dimentional: [s, (j1*seq1len+i1)*seq2len+i2] 
    State3 ****bestH, ****bestP, ****bestMulti;
    // unordered_map<int, State3> *bestH, *bestP, *bestMulti;
    State ****bestM, ****bestM2;
    State3 **bestC;

    // save internal loops score 
    // unordered_map<int, int>*** seq1_internal; 
    // unordered_map<int, int>*** seq2_internal; 
    short ****seq1_internal;
    short ****seq2_internal;

    vector<int> accumulated_regions;

    // bool **seq1_H_pairs, **seq1_P_pairs, **seq1_M_pairs, **seq1_M2_pairs, **seq1_Multi_pairs;
    // bool **seq2_H_pairs, **seq2_P_pairs, **seq2_M_pairs, **seq2_M2_pairs, **seq2_Multi_pairs; 

    vector<unordered_map<int, BoolState>> seq1_H_pairs, seq1_P_pairs, seq1_M_pairs, seq1_M2_pairs, seq1_Multi_pairs;
    vector<unordered_map<int, BoolState>> seq2_H_pairs, seq2_P_pairs, seq2_M_pairs, seq2_M2_pairs, seq2_Multi_pairs;

    pair<string, string> get_hmm_aln(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2);
    pair<string, string> get_hmm_aln_left(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2);
    pair<string, string> get_hmm_aln_right(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2);

    float get_hmm_score(int i1, int j1, int i2, int j2, int s1, bool allowout=false);
    float get_hmm_score_left(int i1, int j1, int i2, int j2, int s1, int s2);
    float get_hmm_score_right(int i1, int j1, int i2, int j2, int s1, int s2, bool allowout=false);

    // prepare
    void prepare(const vector<string> &seqs);

    // i1i2
    int get_key1(int j1, int j2);
    int get_key2(int j1, int j2);

    // backtrace
    tuple<string, string, string, string> get_parentheses(SeqObject& seq1, SeqObject& seq2, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_H(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_P(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_C(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_Multi(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_M2(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_M1(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    
};

#endif // SANKOFF_H