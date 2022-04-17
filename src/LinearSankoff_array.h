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

class SankoffParser {
public:
    int beam;
    int alnbeam;
    BeamAlign hmmalign;
    bool verbose;
    float weight;
    float similarity;
    float aln_viterbi;
    float max_energy_diff;

    int num_seqs;
    int seq1_len;
    int seq2_len;
    int sum_len;
    vector<SeqObject> sequences;

    // int aln_range;
    // int fold_range;
    // int valid_range;
    // vector<int> s2start;

    SankoffParser(int beam_size, int beam_size2, float aln_weight, float energy_diff, bool is_verbose=false);
    ~SankoffParser(){}

    void parse(const vector<string> &seqs, bool limited, const set<pair<int, int>> &allowed_pairs, vector<pair<int, int>> &out_pairs, int num_pairs);
    void parse(const vector<string> &seqs);
    void outside(bool limited, const set<pair<int, int>> &allowed_pairs);

private:
    // state, cost of folding and alignment, three-dimentional: [s, (j1*seq1len+i1)*seq2len+i2] 
    State3 ****bestH, ****bestP, ****bestMulti;
    // unordered_map<int, State3> *bestH, *bestP, *bestMulti;
    State ****bestM, ****bestM2;
    State **bestC;

    vector<int> accumulated_regions;

    bool **seq1_H_pairs, **seq1_P_pairs, **seq1_M_pairs, **seq1_M2_pairs, **seq1_Multi_pairs;
    bool **seq2_H_pairs, **seq2_P_pairs, **seq2_M_pairs, **seq2_M2_pairs, **seq2_Multi_pairs;

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