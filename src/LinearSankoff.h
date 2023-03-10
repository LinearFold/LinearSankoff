/*
 *LinearSankoff.h*
 header file for LinearSankoff.cpp.

 author: Sizhen Li
 created by: 08/2021
*/

#ifndef LINEAR_SANKOFF_H
#define LINEAR_SANKOFF_H

#include <string>
#include <limits>
#include <vector>
#include <set>
#include <unordered_map>

// #include "Utils/utility.h"
// #include "Utils/utility_v.h"

#include "LinearFold.h"
#include "LinearPartition.h"
#include "HMMAlign.h"
#include "check_mem.h"

using namespace std;

#define MIN_CUBE_PRUNING_SIZE 20

// hash for integer pair
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        return (p.first << 15) + p.second; // maximum length: 32768 or Cantor pairing function
    }
};

enum Manner {
  MANNER_NONE = 0,                  // 0: empty

  MANNER_START,                     // 1

  MANNER_H_ALN,                     // 2: hairpin candidate
  MANNER_H_INS1,                    // 3: hairpin candidate
  MANNER_H_INS2,                    // 4: hairpin candidate

  MANNER_HAIRPIN_ALN,               // 5: hairpin
  MANNER_HAIRPIN_INS1,              // 6: hairpin
  MANNER_HAIRPIN_INS2,              // 7: hairpin

  MANNER_SINGLE_ALN,                // 8: single
  MANNER_SINGLE_INS1,               // 9: single
  MANNER_SINGLE_INS2,               //10: single

//   MANNER_STACKING_ALIGN_INTERNAL_LOOP1,
//   MANNER_STACKING_ALIGN_INTERNAL_LOOP2,
//   MANNER_STEM_EXTENSION1,
//   MANNER_STEM_EXTENSION2,

  MANNER_MULTI_ALN,                 // 11: multi = ..M2. [30 restriction on the left and jump on the right]
  MANNER_MULTI_INS1,                // 12: multi = ..M2. [30 restriction on the left and jump on the right]
  MANNER_MULTI_INS2,                // 13: multi = ..M2. [30 restriction on the left and jump on the right]

  MANNER_MULTI_eq_MULTI_plus_U_ALN, // 14: multi = multi + U
  MANNER_MULTI_eq_MULTI_plus_U_INS1,// 15: multi = multi + U
  MANNER_MULTI_eq_MULTI_plus_U_INS2,// 16: multi = multi + U

  MANNER_P_eq_MULTI_ALN,            // 17: P = (multi)
  MANNER_P_eq_MULTI_INS1,           // 18: P = (multi)
  MANNER_P_eq_MULTI_INS2,           // 19: P = (multi)

  MANNER_M2_eq_M_plus_P,            // 20: M2 = M + P
  MANNER_M2_eq_M_plus_P1,           // 21
  MANNER_M2_eq_M_plus_P2,           // 22
  
  MANNER_M_eq_M2,                   // 23: M = M2
  MANNER_M_eq_P,                    // 24: M = P
  MANNER_M_eq_P1,                   // 25: M = P
  MANNER_M_eq_P2,                   // 26: M = P

  MANNER_M_eq_M_plus_U_ALN,         // 27: M = M + U
  MANNER_M_eq_M_plus_U_INS1,        // 28: M = M + U
  MANNER_M_eq_M_plus_U_INS2,        // 29: M = M + U

  MANNER_C_eq_C_plus_U_ALN,         // 30: C = C + U
  MANNER_C_eq_C_plus_U_INS1,        // 31: C = C + U
  MANNER_C_eq_C_plus_U_INS2,        // 32: C = C + U

  MANNER_C_eq_C_plus_P,             // 33: C = C + P
  MANNER_C_eq_C_plus_P1,            // 34: C = C + P
  MANNER_C_eq_C_plus_P2,            // 35: C = C + P
};

struct TraceInfo {
    int split;
    struct {
        char l1;
        int l2;
    } paddings;
};

struct State {
    int j1;
    int i1; 
    int i2;
    float score;
    Manner manner;
    Manner premanner;
    // HMMManner alnmanner; // new
    int seq1foldscore, seq2foldscore;
    aln_value_type alignscore;
    // float deviation;
    // bool check;

    HMMManner startHMMstate, endHMMstate; // TODO: simplify Manner
    // HMMManner startHMMstate;
    // HMMManner endHMMstate;

    TraceInfo trace1;
    TraceInfo trace2;

    State(): j1(-1), i1(-1), i2(-1), manner(MANNER_NONE), premanner(MANNER_NONE), score(VALUE_FMIN), seq1foldscore(VALUE_MIN), seq2foldscore(VALUE_MIN), alignscore(ALN_VALUE_MIN), startHMMstate(HMMMANNER_NONE), endHMMstate(HMMMANNER_NONE) {}; //  cout << "init." << endl;
    // State(value_type s, Manner m): score(s), manner(m) {};

    void init() { // TODO
        j1 = -1;
        i1 = -1;
        i2 = -1;

        manner = MANNER_NONE;
        premanner = MANNER_NONE;

        score = VALUE_FMIN;
        seq1foldscore = VALUE_MIN;
        seq2foldscore = VALUE_MIN;
        alignscore = ALN_VALUE_MIN;
        
        startHMMstate = HMMMANNER_NONE;
        endHMMstate = HMMMANNER_NONE;
    }

    void set(int j1_, int i1_, int i2_, int seq1foldscore_, int seq2foldscore_, Manner premanner_, Manner manner_, float alignscore_, HMMManner start_hmm_m, HMMManner end_hmm_m, float walignscore_) {
        // assert (Fast_Exp(alignscore_) <= 1.0001);
        // check = true;

        j1 = j1_;
        i1 = i1_;
        i2 = i2_;

        seq1foldscore = seq1foldscore_;
        seq2foldscore = seq2foldscore_;
        alignscore = alignscore_;
        // deviation = deviation_;
        score = seq1foldscore + seq2foldscore + walignscore_; // + Fast_Exp(alignscore); 

        premanner = premanner_;
        manner = manner_;
        // alnmanner = alnmanner_;

        startHMMstate = start_hmm_m;
        endHMMstate = end_hmm_m;
    }

    void set(int j1_, int i1_, int i2_, int seq1foldscore_, int seq2foldscore_, Manner premanner_, Manner manner_, int seq1_split_, int seq2_split_, float alignscore_, HMMManner start_hmm_m, HMMManner end_hmm_m, float walignscore_) {
        // assert (Fast_Exp(alignscore_) <= 1.0001);
        // check = true;

        j1 = j1_;
        i1 = i1_;
        i2 = i2_;

        seq1foldscore = seq1foldscore_;
        seq2foldscore = seq2foldscore_;
        alignscore = alignscore_;
        // deviation = deviation_;
        score = seq1foldscore + seq2foldscore + walignscore_; // + Fast_Exp(alignscore);
        
        premanner = premanner_;
        manner = manner_; 
        // alnmanner = alnmanner_;

        trace1.split = seq1_split_;
        trace2.split = seq2_split_;

        startHMMstate = start_hmm_m;
        endHMMstate = end_hmm_m;
    }

    void set(int j1_, int i1_, int i2_, int seq1foldscore_, int seq2foldscore_, Manner premanner_, Manner manner_, char seq1_l1_, int seq1_l2_, char seq2_l1_, int seq2_l2_, float alignscore_, HMMManner start_hmm_m, HMMManner end_hmm_m, float walignscore_) {
        // assert (Fast_Exp(alignscore_) <= 1.0001);
        // check = true;

        j1 = j1_;
        i1 = i1_;
        i2 = i2_;

        seq1foldscore = seq1foldscore_;
        seq2foldscore = seq2foldscore_;
        alignscore = alignscore_;
        // deviation = deviation_;
        score = seq1foldscore + seq2foldscore + walignscore_; // + Fast_Exp(alignscore);
        
        premanner = premanner_;
        manner = manner_;
        // alnmanner = alnmanner_;
        
        trace1.paddings.l1 = seq1_l1_; 
        trace1.paddings.l2 = seq1_l2_;
        trace2.paddings.l1 = seq2_l1_; 
        trace2.paddings.l2 = seq2_l2_;

        startHMMstate = start_hmm_m;
        endHMMstate = end_hmm_m;
    }
};

struct State3 {
    // int j1, i1, i2; 

    State alnobj, ins1obj, ins2obj;

    State3(): alnobj(), ins1obj(), ins2obj() {}; // j1(-1), i1(-1), i2(-1),

    void init() {
        alnobj.init();
        ins1obj.init();
        ins2obj.init();
    }
};

struct SeqObject {
    int seq_len;
    string raw_seq;
    vector<int> nucs;
    vector<int> next_pair[NOTON];

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    SeqObject(): raw_seq(""), seq_len(-1) {};
    SeqObject(string sequence): raw_seq(sequence) {};
};

class BeamSankoffParser {
public:
    float weight;
    int beam;
    int lfbeam;
    int alnbeam;
    int alnm;
    bool use_astar;
    bool add_branch;
    float max_energy_diff;
    bool verbose;

    int seq1_mfe;
    int seq2_mfe;

    BeamAlign hmmalign;
    float similarity;
    float aln_viterbi;

    int num_seqs;
    int seq1_len;
    int seq2_len;
    int sum_len;
    vector<SeqObject> sequences;

    BeamSankoffParser(float aln_weight, int beam_size, int LFbeam, int LAbeam, int LAwidth, bool if_aster, bool if_add_branch, float energy_diff, bool is_verbose=false);
    ~BeamSankoffParser(){}

    void parse(const vector<string> &seqs, bool limited, const set<pair<int, int>> &allowed_pairs, vector<pair<int, int>> &out_pairs, int num_pairs);
    void parse(const vector<string> &seqs);
    void outside(bool limited, const set<pair<int, int>> &allowed_pairs);

private:
    // state, cost of folding and alignment, three-dimentional: [s][j1][h] 
    unordered_map<pair<int, int>, State, pair_hash>*** bestH;
    unordered_map<pair<int, int>, State, pair_hash>*** bestP;
    unordered_map<pair<int, int>, State, pair_hash>*** bestMulti;

    unordered_map<pair<int, int>, State, pair_hash>** bestM; // [s][j1]
    unordered_map<pair<int, int>, State, pair_hash>** bestM2; // [s][j1]

    vector<unordered_map<int, State3> > bestC;

    // single sequence folding outside score as hueristic
    vector<unordered_map<int, int> > seq1_out_H, seq1_out_P, seq1_out_M, seq1_out_M2, seq1_out_Multi;
    vector<unordered_map<int, int> > seq2_out_H, seq2_out_P, seq2_out_M, seq2_out_M2, seq2_out_Multi;

    vector<unordered_map<int, int> > seq1_in_P; // seq1_in_H, seq1_in_P, seq1_in_M, seq1_in_M2, seq1_in_Multi;
    vector<unordered_map<int, int> > seq2_in_P; // seq2_in_H, seq2_in_P, seq2_in_M, seq2_in_M2, seq2_in_Multi;

    // alignment backward scores
    vector<unordered_map<int, float>> aln_fw_score, ins1_fw_score, ins2_fw_score;
    vector<unordered_map<int, float>> aln_bw_score, ins1_bw_score, ins2_bw_score;

    pair<string, string> get_hmm_aln(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2);
    pair<string, string> get_hmm_aln_left(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2);
    pair<string, string> get_hmm_aln_right(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2);

    tuple<string, string, string, string> backtrace_branch_insertion(int i1, int j1, int i2, int j2);
    string backtrace_single_seq(SeqObject *seq, int i, int j, vector<unordered_map<int, int> > &seq_in_P);

    float get_hmm_score(int i1, int j1, int i2, int j2, int s1, bool allowout=false);
    float get_hmm_score_left(int i1, int j1, int i2, int j2, int s1, int s2);
    float get_hmm_score_right(int i1, int j1, int i2, int j2, int s1, int s2, bool allowout=false);

    // prepare
    void prepare(const vector<string> &seqs);

    // i1i2
    // int get_keys(int j1, int i1, int i2);

    // backtrace
    tuple<string, string, string, string> get_parentheses(SeqObject& seq1, SeqObject& seq2, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_H(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_P(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_C(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_Multi(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_M2(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);
    tuple<string, string, string, string> get_parentheses_M1(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign);

    // hzhang: sort keys in each beam to avoid randomness
    // vector<pair<int, State> > keys;
    // vector<vector<pair<float, int>>> sorted_bestM;
    // void sort_keys(unordered_map<int, State3> &map, vector<pair<int,State> > &sorted_keys);
    // void sortM(float threshold,
    //            unordered_map<int, State> &beamstep,
    //            vector<pair<float, int>> &sorted_stepM,
    //            int s,
    //            vector<unordered_map<int, int> > seq1_outside, 
    //            vector<unordered_map<int, int> > seq2_outside);

    // beam prune
    // vector<pair<float, int> > scores;
    // vector<int> invalid_pos;
    // float beam_prune(unordered_map<int, State3> &beamstep, int s, vector<unordered_map<int, int> > seq1_outside, vector<unordered_map<int, int> > seq2_outside);
    // float beam_prune(unordered_map<int, State> &beamstep, int s, vector<unordered_map<int, int> > seq1_outside, vector<unordered_map<int, int> > seq2_outside, bool if_astar);
    // // new datastructure
    vector<tuple<float, int, int, pair<int, int>> > candidates;
    bool beam_prune(unordered_map<pair<int, int>, State, pair_hash> **beststep, int s, vector<unordered_map<int, int> > &seq1_outside, vector<unordered_map<int, int> > &seq2_outside);
    float beam_prune(unordered_map<pair<int, int>, State, pair_hash> *beststep, int s, vector<unordered_map<int, int> > &seq1_outside, vector<unordered_map<int, int> > &seq2_outside);
};

////////
inline pair<int, int> hairpinScore(int i, int j, SeqObject *seq){
    int nuci = seq->nucs[i];
    int nuci1 = (i+1) < seq->seq_len ? seq->nucs[i+1] : -1;

    int jnext = seq->next_pair[nuci][j];
    // if (no_sharp_turn) TODO
    while (jnext - i < 4 && jnext != -1) jnext = seq->next_pair[nuci][jnext];

    // speed up
    if (jnext == -1 || jnext - i - 1 > HAIRPIN_MAX_LEN) 
        return make_pair(-1, 0);
  
    int nucjnext = seq->nucs[jnext];
    int nucjnext_1 = (jnext - 1) > -1 ? seq->nucs[jnext-1] : -1;


// #ifdef lv
    int tetra_hex_tri = -1;
// #ifdef SPECIAL_HP
    if (jnext-i-1 == 4) // 6:tetra
        tetra_hex_tri = seq->if_tetraloops[i];
    else if (jnext-i-1 == 6) // 8:hexa
        tetra_hex_tri = seq->if_hexaloops[i];
    else if (jnext-i-1 == 3) // 5:tri
        tetra_hex_tri = seq->if_triloops[i];

    // cout << "tetra_hex_tri: " << i << " " << j << " " << jnext << " " << tetra_hex_tri <<  endl;

// #endif
    int score = -v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
// #else
//      newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
// #endif
    // this candidate must be the best one at [j, jnext]
    // so no need to check the score
    return make_pair(jnext, score); //  / -100.0
}

inline pair<int, int> multiloopUnpairedScore(int i, int j, SeqObject *seq){
    int nuci = seq->nucs[i];
    int jnext = seq->next_pair[nuci][j];

    // char new_l1 = trace->paddings.l1;
    // int new_l2 = trace->paddings.l2 + jnext - j;
    
    // speed up
    if (jnext != -1) { // && new_l2 < SINGLE_MAX_LEN
        // 1. extend (i, j) to (i, jnext)
        int newscore;
// #ifdef lv
        newscore = -v_score_multi_unpaired(j, jnext - 1);
// #else
//         newscore = state.score + score_multi_unpaired(j, jnext - 1);
// #endif
        return make_pair(jnext, newscore); // / -100.0
    }

    return make_pair(-1, 0);
}

inline int multiloop2Pscore(int i, int j, SeqObject *seq){
    int nuci = seq->nucs[i];
    int nuci1 = (i+1) < seq->seq_len ? seq->nucs[i+1] : -1;
    int nucj = seq->nucs[j];
    int nucj_1 = seq->nucs[j-1];

    int newscore;
// #ifdef lv
    newscore = -v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq->seq_len);
// #else
//     newscore = state.score + score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
// #endif
    // update_if_better(beamP[i], newscore, MANNER_P_eq_MULTI);
    return newscore; // / -100.0;
}


inline int P2PScore(int p, int q, int i, int j, int nucp, int nucp1, int nucq_1, int nucq, int nuci_1, int nuci, int nucj, int nucj1){
    int newscore = -v_score_single(p, q, i, j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1);
//     if (p == i - 1 && q == j + 1) { // helix
// // #ifdef lv
//         newscore = -v_score_single(p, q, i, j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1);
// // #else
// //     newscore = score_helix(nucp, nucp1, nucq_1, nucq) + state.score;
// // #endif
//     } else { // single branch
// // #ifdef lv
//         newscore = -v_score_single(p, q, i, j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1);
// // #else
//     // newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
//     //     precomputed +
//     //     score_single_without_junctionB(p, q, i, j,
//     //                                    nuci_1, nuci, nucj, nucj1) +
//     //     state.score;
// // #endif
//     }
    return newscore; // / -100.0;
}

inline int branch_score(int i, int j, SeqObject *seq){
    int nuci = seq->nucs[i];
    int nuci_1 = (i-1>0)?seq->nucs[i - 1]:-1;
    int nucj = seq->nucs[j];
    int nucj1 = (j+1)<seq->seq_len?seq->nucs[j + 1]:-1;
    
    int newscore;
// #ifdef lv
    newscore = -v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq->seq_len - 1);
// #else
//  newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
// #endif

    return newscore;// / -100.0;
}

inline int external_paired_score(int k, int j, SeqObject *seq){
    int newscore;

    int nuck = seq->nucs[k];
    int nuck1 = seq->nucs[k+1];
    int nucj = seq->nucs[j];
    int nucj1 = (j+1) < seq->seq_len ? seq->nucs[j+1] : -1;

    if (k > 0)
        newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                             nucj, nucj1, seq->seq_len);
    else    
        newscore = - v_score_external_paired(0, j, -1, nuck1,
                                             nucj, nucj1, seq->seq_len);
    
    // cout << "external_paired_score: " << k << " " << j << " " << seq->seq_len << " " << newscore <<  endl; 
    
// #ifdef lv
// #else
//     newscore = score_external_paired(k+1, j, nuck, nuck1,
//                                      nucj, nucj1, seq_length) +
//     prefix_C.score + state.score;
// #endif
    return newscore;// / -100.0; 
}

inline int external_unpaired_score(int j, SeqObject *seq){
    if (j < seq->seq_len) {
// #ifdef lv
        int newscore = -v_score_external_unpaired(j, j);
        return newscore;// / -100.0; 
// #else
//     newscore = score_external_unpaired(j+1, j+1) + beamC.score;
// #endif
    }
    return VALUE_MIN;  
}


inline int multi_unpaired_score2(int i, int j, int p, int q, SeqObject *seq){
    int newscore;
// #ifdef lv
    newscore = - v_score_multi_unpaired(p+1, i-1) - v_score_multi_unpaired(j+1, q-1);
// #else
//     newscore = score_multi_unpaired(p+1, i-1) +
//         score_multi_unpaired(j+1, q-1) + state.score;
// #endif
    return newscore;// / -100.0; 
}

inline int multi_unpaired_score(int j, SeqObject *seq){
    if (j < seq->seq_len) {
// #ifdef lv
        int newscore = - v_score_multi_unpaired(j, j);
// #else
//     newscore = score_multi_unpaired(j + 1, j + 1);
// #endif
        return newscore;
    }
    return VALUE_MIN;// / -100.0; 
}

inline void update_if_better(int i1, int j1, int i2, int j2, State &state, int seq1foldscore, int seq2foldscore, Manner premanner, Manner manner, float alignscore, HMMManner start_hmm_m, HMMManner end_hmm_m,  float wegiht, bool verbose) {
    // if (alignscore <= VALUE_FMIN || alignscore <= LOG_OF_ZERO) return; // TODO
    // if (seq1foldscore <= VALUE_MIN || seq2foldscore <= VALUE_MIN) return;
    // assert (alignscore > LOG_OF_ZERO);
    // assert (j1 >= 0);
    // assert (i1 >= 0);
    // assert (i2 >= 0);
    
    if (state.score <= seq1foldscore + seq2foldscore + wegiht * alignscore) {
        if (verbose) cout << "better and update: "  << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << " " << seq1foldscore << " " << seq2foldscore << " " << alignscore  << " " << premanner << " " << manner << endl;
        state.set(j1, i1, i2, seq1foldscore, seq2foldscore, premanner, manner, alignscore, start_hmm_m, end_hmm_m, wegiht * alignscore);
    }
}

inline void update_if_better(int i1, int j1, int i2, int j2, State &state, int seq1foldscore, int seq2foldscore, Manner premanner, Manner manner, int seq1_split, int seq2_split, float alignscore, HMMManner start_hmm_m, HMMManner end_hmm_m, float wegiht, bool verbose) {
    // if (alignscore <= LOG_OF_ZERO) return; // TODO
    // if (seq1foldscore <= VALUE_MIN || seq2foldscore <= VALUE_MIN) return;
    // assert (alignscore > LOG_OF_ZERO);
    // assert (j1 >= 0);
    // assert (i1 >= 0);
    // assert (i2 >= 0);

    if (state.score <= seq1foldscore + seq2foldscore + wegiht * alignscore) {
        if (verbose) cout << "better and update: "  << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << " " << seq1foldscore << " " << seq2foldscore << " " << alignscore  << " " << premanner << " " << manner << endl;
        state.set(j1, i1, i2, seq1foldscore, seq2foldscore, premanner, manner, seq1_split, seq2_split, alignscore, start_hmm_m, end_hmm_m, wegiht * alignscore);
    }
}

inline void update_if_better(int i1, int j1, int i2, int j2, State &state, int seq1foldscore, int seq2foldscore, Manner premanner, Manner manner, char seq1_l1, int seq1_l2, char seq2_l1, int seq2_l2, float alignscore, HMMManner start_hmm_m, HMMManner end_hmm_m, float wegiht, bool verbose) {
    // if (alignscore <= VALUE_FMIN || alignscore <= LOG_OF_ZERO) return; // TODO
    // if (seq1foldscore <= VALUE_MIN || seq2foldscore <= VALUE_MIN) return;
    // assert (alignscore > LOG_OF_ZERO);
    // assert (j1 >= 0);
    // assert (i1 >= 0);
    // assert (i2 >= 0);
 
    if (state.score <= seq1foldscore + seq2foldscore + wegiht * alignscore) {
        if (verbose) cout << "better and update: "  << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << " " << seq1foldscore << " " << seq2foldscore << " " << alignscore  << " " << premanner << " " << manner << endl;
        state.set(j1, i1, i2, seq1foldscore, seq2foldscore, premanner, manner, seq1_l1, seq1_l2, seq2_l1, seq2_l2, alignscore, start_hmm_m, end_hmm_m, wegiht * alignscore);
    }
}

// temporary: save structures
inline string struc_hairpin_dots(int length){
    if (length <= 0) return "";

    char struc[length + 1];
    for (int i = 0; i < length; i ++)
        struc[i] = '.';
    struc[length] = '\0';
    return string(struc);
}
inline string struc_hairpin_pair(int length){
    if (length <= 0) return "";

    char struc[length + 3];
    struc[0] = '(';
    for (int i = 0; i < length; i ++)
        struc[i + 1] = '.';
    struc[length + 1] = ')';
    struc[length + 2] = '\0';
    return string(struc) ;
}
inline string struc_p2p_dots(string inner, int left, int right){
    left = max(left, 0);
    right = max(right, 0);

    char left_struc[left + 1];
    for (int i = 0; i < left; i ++)
        left_struc[i] = '.';
    left_struc[left] = '\0';

    char right_struc[right + 1];
    for (int i = 0; i < right; i ++)
        right_struc[i] = '.';
    right_struc[right] = '\0';

    return string(left_struc) + inner + string(right_struc);
}
inline string struc_p2p_pair(string inner, int left, int right){
    left = max(left, 0);
    right = max(right, 0);

    char left_struc[left + 2];
    left_struc[0] = '(';
    for (int i = 0; i < left; i ++)
        left_struc[i + 1] = '.';
    left_struc[left + 1] = '\0';

    char right_struc[right + 2];
    for (int i = 0; i < right; i ++)
        right_struc[i] = '.';
    right_struc[right] = ')';
    right_struc[right + 1] = '\0';

    return string(left_struc) + inner + string(right_struc);
}

#endif // LINEAR_SANKOFF_H