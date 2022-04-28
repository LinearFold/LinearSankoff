/*
 *LinearSankoff.cpp*
 The main code for LinearSankoff: Linear-Time simultaneous folding and alignment of two sequences 

 author: Sizhen Li
 created by: 09/2021
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <sys/time.h>

#include "LinearSankoff.h"
// #include "HMMAlign.h"

// #include <sparsehash/dense_hash_map>
// using google::dense_hash_map;      // namespace where class lives by default

// #define dynalign

using namespace std;


#ifdef is_cube_pruning
// cube pruning
bool comparefunc(pair<int,State> a, pair<int,State> b) {
    return a.first > b.first;
}
void BeamSankoffParser::sort_keys(unordered_map<int, State3> &map, vector<pair<int,State> > &sorted_keys) {
    sorted_keys.clear();
    for(auto &item : map) {
        sorted_keys.push_back(make_pair(item.first, item.second.alnobj));
    }
    sort(sorted_keys.begin(), sorted_keys.end(), comparefunc);    
}
void BeamSankoffParser::sortM(float threshold,
                              unordered_map<int, State> &beamstep,
                              vector<pair<float, int>> &sorted_stepM,
                              int s,
                              vector<unordered_map<int, int> > seq1_outside, 
                              vector<unordered_map<int, int> > seq2_outside) {
    sorted_stepM.clear();
    if (threshold == VALUE_FMIN) {
        // no beam pruning before, so scores vector not usable
        for (auto &item : beamstep) {
            State &cand = item.second;
            int i1 = cand.i1;
            int i2 = cand.i2;
            int j1 = cand.j1;
            int j2 = s - j1;

            int k1 = i1 - 1;
            int k2 = i2 - 1;

            // get new score
            float newscore = VALUE_MIN;
            // alignment score
            float alignscore = xlog_mul(xlog_mul(bestC[k1+k2][k1].alignscore, cand.alignscore), hmmalign.trans_probs[bestC[k1+k2][k1].endHMMstate-1][2]);
            float backward_score = hmmalign.bestALN[s][j2].beta;
            alignscore = xlog_mul(alignscore, backward_score);
            // alignscore = xlog_div(alignscore, aln_viterbi);

            // folding score
            int foldingscore = cand.seq1foldscore + cand.seq2foldscore;
            // folding heuristic
            int seq1_out = VALUE_MIN;
            int seq2_out = VALUE_MIN;
            if ((seq1_outside[j1].find(i1) != seq1_outside[j1].end())) {
                seq1_out = seq1_outside[j1][i1];
            }
            if ((seq2_outside[j2].find(i2) != seq2_outside[j2].end())) {
                seq2_out = seq2_outside[j2][i2];
            }

            if (seq1_out == VALUE_MIN || seq2_out == VALUE_MIN) {
                newscore = VALUE_FMIN;
            } else {
                foldingscore += seq1_out + seq2_out;
                newscore = foldingscore + weight * alignscore;
            }

            // float newscore = (k >= 0 ? bestC[k][i1-1].score : 0) + cand.score; // old
            // cout << "sorting M " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << newscore << endl;
            if (newscore > VALUE_FMIN)
                sorted_stepM.push_back(make_pair(newscore, item.first));
        }
    } else {
        for (auto &p : scores) {
            if (p.first >= threshold) sorted_stepM.push_back(p);
        }
    }
    sort(sorted_stepM.begin(), sorted_stepM.end(), greater<pair<float, int>>());
}
#endif

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_H(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    if ( !stk.empty() ) {
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) printf("get_parentheses_H manner at %d, %d, %d, %d: manner %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.manner, state.alignscore, state.seq1foldscore, state.seq2foldscore);

        switch (state.manner) {
            case MANNER_H_ALN:
                {
                    string sub_struc1 = struc_hairpin_pair(j1 - i1 -1);
                    string sub_struc2 = struc_hairpin_pair(j2 - i2 -1);
                    pair<string, string> aln_ret = get_hmm_aln(i1, j1, i2, j2, MANNER_ALN, MANNER_ALN);

                    return make_tuple(sub_struc1, sub_struc2, get<0>(aln_ret), get<1>(aln_ret));
                }
                break;
            case MANNER_H_INS1:
                {
                    string sub_struc1 = struc_hairpin_pair(j1 - i1 -1);
                    string sub_struc2 = struc_hairpin_dots(j2 - i2 -1);
                    pair<string, string> aln_ret = get_hmm_aln(i1, j1, i2+1, j2-1, MANNER_INS1, MANNER_INS1);
                    return make_tuple(sub_struc1, sub_struc2, get<0>(aln_ret), get<1>(aln_ret));
                }
                break;
            case MANNER_H_INS2:
                {
                    string sub_struc1 = struc_hairpin_dots(j1 - i1 -1);
                    string sub_struc2 = struc_hairpin_pair(j2 - i2 -1);
                    pair<string, string> aln_ret = get_hmm_aln(i1+1, j1-1, i2, j2, MANNER_INS2, MANNER_INS2);
                    return make_tuple(sub_struc1, sub_struc2, get<0>(aln_ret), get<1>(aln_ret));
                }         
                break;
            
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                assert(false);
        }

    }
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_Multi(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    // Multi -> M2
    if ( !stk.empty() ) {
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) printf("get_parentheses_Multi manner at %d, %d, %d, %d: manner %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.manner, state.alignscore, state.seq1foldscore, state.seq2foldscore);

        int seq1_l1 = state.trace1.paddings.l1;
        int seq1_l2 = state.trace1.paddings.l2;

        int seq2_l1 = state.trace2.paddings.l1;
        int seq2_l2 = state.trace2.paddings.l2;

        int p1 = i1 + seq1_l1; // i <= p < q <= j
        int q1 = j1 - seq1_l2;

        int p2 = i2 + seq2_l1;
        int q2 = j2 - seq2_l2;

        // stk.push(make_tuple(p1, q1, p2, q2, bestM2[q1+q2][q2][make_pair(p1, p2)]));
        // stk.push(make_tuple(p1, q1, p2, q2, bestM2[get_keys(p1, q1)][mapping(p1, p2)][mapping(q1, q2)]));
        stk.push(make_tuple(p1, q1, p2, q2, bestM2[q1+q2][get_keys(q1, p1, p2)]));
        tuple<string, string, string, string> m2result = get_parentheses_M2(seq1, seq2, stk, hmmalign);

        pair<string, string> pre_aln, post_aln;
        string seq1_struc, seq2_struc, seq1_aln, seq2_aln;
        switch(state.manner) {
            case MANNER_MULTI_eq_MULTI_plus_U_ALN: case MANNER_MULTI_ALN:
                // pre_aln = get_hmm_aln(i1, p1-1, i2, p2-1, MANNER_ALN, state.startHMMstate);
                // post_aln = get_hmm_aln(q1+1, j1, q2+1, j2, state.endHMMstate, MANNER_ALN);

                pre_aln = get_hmm_aln_left(i1, p1, i2, p2, MANNER_ALN, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1, q2, j2, MANNER_ALN, MANNER_ALN);

                seq1_struc = struc_p2p_pair(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_pair(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            case MANNER_MULTI_eq_MULTI_plus_U_INS1: case MANNER_MULTI_INS1:
                // pre_aln = get_hmm_aln(i1, p1-1, i2+1, p2-1, MANNER_INS1, state.startHMMstate);
                // post_aln = get_hmm_aln(q1+1, j1, q2+1, j2-1, state.endHMMstate, MANNER_INS1);

                pre_aln = get_hmm_aln_left(i1, p1, i2+1, p2, MANNER_INS1, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1, q2, j2-1, MANNER_ALN, MANNER_INS1);

                seq1_struc = struc_p2p_pair(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_dots(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            case MANNER_MULTI_eq_MULTI_plus_U_INS2: case MANNER_MULTI_INS2:
                // pre_aln = get_hmm_aln(i1+1, p1-1, i2, p2-1, MANNER_INS2, state.startHMMstate);
                // post_aln = get_hmm_aln(q1+1, j1-1, q2+1, j2, state.endHMMstate, MANNER_INS2);

                pre_aln = get_hmm_aln_left(i1+1, p1, i2, p2, MANNER_INS2, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1-1, q2, j2, MANNER_ALN, MANNER_INS2);

                seq1_struc = struc_p2p_dots(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_pair(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            default:
                printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                assert(false);
        }

        seq1_aln = get<0>(pre_aln) + get<2>(m2result) + get<0>(post_aln);
        seq2_aln = get<1>(pre_aln) + get<3>(m2result) + get<1>(post_aln);

        // cout << seq1_struc << endl;
        // cout << seq2_struc << endl;

        return make_tuple(seq1_struc, seq2_struc, seq1_aln, seq2_aln);
    }

    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_M2(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    // M2 -> M+P
    if (!stk.empty()){
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) printf("get_parentheses_M2 manner at %d, %d, %d, %d: manner %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.manner, state.alignscore, state.seq1foldscore, state.seq2foldscore);

        int k1 = state.trace1.split;
        int k2 = state.trace2.split; // N.B.

        stk.push(make_tuple(i1, k1, i2, k2, bestM[k1+k2][get_keys(k1, i1, i2)])); // push M1
        tuple<string, string, string, string> pre_result = get_parentheses_M1(seq1, seq2, stk, hmmalign);

        stack<tuple<int, int, int, int, State>> stk2; // push P
        // stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestPALN[j1+j2][get_keys(j1, k1+1, k2+1)])); 
        stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][get_keys(j1, k1+1, k2+1)].alnobj)); 
        tuple<string, string, string, string> result = get_parentheses_P(seq1, seq2, stk2, hmmalign);

        return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
    }
    
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_M1(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    // M->P, M->M2, M->M+U
    if (!stk.empty()){
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) printf("get_parentheses_M1 manner at %d, %d, %d, %d: manner %d %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.premanner, state.manner, state.alignscore,  state.seq1foldscore, state.seq2foldscore);
        
        tuple<string, string, string, string> result;
        switch(state.manner) {
            case MANNER_M_eq_M_plus_U_ALN:
                if (verbose) cout << "MANNER_M_eq_M_plus_U_ALN" << endl;
                stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[j1-1+j2-1][get_keys(j1-1, i1, i2)]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result)+".", get<1>(result)+".", get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+seq2.raw_seq.at(j2));
                break;
            case MANNER_M_eq_M_plus_U_INS1:
                if (verbose) cout << "MANNER_M_eq_M_plus_U_INS1" << endl;
                stk.push(make_tuple(i1, j1-1, i2, j2, bestM[j1-1+j2][get_keys(j1-1, i1, i2)])); // [j1+j2-1][j2][make_pair(i1, i2)]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result)+".", get<1>(result), get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+"-");
                break;
            case MANNER_M_eq_M_plus_U_INS2:
                if (verbose) cout << "MANNER_M_eq_M_plus_U_INS2" << endl;
                stk.push(make_tuple(i1, j1, i2, j2-1, bestM[j1+j2-1][get_keys(j1, i1, i2)])); // [j1+j2-1][j2-1][make_pair(i1, i2)]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result), get<1>(result)+".", get<2>(result)+"-", get<3>(result)+seq2.raw_seq.at(j2));
                break;
            case MANNER_M_eq_M2:
                if (verbose) cout << "MANNER_M_eq_M2" << endl;
                stk.push(make_tuple(i1, j1, i2, j2, bestM2[j1+j2][get_keys(j1, i1, i2)])); // [j1+j2][j2][make_pair(i1, i2)]));
                return get_parentheses_M2(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_M_eq_P:
                if (verbose) cout << "MANNER_M_eq_P" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestPALN[j1+j2][get_keys(j1, i1, i2)]));
                stk.push(make_tuple(i1, j1, i2, j2, bestP[j1+j2][get_keys(j1, i1, i2)].alnobj));
                return get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;  
            default:
                printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                assert(false);
        }
    }
    
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_P(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    // P->P, H->P, Multi->P
    if ( !stk.empty() ) {
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) printf("get_parentheses_P manner at %d, %d, %d, %d: manner %d %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.premanner, state.manner, state.alignscore, state.seq1foldscore, state.seq2foldscore);

        // H to P or Multi to P
        // the most inside base pair in stacking
        switch(state.manner) {
            case MANNER_HAIRPIN_ALN:
                if (verbose) cout << "MANNER_HAIRPIN_ALN" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestHALN[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][get_keys(j1, i1, i2)].alnobj)); 
                
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_HAIRPIN_INS1:
                if (verbose) cout << "MANNER_HAIRPIN_INS1" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestHINS1[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][get_keys(j1, i1, i2)].ins1obj)); 
                
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_HAIRPIN_INS2:
                if (verbose) cout << "MANNER_HAIRPIN_INS2" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestHINS2[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][get_keys(j1, i1, i2)].ins2obj)); 
                
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_ALN:
                if (verbose) cout << "MANNER_P_eq_MULTI_ALN" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestMultiALN[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][get_keys(j1, i1, i2)].alnobj)); 
                
                return get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_INS1:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS1" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestMultiINS1[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][get_keys(j1, i1, i2)].ins1obj)); 
                
                return get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_INS2:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS2" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestMultiINS2[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][get_keys(j1, i1, i2)].ins2obj)); 
                
                return get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            
            // default:
            //     printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
            //     assert(false);
        }
        
        // P2P
        int seq1_l1 = state.trace1.paddings.l1;
        int seq1_l2 = state.trace1.paddings.l2;

        int seq2_l1 = state.trace2.paddings.l1;
        int seq2_l2 = state.trace2.paddings.l2;

        int p1 = i1 + seq1_l1; // i <= p < q <= j
        int q1 = j1 - seq1_l2;

        int p2 = i2 + seq2_l1;
        int q2 = j2 - seq2_l2;

        // if (verbose) printf("get_parentheses_P manner at %d, %d, %d, %d: manner %d %d\n", p1, q1, p2, q2, state.premanner, state.manner);

        // get inner result 
        tuple<string, string, string, string> result;
        switch(state.premanner)
        {
            case MANNER_SINGLE_ALN: 
                if (verbose) cout << "MANNER_SINGLE_ALN" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestPALN[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].alnobj)); 
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_SINGLE_INS1:
                if (verbose) cout << "MANNER_SINGLE_INS1" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestPINS1[q1+q2][get_keys(q1, p1, p2)]));
                stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins1obj));
                result = get_parentheses_P(seq1, seq2, stk, hmmalign); 
                break;

            case MANNER_SINGLE_INS2:
                if (verbose) cout << "MANNER_SINGLE_INS2" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestPINS2[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins2obj)); 
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_ALN: 
                if (verbose) cout << "MANNER_HAIRPIN_ALN" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestHALN[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][get_keys(q1, p1, p2)].alnobj)); 
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_INS1:
                if (verbose) cout << "MANNER_HAIRPIN_INS1" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestHINS1[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][get_keys(q1, p1, p2)].ins1obj)); 
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_INS2:
                if (verbose) cout << "MANNER_HAIRPIN_INS2" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestHINS2[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][get_keys(q1, p1, p2)].ins2obj)); 
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_ALN: 
                if (verbose) cout << "MANNER_P_eq_MULTI_ALN" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestMultiALN[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].alnobj));  
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_INS1:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS1" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestMultiINS1[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].ins1obj)); 
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_INS2:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS2" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestMultiINS2[q1+q2][get_keys(q1, p1, p2)])); 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].ins2obj)); 
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;

            default:
                printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                assert(false);
        }

        // get pre/post alignment and structure
        pair<string, string> pre_aln, post_aln;
        string seq1_struc, seq2_struc;
        {
            int aln_p1 = p1;
            int aln_p2 = p2;
            int aln_q1 = q1;
            int aln_q2 = q2;

            HMMManner hmm_pre_end, hmm_post_start; // = state.startHMMstate;
            // HMMManner hmm_post_start = state.endHMMstate;

            switch (state.premanner)
            {
                case MANNER_SINGLE_ALN: case MANNER_HAIRPIN_ALN: case MANNER_P_eq_MULTI_ALN: 
                    {
                        hmm_pre_end = MANNER_ALN;
                        hmm_post_start = MANNER_ALN;
                    }
                    break;
                case MANNER_SINGLE_INS1: case MANNER_HAIRPIN_INS1: case MANNER_P_eq_MULTI_INS1:
                    {
                        // aln_p2++;
                        // aln_q2--; 

                        hmm_pre_end = MANNER_INS1;
                        hmm_post_start = MANNER_INS1;
                        break;
                    }
                case MANNER_SINGLE_INS2: case MANNER_HAIRPIN_INS2: case MANNER_P_eq_MULTI_INS2:
                    {
                        // aln_p1++;
                        // aln_q1--; // N.B.
                        
                        hmm_pre_end = MANNER_INS2;
                        hmm_post_start = MANNER_INS2;

                        break;
                    }
                default:
                    printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                    assert(false);
            }

            switch (state.manner) {
                case MANNER_SINGLE_ALN: 
                    {
                        // pre_aln = get_hmm_aln(i1, aln_p1-1, i2, aln_p2-1, MANNER_ALN, hmm_pre_end);
                        // post_aln = get_hmm_aln(aln_q1+1, j1, aln_q2+1, j2, hmm_post_start, MANNER_ALN);
                        
                        pre_aln = get_hmm_aln_left(i1, aln_p1, i2, aln_p2, MANNER_ALN, hmm_pre_end);
                        post_aln = get_hmm_aln_right(aln_q1, j1, aln_q2, j2, hmm_post_start, MANNER_ALN);
                        
                        seq1_struc = struc_p2p_pair(get<0>(result), seq1_l1-1, seq1_l2-1);
                        seq2_struc = struc_p2p_pair(get<1>(result), seq2_l1-1, seq2_l2-1);
                        break;
                    }
                case MANNER_SINGLE_INS1:
                    {
                        // pre_aln = get_hmm_aln(i1, aln_p1-1, i2+1, aln_p2-1, MANNER_INS1, hmm_pre_end);
                        // post_aln = get_hmm_aln(aln_q1+1, j1, aln_q2+1, j2-1, hmm_post_start, MANNER_INS1);

                        pre_aln = get_hmm_aln_left(i1, aln_p1, i2+1, aln_p2, MANNER_INS1, hmm_pre_end);
                        post_aln = get_hmm_aln_right(aln_q1, j1, aln_q2, j2-1, hmm_post_start, MANNER_INS1);

                        seq1_struc = struc_p2p_pair(get<0>(result), seq1_l1-1, seq1_l2-1);
                        seq2_struc = struc_p2p_dots(get<1>(result), seq2_l1-1, seq2_l2-1);
                        break;
                    }
                case MANNER_SINGLE_INS2:
                    {
                        // pre_aln = get_hmm_aln(i1+1, aln_p1-1, i2, aln_p2-1, MANNER_INS2, hmm_pre_end);
                        // post_aln = get_hmm_aln(aln_q1+1, j1-1, aln_q2+1, j2, hmm_post_start, MANNER_INS2);

                        pre_aln = get_hmm_aln_left(i1+1, aln_p1, i2, aln_p2, MANNER_INS2, hmm_pre_end);
                        post_aln = get_hmm_aln_right(aln_q1, j1-1, aln_q2, j2, hmm_post_start, MANNER_INS2);

                        seq1_struc = struc_p2p_dots(get<0>(result), seq1_l1-1, seq1_l2-1);
                        seq2_struc = struc_p2p_pair(get<1>(result), seq2_l1-1, seq2_l2-1);
                        break;
                    }
                
                default:  // MANNER_NONE or other cases
                    printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                    assert(false);
            }
        }

        // combine inside and outside alignment/structure
        string seq1_aln = get<0>(pre_aln) + get<2>(result) + get<0>(post_aln);
        string seq2_aln = get<1>(pre_aln) + get<3>(result) + get<1>(post_aln);

        // cout << seq1_struc <<  " " << seq2_struc << endl;
        // cout << seq1_aln << " " << seq2_aln << endl;

        return make_tuple(seq1_struc, seq2_struc, seq1_aln, seq2_aln);
    }
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_C(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    if ( !stk.empty() ) {
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) cout << "get_parentheses_C: " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.manner << " " << state.premanner << " " << state.alignscore << " " << state.seq1foldscore << " " << state.seq2foldscore<< endl;

        if (j1 == 0 && j2 == 0) return make_tuple("", "", "", "");

        switch(state.manner) {
            case MANNER_C_eq_C_plus_P:
                {   
                    int k1 = state.trace1.split;
                    int k2 = state.trace2.split;

                    if (verbose) cout <<  "MANNER_C_eq_C_plus_P " << k1 << " " << k2 << endl;

                    HMMManner manner = state.startHMMstate;
                    switch (manner)
                    {
                        case MANNER_ALN:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].alnobj));
                            break;
                        case MANNER_INS1:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins1obj));
                            break;
                        case MANNER_INS2:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1])); // push C
                    tuple<string, string, string, string> pre_result = get_parentheses_C(seq1, seq2, stk, hmmalign);

                    stack<tuple<int, int, int, int, State>> stk2; // push P
                    // stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestPALN[j1+j2][get_keys(j1, k1+1, k2+1)])); 
                    stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][get_keys(j1, k1+1, k2+1)].alnobj)); 
                    tuple<string, string, string, string> result = get_parentheses_P(seq1, seq2, stk2, hmmalign);

                    return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                    break;
                }
                
            case MANNER_C_eq_P:
                {
                    if (verbose) cout <<  "MANNER_C_eq_P" << endl;
                    // stk.push(make_tuple(i1+1, j1, i2+1, j2, bestPALN[j1+j2][get_keys(j1, i1+1, i2+1)])); // N.B.
                    stk.push(make_tuple(i1+1, j1, i2+1, j2, bestP[j1+j2][get_keys(j1, i1+1, i2+1)].alnobj)); // N.B.
                    
                    return get_parentheses_P(seq1, seq2, stk, hmmalign);
                    break;
                }
            case MANNER_C_eq_C_plus_U_ALN:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_ALN" << endl;
                    HMMManner manner = state.startHMMstate;
                    switch (manner)
                    {
                        case MANNER_ALN:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1+j2-1][j1-1].alnobj));
                            break;
                        case MANNER_INS1:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1+j2-1][j1-1].ins1obj));
                            break;
                        case MANNER_INS2:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1+j2-1][j1-1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1+j2-1][j1-1]));
                    tuple<string, string, string, string> result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    return make_tuple(get<0>(result)+".", get<1>(result)+".", get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+seq2.raw_seq.at(j2));
                    break;
                }
            case MANNER_C_eq_C_plus_U_INS1:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_INS1" << endl;
                    HMMManner manner = state.startHMMstate;
                    switch (manner)
                    {
                        case MANNER_ALN:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1+j2][j1-1].alnobj));
                            break;
                        case MANNER_INS1:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1+j2][j1-1].ins1obj));
                            break;
                        case MANNER_INS2:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1+j2][j1-1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1+j2][j1-1]));
                    tuple<string, string, string, string> result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    return make_tuple(get<0>(result)+".", get<1>(result), get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+"-");
                    break;
                }
            case MANNER_C_eq_C_plus_U_INS2:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_INS2" << endl;
                    HMMManner manner = state.startHMMstate;
                    switch (manner)
                    {
                        case MANNER_ALN:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1+j2-1][j1].alnobj));
                            break;
                        case MANNER_INS1:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1+j2-1][j1].ins1obj));
                            break;
                        case MANNER_INS2:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1+j2-1][j1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1+j2-1][j1]));
                    tuple<string, string, string, string> result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    return make_tuple(get<0>(result), get<1>(result)+".", get<2>(result)+"-", get<3>(result)+seq2.raw_seq.at(j2)); 
                    break;
                }
                
            default:
                printf("wrong manner at %d, %d, %d, %d: manner %d\n", i1, j1, i2, j2, state.manner); fflush(stdout);
                assert(false);
        }

    }
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> BeamSankoffParser::get_parentheses(SeqObject& seq1, SeqObject& seq2, BeamAlign &hmmalign) {
    stack<tuple<int, int, int, int, State>> stk;
    HMMManner manner = bestC[seq1.seq_len + seq2.seq_len][seq1.seq_len].alnobj.startHMMstate;
    switch (manner)
    {
        case MANNER_ALN:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1.seq_len - 1 + seq2.seq_len-1][seq1.seq_len-1].alnobj));
            break;
        case MANNER_INS1:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1.seq_len - 1 + seq2.seq_len-1][seq1.seq_len-1].ins1obj));
            break;
        case MANNER_INS2:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1.seq_len - 1 + seq2.seq_len-1][seq1.seq_len-1].ins2obj));
            break;
        
        default:
            break;
    }
    // stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1.seq_len - 1 + seq2.seq_len-1][seq1.seq_len-1]));

    return get_parentheses_C(seq1, seq2, stk, hmmalign);
}

float BeamSankoffParser::get_hmm_score(int i1, int j1, int i2, int j2, int s1, bool allowout){
    // (...)
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return xlog(0.0);
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return xlog(0.0);

    // cout << "get_hmm_score " << s1 << " " << i1 << " " << j1 << " " << i2 << " " << j2 << endl; // " " << hmmalign.all_local_scores[s1][i1][i2][j1][j2]<< endl;

    // return hmmalign.all_local_scores[s1][i1][i2-hmmalign.low_bounds[i1]][j1-i1+1][j2-i2+1];
    return hmmalign.all_local_scores[s1][i1][i2][j1][j2];
}

float BeamSankoffParser::get_hmm_score_left(int i1, int j1, int i2, int j2, int s1, int s2){
    // (-->(
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return xlog(0.0);
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return xlog(0.0);

    // cout << "get_hmm_score_left " << s1 << " " << s2 << " " << i1 << " " << j1 << " " << i2 << " " << j2  << endl; //" " << hmmalign.left_local_scores[s1][s2][i1][i2][j1][j2]<< endl;

    // return hmmalign.left_local_scores[s1][s2][i1][i2-hmmalign.low_bounds[i1]][j1-i1+1][j2-i2+1];
    return hmmalign.left_local_scores[s1][s2][i1][i2][j1][j2];
}

float BeamSankoffParser::get_hmm_score_right(int i1, int j1, int i2, int j2, int s1, int s2, bool allowout){
    // )-->)
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return xlog(0.0);
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return xlog(0.0);

    if ((j1-i1+1 > 30) || (j2-i2+1 > 30))
        return hmmalign.viterbi_path_local_right(i1, j1, i2, j2, static_cast<HMMManner>(s1+1), static_cast<HMMManner>(s2+1));

    // cout << "get_hmm_score_right " << s1 << " " << s2 << " " << i1 << " " << j1 << " " << i2 << " " << j2 << endl; // " " << hmmalign.right_local_scores[s1][s2][i1][i2][j1][j2] << endl;

    // return hmmalign.right_local_scores[s1][s2][i1][i2-hmmalign.low_bounds[i1]][j1-i1+1][j2-i2+1];
    return hmmalign.right_local_scores[s1][s2][i1][i2][j1][j2];
}

pair<string, string> BeamSankoffParser::get_hmm_aln(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2){
    // cout << "hairpin: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;

    vector<char> aln1, aln2;

    hmmalign.viterbi_path_local(i1, j1, i2, j2, s1, s2);
    hmmalign.traceback2(i1, j1, i2, j2, aln1, aln2, s2);

    reverse(aln1.begin(), aln1.end());
    reverse(aln2.begin(), aln2.end());
    // for(auto e : aln1) cout << e;
    // cout << endl;
    // for(auto e : aln2) cout << e;
    // cout << endl;

    string seq1aln(aln1.begin(), aln1.end());
    string seq2aln(aln2.begin(), aln2.end());

    // cout << seq1aln << " " << seq2aln << endl;

    aln1.clear();
    aln2.clear();

    return make_pair(seq1aln, seq2aln);
}

pair<string, string> BeamSankoffParser::get_hmm_aln_left(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2){
    // cout << "left: " << i1 << " " << j1 << " " << i2 << " " << j2 <<  " " << s1 << " " << s2 << endl;

    vector<char> aln1, aln2;

    hmmalign.viterbi_path_local_left(i1, j1, i2, j2, s1, s2);
    hmmalign.traceback2(i1, j1, i2, j2, aln1, aln2, s2);

    reverse(aln1.begin(), aln1.end());
    reverse(aln2.begin(), aln2.end());
    // for(auto e : aln1) cout << e;
    // cout << endl;
    // for(auto e : aln2) cout << e;
    // cout << endl;

    string seq1aln, seq2aln;
    int count = 0;
    for(char e : aln1) {
        if (count == aln1.size() - 1)
            break;
        seq1aln.push_back(e);
        count++;
    }
    count = 0;
    for(char e : aln2) {
       if (count == aln2.size() - 1)
            break;
        seq2aln.push_back(e);
        count++;
    }
    
    // cout << seq1aln << " " << seq2aln << endl;

    aln1.clear();
    aln2.clear();

    return make_pair(seq1aln, seq2aln);
}

pair<string, string> BeamSankoffParser::get_hmm_aln_right(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2){
    // cout << "right: " << i1 << " " << j1 << " " << i2 << " " << j2 <<  " " << s1 << " " << s2 << endl;

    vector<char> aln1, aln2;

    hmmalign.viterbi_path_local_right(i1, j1, i2, j2, s1, s2, false);
    hmmalign.traceback2(i1, j1, i2, j2, aln1, aln2, s2);

    reverse(aln1.begin(), aln1.end());
    reverse(aln2.begin(), aln2.end());
    // for(auto e : aln1) cout << e;
    // cout << endl;
    // for(auto e : aln2) cout << e;
    // cout << endl;

    string seq1aln, seq2aln;
    int count = 0;
    for(char e : aln1) {
        if (count > 0)
            seq1aln.push_back(e);
        count++;
    }
    count = 0;
    for(char e : aln2) {
        if (count > 0)
            seq2aln.push_back(e);
        count++;
    }

    // cout << seq1aln << " " << seq2aln << endl;

    aln1.clear();
    aln2.clear();

    return make_pair(seq1aln, seq2aln);
}

int BeamSankoffParser::get_keys(int j1, int i1, int i2){
    return (j1 * seq1_len + i1) * seq2_len + i2;
}

// beam prune
unsigned long quickselect_partition(vector<tuple<float, int, int> >& scores, unsigned long lower, unsigned long upper) {
    float pivot = get<0>(scores[upper]);
    // cout << "pivot: " << pivot << endl;
    while (lower < upper) {
        // cout << "lower/uppeer: " << lower << " " << upper << endl;
        while (get<0>(scores[lower]) < pivot) ++lower;
        while (get<0>(scores[upper]) > pivot) --upper;
        if (get<0>(scores[lower]) == get<0>(scores[upper])) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}
// in-place quick-select
float quickselect(vector<tuple<float, int, int> >& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return get<0>(scores[lower]);
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return get<0>(scores[split]);
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

// beam prune
// unsigned long quickselect_partition(vector<pair<float, int> >& scores, unsigned long lower, unsigned long upper) {
//     float pivot = scores[upper].first;
//     while (lower < upper) {
//         while (scores[lower].first < pivot) ++lower;
//         while (scores[upper].first > pivot) --upper;
//         if (scores[lower].first == scores[upper].first) ++lower;
//         else if (lower < upper) swap(scores[lower], scores[upper]);

//     }
//     return upper;
// }
// // in-place quick-select
// float quickselect(vector<pair<float, int> >& scores, unsigned long lower, unsigned long upper, unsigned long k) {
//     if ( lower == upper ) return scores[lower].first;
//     unsigned long split = quickselect_partition(scores, lower, upper);
//     unsigned long length = split - lower + 1;
//     if (length == k) return scores[split].first;
//     else if (k  < length) return quickselect(scores, lower, split-1, k);
//     else return quickselect(scores, split+1, upper, k - length);
// }

float BeamSankoffParser::beam_prune(unordered_map<int, State3> &beamstep, int s, vector<unordered_map<int, int> > seq1_outside, vector<unordered_map<int, int> > seq2_outside){
    scores.clear();
    for (auto &item : beamstep) {
        State &cand = item.second.alnobj;
        int j1 = cand.j1;
        int i1 = cand.i1;
        int i2 = cand.i2;
        int j2 = s - j1;

        int k1 = i1 - 1;
        int k2 = i2 - 1;
        
        // assert (j1 > 0); // lisiz TODO: DEBUG
        // assert (i1 > 0); // lisiz TODO: DEBUG
        // assert (i2 > 0); // lisiz TODO: DEBUG
        // assert (k1 * k2 >= 0);

        // for tmRNA improve accuracy slightly 
        // if (cand.manner == MANNER_HAIRPIN_ALN) continue;

        State& prefix_C = bestC[k1+k2][k1].alnobj;
        HMMManner forward_endmanner = prefix_C.endHMMstate;
        if (forward_endmanner == HMMMANNER_NONE) {
            scores.push_back(make_pair(LOG_OF_ZERO, item.first));
            continue;
        }

        float newscore; // final score
        if (use_astar) {
            // + alignment & folding heuristic 
            int foldingscore;
            float alignscore, forward_score, backward_score;
    
#ifdef dynalign
            // alignment score
            int alignscore = prefix_C.alignscore + cand.alignscore + abs(seq2_len-j2-seq1_len+j1);
#else       
            // alignment forward score
            forward_score = prefix_C.alignscore; 
            alignscore = xlog_mul(xlog_mul(forward_score, cand.alignscore), hmmalign.trans_probs[forward_endmanner-1][2]);
            // alignment backward/heuristic score
            backward_score = aln_backward_score[j1][j2];      
            alignscore = xlog_mul(alignscore, backward_score);
#endif      
            // folding score
            foldingscore = cand.seq1foldscore + cand.seq2foldscore;
            // folding heuristic (DEBUG)
            // assert (seq1_outside[j1].find(i1) != seq1_outside[j1].end());
            // assert (seq2_outside[j2].find(i2) != seq2_outside[j2].end());
            foldingscore += seq1_outside[j1][i1] + seq2_outside[j2][i2];

            // final score
            if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN)
                newscore = LOG_OF_ZERO;
            else
                newscore = foldingscore + weight * alignscore;
        } 
        else {
            // original
            int foldingscore = (prefix_C.seq1foldscore + prefix_C.seq2foldscore) + (cand.seq1foldscore + cand.seq2foldscore);
            float alignscore = xlog_mul(xlog_mul(prefix_C.alignscore, cand.alignscore), hmmalign.trans_probs[forward_endmanner-1][2]);
            
            // deviation
            // newscore = foldingscore + weight * xlog_div(alignscore, viterbi_score); 

            if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN)
                newscore = LOG_OF_ZERO;
            else
                newscore = foldingscore + weight * alignscore;
        }
        
        scores.push_back(make_pair(newscore, item.first));
    }

    if (scores.size() <= beam) return VALUE_FMIN;
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}
float BeamSankoffParser::beam_prune(unordered_map<int, State> &beamstep, int s, vector<unordered_map<int, int> > seq1_outside, vector<unordered_map<int, int> > seq2_outside){
    scores.clear();
    // invalid_pos.clear();
    for (auto &item : beamstep) {
        State &cand = item.second;
        int j1 = cand.j1;
        int i1 = cand.i1;
        int i2 = cand.i2;
        int j2 = s - j1;
        int k1 = i1 - 1;
        int k2 = i2 - 1;

        // assert (j1 > 0); // lisiz TODO: DEBUG
        // assert (i1 > 0); // lisiz TODO: DEBUG
        // assert (i2 > 0); // lisiz TODO: DEBUG
        // assert (k1 * k2 >= 0);

        State& prefix_C = bestC[k1+k2][k1].alnobj;
        HMMManner forward_endmanner = prefix_C.endHMMstate;
        if (forward_endmanner == HMMMANNER_NONE) {
            scores.push_back(make_pair(LOG_OF_ZERO, item.first));
            continue;
        }

        float newscore; // final score
        if (use_astar) {
                ////////////////// + alignment & folding heuristic 
            float alignscore, forward_score, backward_score;
            int foldingscore, seq1_out, seq2_out;
        
            // alignment score
    #ifdef dynalign
            int alignscore = prefix_C.alignscore + cand.alignscore + abs(seq2_len-j2-seq1_len+j1);
    #else       
            // add forward score
            forward_score = prefix_C.alignscore; 
            alignscore = xlog_mul(xlog_mul(forward_score, cand.alignscore), hmmalign.trans_probs[forward_endmanner-1][2]);
            // backward/heuristic score
            switch (cand.endHMMstate)
            {
                case MANNER_ALN:
                    backward_score = aln_backward_score[j1][j2];
                    break;
                case MANNER_INS1:
                    backward_score = ins1_backward_score[j1][j2]; 
                    break;
                case MANNER_INS2:
                    backward_score = ins2_backward_score[j1][j2]; 
                    break;          
                default:
                    backward_score = aln_backward_score[j1][j2];
                    break;
            }
            alignscore = xlog_mul(alignscore, backward_score);
#endif
            // folding score
            foldingscore = cand.seq1foldscore + cand.seq2foldscore;
            // folding heuristic DEBUG
            // assert (seq1_outside[j1].find(i1) != seq1_outside[j1].end());
            // assert (seq2_outside[j2].find(i2) != seq2_outside[j2].end());
            foldingscore += seq1_outside[j1][i1] + seq2_outside[j2][i2];

            // final score
            if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN) // DEBUG
                newscore = LOG_OF_ZERO;
            else
                newscore = foldingscore + weight * alignscore;
        } else {
            // original
            int foldingscore = (prefix_C.seq1foldscore + prefix_C.seq2foldscore) + (cand.seq1foldscore + cand.seq2foldscore);
            float alignscore = xlog_mul(xlog_mul(prefix_C.alignscore, cand.alignscore), hmmalign.trans_probs[forward_endmanner-1][2]);
            
            // deviation
            // newscore = foldingscore + weight * xlog_div(alignscore, viterbi_score);

            if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN) { // DEBUG
                invalid_pos.push_back(item.first);
                continue;
            }

            newscore = foldingscore + weight * alignscore;
        }

        scores.push_back(make_pair(newscore, item.first));
    }

    if (scores.size() <= beam) return VALUE_FMIN;
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

/*
void BeamSankoffParser::outside(bool limited, const set<pair<int, int>> &allowed_pairs){
    SeqObject* seq1 = &sequences[0];
    SeqObject* seq2 = &sequences[1];

    // prepare
    // bestC_beta[sum_len - 2][seq1_len - 1].set(seq1_len-1, 0, 0, 0, 0, MANNER_NONE, MANNER_NONE, xlog(1.0), weight * xlog(1.0), MANNER_ALN, MANNER_ALN);

    // from right to left
    // verbose = true; // debug
    for (int s = sum_len - 2; s >= 0; --s) {
        if (s % 100 == 0) cout << "s: " << s << endl;
        // important otherwise cannot iterate    
        unordered_map<int, State3>& beamH = bestH[s];
        unordered_map<int, State3>& beamMulti = bestMulti[s];
        unordered_map<int, State3>& beamP = bestP[s];
        unordered_map<int, State>& beamM2 = bestM2[s];
        unordered_map<int, State>& beamM = bestM[s];
        unordered_map<int, State>& beamC = bestC[s];

        unordered_map<int, State3>& beamH_beta = bestH_beta[s];
        unordered_map<int, State3>& beamMulti_beta = bestMulti_beta[s];
        unordered_map<int, State3>& beamP_beta = bestP_beta[s];
        unordered_map<int, State>& beamM2_beta = bestM2_beta[s];
        unordered_map<int, State>& beamM_beta = bestM_beta[s];
        unordered_map<int, State>& beamC_beta = bestC_beta[s];

        // beam of C
        // C + U
        if (verbose) cout << "beam of C" << endl;
        {   
            for (auto &item : beamC) {
                int j1 = item.first;
                int j2 = s - j1;
                if (j2 < 0) continue; // debug
                State& state = item.second;
                State& state_beta = beamC_beta[j1];
                HMMManner current_manner = state.endHMMstate;

                int newscore1 = external_unpaired_score(j1 + 1, seq1);
                int newscore2 = external_unpaired_score(j2 + 1, seq2);

                if (verbose) cout << "C+U: " << j1 << " "  << j2  << " "<<  j1 << " "  << j2  << endl;
                if (verbose) cout << current_manner << " " << state.alignscore << endl;

                float trans_emit_prob, alignscore;

                // (seq1_len-1,seq1_len-1) ALN/INS1/INS2 => (seq1_len,seq2_len) ALN without emission prob.
                if (j1==seq1_len-1 && j2==seq2_len-1) {
                    alignscore = hmmalign.trans_probs[current_manner-1][2];
                    update_if_better(0, j1, 0, j2, state_beta, // bestC[j1+1][j2+1],
                                     0, 0,
                                     MANNER_NONE, MANNER_C_eq_C_plus_U_ALN,
                                     alignscore, current_manner, MANNER_ALN, weight, verbose); // seqs start from 0, j start from 1
                }

                if (j1 < seq1_len - 1 && j2 < seq2_len - 1){ // ALN
                    if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]) { // hmm constraint
                        // get next state
                        State& nextstate_beta = bestC_beta[s+2][j1+1];
                        if (nextstate_beta.startHMMstate == MANNER_ALN) {
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(current_manner, nextstate_beta.startHMMstate, j1+1, j2+1, true);
                            alignscore = xlog_mul(nextstate_beta.alignscore, trans_emit_prob);
                            
                            // cout << "ALN: " << nextstate_beta.alignscore << " " << trans_emit_prob << " " << alignscore << " " << state.alignscore << endl;
                            
                            update_if_better(0, j1, 0, j2, state_beta, // bestC[j1+1][j2+1],
                                            nextstate_beta.seq1foldscore + newscore1,
                                            nextstate_beta.seq2foldscore + newscore2,
                                            nextstate_beta.manner, MANNER_C_eq_C_plus_U_ALN,
                                            alignscore, current_manner, nextstate_beta.endHMMstate, weight, verbose); // seqs start from 0, j start from 1
                        }
                    }
                } 

                if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) { // INS1
                    if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) { // hmm constraint
                        // get next state
                        State& nextstate_beta = bestC_beta[s+1][j1+1];
                        if (nextstate_beta.startHMMstate == MANNER_INS1) {
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(current_manner, nextstate_beta.startHMMstate, j1+1, j2, true);
                            alignscore = xlog_mul(nextstate_beta.alignscore, trans_emit_prob);
                            
                            // cout << "INS1: " << nextstate_beta.alignscore << " " << trans_emit_prob << " " << alignscore << " " << state.alignscore << endl;

                            update_if_better(0, j1, 0, j2, state_beta, // bestC[j1+1][j2],
                                            nextstate_beta.seq1foldscore + newscore1,
                                            nextstate_beta.seq2foldscore,
                                            nextstate_beta.manner, MANNER_C_eq_C_plus_U_INS1,
                                            alignscore, current_manner, nextstate_beta.endHMMstate, weight, verbose);   
                        }
                    }
                } 

                if (j2 < seq2->seq_len - 1 && j1 <= seq1->seq_len - 1) { // INS2
                    if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) { // hmm constraint
                        // get next state
                        State& nextstate_beta = bestC_beta[s+1][j1];
                        if (nextstate_beta.startHMMstate == MANNER_INS2) {
                            // cout << "INS2: " << nextstate_beta.alignscore << " " << trans_emit_prob << " " << alignscore << " " << state.alignscore << endl;

                            trans_emit_prob = hmmalign.get_trans_emit_prob0(current_manner, nextstate_beta.startHMMstate, j1, j2+1, true);
                            alignscore = xlog_mul(nextstate_beta.alignscore, trans_emit_prob);
                        
                            update_if_better(0, j1, 0, j2, state_beta, // bestC[j1][j2+1],
                                            nextstate_beta.seq1foldscore,
                                            nextstate_beta.seq2foldscore + newscore2,
                                            nextstate_beta.manner, MANNER_C_eq_C_plus_U_INS2,
                                            alignscore, current_manner, nextstate_beta.endHMMstate, weight, verbose);   
                        }
                    }
                }
            } // for loop j1
        } // beam of C

        if (s == 0) break;

        // beam of M
        // for every state in M[j]
        //   1. M = M + unpaired
        if (verbose) cout << "beam of M" << endl;
        {
            for (auto &item : beamM) {
                State& state = item.second;
                int j1 = state.j1;
                int i1 = state.i1;
                int i2 = state.i2;
                int j2 = s - j1;

                // get beta state
                State& state_beta = beamM_beta[item.first];
                HMMManner current_manner = state.endHMMstate;
                
                if (verbose) cout << "M = M + U " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;

                float trans_emit_prob, alignscore;
                if (j1 < seq1->seq_len - 1 && j2 < seq2->seq_len - 1){
                    if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]) { // hmm constraint
                        // get next state
                        State& nextstate_beta = bestM_beta[s+2][get_keys(j1+1, i1, i2)];
                        if (nextstate_beta.startHMMstate == MANNER_ALN) {
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(current_manner, nextstate_beta.startHMMstate, j1+1, j2+1, true);
                            alignscore = xlog_mul(nextstate_beta.alignscore, trans_emit_prob);

                            update_if_better(i1, j1, i2, j2, state_beta,
                                            nextstate_beta.seq1foldscore,
                                            nextstate_beta.seq2foldscore,
                                            nextstate_beta.manner, MANNER_M_eq_M_plus_U_ALN,
                                            alignscore, current_manner, nextstate_beta.endHMMstate, weight, verbose);
                        }
                    }
                }

                if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) {
                    // hmm constraint
                    if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) {
                        // get next state
                        State& nextstate_beta = bestM_beta[s+1][get_keys(j1+1, i1, i2)];
                        if (nextstate_beta.startHMMstate == MANNER_INS1) {
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(current_manner, nextstate_beta.startHMMstate, j1+1, j2, true);
                            alignscore = xlog_mul(nextstate_beta.alignscore, trans_emit_prob);

                            update_if_better(i1, j1, i2, j2, state_beta,
                                            nextstate_beta.seq1foldscore,
                                            nextstate_beta.seq2foldscore,
                                            nextstate_beta.manner, MANNER_M_eq_M_plus_U_INS1,
                                            alignscore, current_manner, nextstate_beta.endHMMstate, weight, verbose);   
                        }
                    }
                }

                if (j1 <= seq1->seq_len - 1 && j2 < seq2->seq_len - 1) {
                    // hmm constraint
                    if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) {
                        // get next state
                        State& nextstate_beta = bestM_beta[s+1][get_keys(j1, i1, i2)];
                        if (nextstate_beta.startHMMstate == MANNER_INS2) {
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(current_manner, nextstate_beta.startHMMstate, j1, j2+1, true);
                            alignscore = xlog_mul(nextstate_beta.alignscore, trans_emit_prob);

                            update_if_better(i1, j1, i2, j2, state_beta,
                                            nextstate_beta.seq1foldscore,
                                            nextstate_beta.seq2foldscore,
                                            nextstate_beta.manner, MANNER_M_eq_M_plus_U_INS2,
                                            alignscore, current_manner, nextstate_beta.endHMMstate, weight, verbose);  
                        } 
                    }
                }
            } // for loop j1
        }// beam of M

        // beam of M2
        // for every state in M2[j]
        //   1. multi-loop  (by extending M2 on the left)
        //   2. M = M2
        if (verbose) cout << "beam of M2" << endl;
        {
            for (auto &item : beamM2) {
                State& state = item.second;
                int j1 = state.j1;
                int i1 = state.i1;
                int i2 = state.i2;
                int j2 = s - j1;

                // get beta state
                State& state_beta = beamM2_beta[item.first];

                // 1. multi-loop
                {
                    for (int p1 = i1-1; p1 >= max(i1 - SINGLE_MAX_LEN, 1); --p1) {
                        int nucp1 = seq1->nucs[p1];
                        int q1 = seq1->next_pair[nucp1][j1];
                        
                        if (q1 != -1 && ((i1 - p1 - 1) <= SINGLE_MAX_LEN)) {
                            // the current shape is p..i M2 j ..q
                            // for (int p2 = i2-1; p2 >= max(i2 - MAX_LOOP_LEN, 1); --p2) {
                            for (int p2 = min(i2-1, hmmalign.up_bounds[p1]); p2 >= max(max(1, i2 - SINGLE_MAX_LEN), hmmalign.low_bounds[p1]); --p2) {
                                // hmm constraint
                                // if (p2 < hmmalign.low_bounds[p1] || p2 > hmmalign.up_bounds[p1]) continue;
                                
                                int nucp2 = seq2->nucs[p2];
                                int q2 = seq2->next_pair[nucp2][j2];

                                // hmm constraint, speed up
                                // hmmalign.low_bounds[q1] - 2: q2 could be hmmalign.low_bounds[q1] - 1
                                if (q2 < max(j2, hmmalign.low_bounds[q1]))
                                    q2 = seq2->next_pair[nucp2][max(j2, hmmalign.low_bounds[q1] - 1)];

                                if (q2 != -1 && ((i2 - p2 - 1) <= SINGLE_MAX_LEN)) {
                                    int newscore1 = multi_unpaired_score2(i1, j1, p1, q1, seq1);
                                    int newscore2 = multi_unpaired_score2(i2, j2, p2, q2, seq2);

                                    if (verbose) cout << "M2 to Multi: " << i1 << " " << j1 << " " << i2 << " " << j2  << " " << p1 << " " << q1 << " " << p2 << " " << q2 << endl;

                                    float pre_align_trans, post_align_trans, alignscore;
                                    State3& nextstate3_beta = bestMulti_beta[q1+q2][get_keys(q1, p1, p2)];
                                    // update bestMultiALN
                                    // align cost = (p, i) + (i, j) + (j, q);
                                    {
                                        // get next state
                                        State& nextstate_beta = nextstate3_beta.alnobj;

                                        // left/right alignment score
                                        // Note: (->(...)<-) do not include the innermost alignscore
                                        pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, 2);
                                        post_align_trans = get_hmm_score_right(j1, q1, j2, q2, 2, 2);
                                        
                                        // sum of alignment score
                                        alignscore = xlog_mul(pre_align_trans, nextstate_beta.alignscore);
                                        alignscore = xlog_mul(alignscore, post_align_trans);
                                        // cout << pre_align_trans << " " << nextstate_beta.alignscore << " " << post_align_trans << endl;

                                        update_if_better(i1, j1, i2, j2, state_beta,
                                                         nextstate_beta.seq1foldscore + newscore1,
                                                         nextstate_beta.seq2foldscore + newscore2, 
                                                         nextstate_beta.manner, MANNER_MULTI_ALN,
                                                         static_cast<char>(i1 - p1), q1 - j1,
                                                         static_cast<char>(i2 - p2), q2 - j2,
                                                         alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                    }
                                    // update bestMultiINS1
                                    {
                                        // get next state
                                        State& nextstate_beta = nextstate3_beta.ins1obj;

                                        // left/right alignment score:
                                        // (p1, q1) is inserted, do not include (p2, q2) 
                                        pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, 2);
                                        post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, 2, 0);
                                        
                                        // sum of alignment score
                                        alignscore = xlog_mul(pre_align_trans, nextstate_beta.alignscore);
                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                        update_if_better(i1, j1, i2, j2, state_beta,
                                                         nextstate_beta.seq1foldscore + newscore1,
                                                         nextstate_beta.seq2foldscore + newscore2, 
                                                         nextstate_beta.manner, MANNER_MULTI_INS1,
                                                         static_cast<char>(i1 - p1), q1 - j1,
                                                         static_cast<char>(i2 - p2), q2 - j2,
                                                         alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                    }
                                    // update bestMultiINS2
                                    {
                                        // get next state
                                        State& nextstate_beta = nextstate3_beta.ins2obj;

                                        // left/right alignment score:
                                        // (p2, q2) is inserted, do not include (p1, q1) 
                                        pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, 2);
                                        post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, 2, 1);
                                        
                                        // sum of alignment score
                                        alignscore = xlog_mul(pre_align_trans, nextstate_beta.alignscore);
                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                        update_if_better(i1, j1, i2, j2, state_beta, // bestMultiINS2[q1+q2][get_keys(q1, p1, p2)],
                                                         nextstate_beta.seq1foldscore + newscore1,
                                                         nextstate_beta.seq2foldscore + newscore2, 
                                                         nextstate_beta.manner, MANNER_MULTI_INS2,
                                                         static_cast<char>(i1 - p1), q1 - j1,
                                                         static_cast<char>(i2 - p2), q2 - j2,
                                                         alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                    }                          
                                } // q2
                            } // p2
                        } // q1
                    } // p1
                }// 1. multi-loop

                // 2. M = M2
                {
                    if (verbose) cout << "M2 to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                    State& nextstate_beta = beamM_beta[item.first];
                    // TODO: check manner, M must be M + P ??
                    update_if_better(i1, j1, i2, j2, state_beta,
                                     nextstate_beta.seq1foldscore,
                                     nextstate_beta.seq2foldscore,
                                     nextstate_beta.manner, MANNER_M_eq_M2,
                                     nextstate_beta.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                }// 2. M = M2
            }// end beamM2
        }// beam of M2
     
        // beam of P
        // for every state in P[j]
        //   1. generate new helix/bulge
        //   2. M = P
        //   3. M2 = M + P
        //   4. C = C + P
        if (verbose) cout << "beam of P" << endl;
        {
            for (auto &item : beamP) {
                for (int m=0; m<3; m++){
                    State *stateptr, *state_betaptr;
                    switch (m)
                    {
                        case 0:
                            stateptr = &item.second.ins1obj;
                            state_betaptr = &beamP_beta[item.first].ins1obj;
                            break;
                        case 1:
                            stateptr = &item.second.ins2obj;
                            state_betaptr = &beamP_beta[item.first].ins2obj;
                            break;
                        case 2:
                            stateptr = &item.second.alnobj;
                            state_betaptr = &beamP_beta[item.first].alnobj;
                            break;
                        default:
                            break;
                    }
                    State& state = *stateptr;
                    State& state_beta = *state_betaptr;

                    // State state = item.second;
                    int j1 = state.j1;
                    int i1 = state.i1;
                    int i2 = state.i2;
                    int j2 = s - j1;
                    if (j1 < 0 || i1 < 0 || i2 < 0) continue; // lisiz TODO: DEBUG

                    // multilign
                    if (limited) {
                        if (allowed_pairs.find(make_pair(i1, j1)) == allowed_pairs.end()) continue;
                    }

                    // 1. generate new helix / single_branch
                    // new state is of shape p..i..j..q
                    // Note: p >= 1
                    if (i1 > 1 && j1 < seq1->seq_len - 1){
                        if (i2 > 1 && j2 < seq2->seq_len - 1){
                            int nuci1 = seq1->nucs[i1];
                            int nuci1_1 = seq1->nucs[i1 - 1];
                            int nucj1 = seq1->nucs[j1];
                            int nucj1p1 = seq1->nucs[j1 + 1];

                            int nuci2 = seq2->nucs[i2];
                            int nuci2_1 = seq2->nucs[i2 - 1];
                            int nucj2 = seq2->nucs[j2];
                            int nucj2p1 = seq2->nucs[j2 + 1];

                            // p1<i1, q1>j1; p2<i2, q2>j2
                            {
                                // for (int p1 = i1 - 1; p1 >= max(i1 - SINGLE_MAX_LEN, 1); --p1) {
                                for (int p1 = i1; p1 >= max(i1 - MAX_LOOP_LEN, 1); --p1) {
                                    int nucp1 = seq1->nucs[p1];
                                    int nucp1p1 = seq1->nucs[p1 + 1]; // hzhang: move here
                                    // int q1 = seq1->next_pair[nucp1][j1];
                                    int q1 = p1 == i1 ? j1 : seq1->next_pair[nucp1][j1];

                                    while (q1 != -1 && ((i1 - p1) + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                                        if (m == 1) { // INS2 Noted
                                            if (p1 != i1 || j1 != q1) break;
                                        } else {
                                            if (p1 == i1 || j1 == q1) break;
                                        }
                                        
                                        int nucq1 = seq1->nucs[q1];
                                        int nucq1_1 = seq1->nucs[q1 - 1];

                                        int p2p1 = 0; 
                                        if (p1 < i1 && q1 > j1) p2p1 = P2PScore(p1,q1,i1,j1,nucp1,nucp1p1,nucq1_1,nucq1,nuci1_1,nuci1,nucj1,nucj1p1); 

                                        // for (int p2 = min(i2, hmmalign.up_bounds[p1]); p2 >= max(i2 - MAX_LOOP_LEN, hmmalign.low_bounds[p1]); --p2) { // seq2 for loop
                                        for (int p2 = min(i2, hmmalign.up_bounds[p1]); p2 >= max(max(1, i2 - MAX_LOOP_LEN), hmmalign.low_bounds[p1]); --p2) { // seq2 for loop
                                            // if (p2 < hmmalign.low_bounds[p1] || p2 > hmmalign.up_bounds[p1]) continue;
                                            int nucp2 = seq2->nucs[p2];
                                            int nucp2p1 = seq2->nucs[p2 + 1]; // hzhang: move here
                                            int q2 = p2 == i2 ? j2 : seq2->next_pair[nucp2][j2];

                                            // speed up
                                            // q2 could be j2 or hmmalign.low_bounds[q1]-1
                                            if ((max(j2, hmmalign.low_bounds[q1]) - 1) > q2)
                                                q2 = seq2->next_pair[nucp2][max(j2, hmmalign.low_bounds[q1]) - 1]; // TODO

                                            while (q2 <= (hmmalign.up_bounds[q1]) && q2 != -1 && ((i2 - p2) + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                                if (m == 0) { // INS1 Noted
                                                    if (p2 != i2 || j2 != q2) break;
                                                } else {
                                                    if (p2 == i2 || j2 == q2) break;
                                                }

                                                // TODO: redundant calculation
                                                int nucq2 = seq2->nucs[q2];
                                                int nucq2_1 = seq2->nucs[q2 - 1];

                                                int p2p2 = 0; 
                                                if (p2 < i2 && q2 > j2) p2p2 = P2PScore(p2,q2,i2,j2,nucp2,nucp2p1,nucq2_1,nucq2,nuci2_1,nuci2,nucj2,nucj2p1);

                                                // if (verbose) 
                                                // if (m==2 && i1==2 && j1==115 && i2==1 && j2==115)
                                                //     cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;
                                                
                                                // get next state3
                                                int key = get_keys(q1, p1, p2);
                                                if (bestP_beta[q1+q2].find(key) != bestP_beta[q1+q2].end()) {
                                                    State3& nextstate3_beta = bestP_beta[q1+q2][key];
                                                    pair<float, HMMManner> pre_alignscore, post_alignscore;
                                                    float pre_align_trans, post_align_trans, alignscore;
                                                    // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj), si == sj == ALN/INS1/INS2
                                                    {   
                                                        // get next state
                                                        State& nextstate_beta = nextstate3_beta.alnobj;
                                                                                    
                                                        // si == sj == ALN
                                                        // align cost = (p, i) + (i, j) + (j, q)
                                                        // Note: (->(...)<-) do not include outsidemost alignscore
                                                        pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                                                        post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                                                        
                                                        alignscore = xlog_mul(pre_align_trans, nextstate_beta.alignscore);
                                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                                        update_if_better(i1, j1, i2, j2, state_beta,
                                                                        nextstate_beta.seq1foldscore + p2p1,
                                                                        nextstate_beta.seq2foldscore + p2p2,
                                                                        nextstate_beta.manner, MANNER_SINGLE_ALN, 
                                                                        static_cast<char>(i1 - p1), q1 - j1,
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                                    }
                                                        
                                                    // if (p2+1-hmmalign.low_bounds[p1] >=0)
                                                    {
                                                        // get next state
                                                        State& nextstate_beta = nextstate3_beta.ins1obj;

                                                        // si == sj == INS1
                                                        // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                                        // Note: (->(...)<-) do not include outsidemost alignscore
                                                        pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                                                        post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0);
                                                        
                                                        alignscore = xlog_mul(pre_align_trans, nextstate_beta.alignscore);
                                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                                        // debug: 5 112 4 112 2 
                                                        // if (m==2 && i1==2 && j1==115 && i2==1 && j2==115){
                                                        //     // 4 113 3 113 
                                                        //     // if (p1==4 && q1==113 && p2==3 && q2==113) {
                                                        //     cout << state_beta.score << " " << state_beta.seq1foldscore << " " << state_beta.seq2foldscore << " " << state_beta.alignscore << endl;
                                                        //     cout << "0" << " " << nextstate_beta.seq1foldscore + p2p1 + nextstate_beta.seq2foldscore + p2p2 + weight*alignscore << " ";
                                                        //     cout << " " << nextstate_beta.seq1foldscore + p2p1 << " " <<nextstate_beta.seq2foldscore + p2p2 << " " << alignscore << endl;
                                                        //     cout << pre_align_trans << " " <<  nextstate_beta.alignscore << " " << post_align_trans << endl;
                                                        //     // }
                                                        // }

                                                        update_if_better(i1, j1, i2, j2, state_beta, // bestPINS1[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                        nextstate_beta.seq1foldscore + p2p1,
                                                                        nextstate_beta.seq2foldscore + p2p2,
                                                                        nextstate_beta.manner, MANNER_SINGLE_INS1, 
                                                                        static_cast<char>(i1 - p1), q1 - j1,
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                                    }
                                                
                                                    // if (p2-hmmalign.low_bounds[p1+1] >= 0)
                                                    {
                                                        // get next state
                                                        State& nextstate_beta = nextstate3_beta.ins2obj;

                                                        // si == sj == INS2
                                                        // align cost = (p, i) + (i, j) + (j, q);
                                                        // Note: (->(...)<-) do not include outsidemost alignscore
                                                        pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, m);
                                                        post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, m, 1);
                                                        
                                                        alignscore = xlog_mul(pre_align_trans, nextstate_beta.alignscore);
                                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                                        // debug: 5 112 4 112 2 
                                                        // if (m==2 && i1==2 && j1==115 && i2==1 && j2==115){
                                                        //     // 4 113 3 113 
                                                        //     // if (p1==4 && q1==113 && p2==3 && q2==113) {
                                                        //     cout << state_beta.score << " " << state_beta.seq1foldscore << " " << state_beta.seq2foldscore << " " << state_beta.alignscore << endl;
                                                        //     cout << "1" << " " << nextstate_beta.seq1foldscore + p2p1 + nextstate_beta.seq2foldscore + p2p2 + weight*alignscore << " ";
                                                        //     cout << " " << nextstate_beta.seq1foldscore + p2p1 << " " <<nextstate_beta.seq2foldscore + p2p2 << " " << alignscore << endl;
                                                        //     cout << pre_align_trans << " " <<  nextstate_beta.alignscore << " " << post_align_trans << endl;
                                                        //     // }
                                                        // }

                                                        update_if_better(i1, j1, i2, j2, state_beta, // bestPINS2[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                        nextstate_beta.seq1foldscore + p2p1,
                                                                        nextstate_beta.seq2foldscore + p2p2,
                                                                        nextstate_beta.manner, MANNER_SINGLE_INS2, 
                                                                        static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                                    }
                                                }
                                                q2 = seq2->next_pair[nucp2][q2];
                                            } // while loop enumerate q2
                                        } // for loop enumerate p2
                                        q1 = seq1->next_pair[nucp1][q1];
                                    } // while loop enumerate q1
                                } // for loop enumerate p1
                            } // p1<i1, q1>j1; p2<i2, q2>j2
                        }  
                    } // 1. generate new helix / single_branch

                    if (m != 2) continue;

                    // 2. M = P
                    // accessible pairs in multiloop must be aligned
                    {   
                        if (verbose) cout << "P to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                        int newscore1 = branch_score(i1, j1, seq1);
                        int newscore2 = branch_score(i2, j2, seq2);

                        // get next state
                        State& nextstate_beta = beamM_beta[item.first];
                        // TODO: check manner, must be M2->M

                        update_if_better(i1, j1, i2, j2, state_beta, // bestM[key1][newi2][newj2], // beamM[j2][make_pair(i1, i2)], 
                                         nextstate_beta.seq1foldscore + newscore1,
                                         nextstate_beta.seq2foldscore + newscore2, 
                                         nextstate_beta.manner, MANNER_M_eq_P,
                                         nextstate_beta.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                    } // 2. M = P

                    // 3. M2 = M + P 
                    // accessible pairs in multiloop must be aligned
                    // TODO: use_cube_pruning
                    {
                        int k1 = i1 - 1;
                        int k2 = i2 - 1;
                        // int mnewk2 = mapping(k1, k2);

                        if (k1 > 1 && k2 > 1 && !bestM[k1 + k2].empty()) {
                            int newscore1 = branch_score(i1, j1, seq1);
                            int newscore2 = branch_score(i2, j2, seq2); 

                            float alignscore;
                            for (auto &m : bestM[k1 + k2]) {
                                State& mstate = m.second;
                                int newj1 = mstate.j1;
                                if (newj1 != k1) continue;
                                int newi1 = mstate.i1;
                                int newi2 = mstate.i2;

                                State& mstate_beta = bestM_beta[k1+k2][m.first];

                                if (verbose) cout << "M2=M+P: " << i1 << " " << j1 << " "  << i2 << " " << j2 <<  " " << newi1 << " " << j1 << " "  << newi2 << " " << j2  << endl;

                                State& m2state_beta = beamM2_beta[get_keys(j1, newi1, newi2)];
                                float alignscore;
                                if (m2state_beta.manner != 0) {
                                    // update mstate_beta
                                    alignscore = xlog_mul(m2state_beta.alignscore, state.alignscore);
                                    // alignscore = xlog_mul(alignscore, hmmalign.trans_probs[2][2]);
                                    alignscore = xlog_mul(alignscore, hmmalign.trans_probs[mstate.endHMMstate-1][2]);
                                    // cout << "M+P update M: " << alignscore << endl;
                                    update_if_better(newi1, k1, newi2, k2, mstate_beta, // bestM2[get_keys(newi1, j1)][mnewi2][newj2], // beamM2[j2][make_pair(newi1, newi2)], // Note: not i1 i2 but newi1 newi2new
                                                     m2state_beta.seq1foldscore + newscore1 + state.seq1foldscore,
                                                     m2state_beta.seq2foldscore + newscore2 + state.seq2foldscore,
                                                     mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                     alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                    // cout << bestM_beta[k1+k2][m.first].alignscore << " " << m.first << " " << get_keys(k1, newi1, newi2) << " " << bestM_beta[k1+k2][m.first].startHMMstate << endl;

                                    // update state_beta
                                    alignscore = xlog_mul(m2state_beta.alignscore, mstate.alignscore);
                                    // alignscore = xlog_mul(alignscore, hmmalign.trans_probs[2][2]);
                                    alignscore = xlog_mul(alignscore, hmmalign.trans_probs[mstate.endHMMstate-1][2]);
                                    // cout << "M+P update P: " << alignscore << endl;
                                    update_if_better(i1, j1, i2, j2, state_beta, // bestM2[get_keys(newi1, j1)][mnewi2][newj2], // beamM2[j2][make_pair(newi1, newi2)], // Note: not i1 i2 but newi1 newi2new
                                                     m2state_beta.seq1foldscore + newscore1 + mstate.seq1foldscore,
                                                     m2state_beta.seq2foldscore + newscore2 + mstate.seq2foldscore,
                                                     mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                     alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                    // cout << state_beta.seq1foldscore << " " << state_beta.seq2foldscore << " " << state_beta.alignscore << endl;
                                    // cout <<  m2state_beta.seq1foldscore + newscore1 + mstate.seq1foldscore << " " << m2state_beta.seq2foldscore + newscore2 + mstate.seq2foldscore << " " << alignscore << endl;
                                    // cout << m2state_beta.alignscore << " " <<  mstate.alignscore << endl;
                                }
                            }
                        }
                    } // 3. M2 = M + P 

                    // 4. C = C + P
                    // external pairs must be aligned
                    {   
                        int k1 = i1 - 1;
                        int k2 = i2 - 1;
                        if (k1 < 0 && k2 < 0) continue;

                        int newscore1 = external_paired_score(k1, j1, seq1);
                        int newscore2 = external_paired_score(k2, j2, seq2);

                        State& state_C = bestC[k1+k2][k1]; // bestC[k1][k2];
                        State& state_C_beta = bestC_beta[k1+k2][k1]; // bestC[k1][k2];
                        State& state_C_beta_out = beamC_beta[j1];
                        if (state_C_beta_out.startHMMstate != MANNER_ALN) continue;
                        if (state_C.endHMMstate == HMMMANNER_NONE) continue;
                        
                        if (verbose) cout << "C+P: "<< i1 << " " << j1 << " "  << i2 << " " << j2 << endl;  
                        if (verbose) cout << "state_C.endHMMstate: " << state_C.endHMMstate << endl; // debug   

                        float alignscore1, alignscore2;
                        alignscore1 = xlog_mul(state_C_beta_out.alignscore, state.alignscore); 
                        alignscore1 = xlog_mul(alignscore1, hmmalign.trans_probs[state_C.endHMMstate-1][2]);
                       
                        alignscore2 = xlog_mul(state_C_beta_out.alignscore, state_C.alignscore);
                        alignscore2 = xlog_mul(alignscore2, hmmalign.trans_probs[state_C.endHMMstate-1][2]);
                        alignscore2 = xlog_mul(alignscore2, hmmalign.trans_probs[2][state_C_beta_out.startHMMstate-1]);

                        // debug 
                        // cout << state_C.score << " " << state_C.seq1foldscore << " " << state_C.seq1foldscore << " " << state_C.alignscore << endl;
                        // cout << state.score << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << endl;

                        // update beta C
                        // if (i1 == 377 && j1 == 442 && i2 == 371 && j2 == 437) cout << "update C" << endl;
                        update_if_better(0, k1, 0, k2, state_C_beta, // bestC[j1][j2],
                                        state_C_beta_out.seq1foldscore + state.seq1foldscore + newscore1,
                                        state_C_beta_out.seq2foldscore + state.seq2foldscore + newscore2, 
                                        state_C_beta_out.manner, MANNER_C_eq_C_plus_P,
                                        k1, k2,
                                        alignscore1, state_C.endHMMstate, state_C_beta_out.endHMMstate, weight, verbose);

                        // update beta P
                        // if (i1 == 377 && j1 == 442 && i2 == 371 && j2 == 437) cout << "update P" << endl;
                        update_if_better(i1, j1, i2, j2, state_beta, // bestC[j1][j2],
                                        state_C_beta_out.seq1foldscore + state_C.seq1foldscore + newscore1,
                                        state_C_beta_out.seq2foldscore + state_C.seq2foldscore + newscore2, 
                                        state_C_beta_out.manner, MANNER_C_eq_C_plus_P,
                                        k1, k2,
                                        alignscore2, MANNER_ALN, MANNER_ALN, weight, verbose);

                        // verbose = false; // debug
                    } // 4. C = C + P
                } // end beamP
            } // loop ALN/INS1/INS2
        }// beam of P

        // beam of Multi
        // for every state in Multi[j]
        //   1. extend (i, j) to (i, jnext)
        //   2. generate P (i, j)
        if (verbose) cout << "beam of Multi" << endl;
        {
            for (auto &item : beamMulti) {
                for (int m=0; m<3; m++){
                    State *stateptr, *state_betaptr;
                    switch (m)
                    {
                        case 0:
                            stateptr = &item.second.ins1obj;
                            state_betaptr = &beamMulti_beta[item.first].ins1obj;
                            break;
                        case 1:
                            stateptr = &item.second.ins2obj;
                            state_betaptr = &beamMulti_beta[item.first].ins2obj; 
                            break;
                        case 2:
                            stateptr = &item.second.alnobj;
                            state_betaptr = &beamMulti_beta[item.first].alnobj;
                            break;
                        
                        default:
                            break;
                    }
                    State& state = *stateptr;
                    State& state_beta = *state_betaptr;

                    // State state = item.second;
                    int j1 = state.j1;
                    int i1 = state.i1;
                    int i2 = state.i2;
                    int j2 = s - j1; 
                    if (i1 < 0) continue; // TODO: debug

                    // 1. extend (i, j) to (i, jnext)
                    {
                        tuple<int, int, char, int> result1 = multiloopUnpairedScore(i1, j1, seq1, &state.trace1);
                        tuple<int, int, char, int> result2 = multiloopUnpairedScore(i2, j2, seq2, &state.trace2);

                        int j1next = get<0>(result1);
                        int j2next = get<0>(result2);

                        int new_seq1_l2 = get<3>(result1);
                        int new_seq2_l2 = get<3>(result2);

                        if (verbose) cout << "beamMulti jnext: " << m << " "  << i1  <<  " " << j1 << " " << i2 << " " << j2 << " "  << j1next  <<  " " << j2next << endl; 
                        
                        // State3& nextstate3_beta = bestMulti_beta[j1next+j2next][get_keys(j1next, i1, i2)];

                        tuple<float, HMMManner, HMMManner, string, string> alnret;
                        // pair<float, HMMManner> alignscore;
                        float alignscore;
                        int key, news;
                        if (m == 2) {
                            if (j1next != -1 && j2next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1next+j2next][get_keys(j1next, i1, i2)].alnobj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1next, j2, j2next, 2, 2); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore + get<1>(result1),
                                                nextstate_beta.seq2foldscore + get<1>(result2),
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                nextstate_beta.trace1.paddings.l1, new_seq1_l2,
                                                nextstate_beta.trace2.paddings.l1, new_seq2_l2,
                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                            }
                            if (j1next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1next+j2][get_keys(j1next, i1, i2)].alnobj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1next, j2, j2, 2, 2); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore + get<1>(result1),
                                                nextstate_beta.seq2foldscore,
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                nextstate_beta.trace1.paddings.l1, new_seq1_l2,
                                                nextstate_beta.trace2.paddings.l1, nextstate_beta.trace2.paddings.l2,
                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                            }
                            if (j2next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1+j2next][get_keys(j1, i1, i2)].alnobj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1, j2, j2next, 2, 2); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore,
                                                nextstate_beta.seq2foldscore + get<1>(result2),
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                nextstate_beta.trace1.paddings.l1, nextstate_beta.trace1.paddings.l2,
                                                nextstate_beta.trace2.paddings.l1, new_seq2_l2,
                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                            }
                        }
                                
                        else if (m == 0) {
                            if (j1next != -1 && j2next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1next+j2next][get_keys(j1next, i1, i2)].ins1obj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1next, j2, j2next-1, 0, 0); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore + get<1>(result1),
                                                nextstate_beta.seq2foldscore + get<1>(result2),
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                nextstate_beta.trace1.paddings.l1, new_seq1_l2,
                                                nextstate_beta.trace2.paddings.l1, new_seq2_l2,
                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                            }
                            if (j1next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1next+j2][get_keys(j1next, i1, i2)].ins1obj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1next, j2, j2-1, 0, 0); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore + get<1>(result1),
                                                nextstate_beta.seq2foldscore,
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                nextstate_beta.trace1.paddings.l1, new_seq1_l2,
                                                nextstate_beta.trace2.paddings.l1, nextstate_beta.trace2.paddings.l2,
                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                            }
                            if (j2next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1+j2next][get_keys(j1, i1, i2)].ins1obj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1, j2, j2next-1, 0, 0); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore,
                                                nextstate_beta.seq2foldscore + get<1>(result2),
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                nextstate_beta.trace1.paddings.l1, nextstate_beta.trace1.paddings.l2,
                                                nextstate_beta.trace2.paddings.l1, new_seq2_l2,
                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                            }
                        }
                                
                        else if (m == 1) {
                            if (j1next != -1 && j2next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1next+j2next][get_keys(j1next, i1, i2)].ins2obj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1next-1, j2, j2next, 1, 1); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore + get<1>(result1),
                                                nextstate_beta.seq2foldscore + get<1>(result2),
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                nextstate_beta.trace1.paddings.l1, new_seq1_l2,
                                                nextstate_beta.trace2.paddings.l1, new_seq2_l2,
                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                            }
                            if (j1next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1next+j2][get_keys(j1next, i1, i2)].ins2obj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1next-1, j2, j2, 1, 1); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore + get<1>(result1),
                                                nextstate_beta.seq2foldscore,
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                nextstate_beta.trace1.paddings.l1, new_seq1_l2,
                                                nextstate_beta.trace2.paddings.l1, nextstate_beta.trace2.paddings.l2,
                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                            }
                            if (j2next != -1) {
                                // hmm constraint
                                // if (j2next < hmmalign.low_bounds[j1next] || j2next > hmmalign.up_bounds[j1next]) continue;

                                // get next state
                                State& nextstate_beta = bestMulti_beta[j1+j2next][get_keys(j1, i1, i2)].ins2obj;

                                // )<-)
                                alignscore = get_hmm_score_right(j1, j1-1, j2, j2next, 1, 1); // true
                                alignscore = xlog_mul(nextstate_beta.alignscore, alignscore);
            
                                update_if_better(i1, j1, i2, j2, state_beta,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                nextstate_beta.seq1foldscore,
                                                nextstate_beta.seq2foldscore + get<1>(result2),
                                                state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                nextstate_beta.trace1.paddings.l1, nextstate_beta.trace1.paddings.l2,
                                                nextstate_beta.trace2.paddings.l1, new_seq2_l2,
                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                            }
                        }
                    } // 1. extend (i, j) to (i, jnext)

                    // 2. generate P (i, j)
                    {
                        if (verbose) cout << "Multi to P: " << m <<  " " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                        int score1 = multiloop2Pscore(i1, j1, seq1);
                        int score2 = multiloop2Pscore(i2, j2, seq2);

                        State3& nextstate3_beta = beamP_beta[item.first];

                        if (m == 2){
                            State& nexstate_beta = nextstate3_beta.alnobj;
                            update_if_better(i1, j1, i2, j2, state_beta, // beamPALN[get_keys(j1, i1, i2)],
                                             nexstate_beta.seq1foldscore + score1, 
                                             nexstate_beta.seq2foldscore + score2,
                                             nexstate_beta.manner, MANNER_P_eq_MULTI_ALN,
                                             nexstate_beta.alignscore, nexstate_beta.startHMMstate, nexstate_beta.endHMMstate, weight, verbose);
                            // cout << state_beta.seq1foldscore << " " << state_beta.seq2foldscore << " " << state_beta.alignscore << endl;
                            // cout << nexstate_beta.seq1foldscore << " " << nexstate_beta.seq2foldscore << " " << nexstate_beta.alignscore << endl;
                            // cout << score1 << " " << score2 << endl;
                            // cout << bestMulti_beta[s][item.first].alnobj.alignscore << endl;
                        }
                            
                        else if (m == 0){
                            State& nexstate_beta = nextstate3_beta.ins1obj;
                            update_if_better(i1, j1, i2, j2, state_beta, // beamPINS1[get_keys(j1, i1, i2)],
                                             nexstate_beta.seq1foldscore + score1, 
                                             nexstate_beta.seq2foldscore + score2,
                                             nexstate_beta.manner, MANNER_P_eq_MULTI_INS1,
                                             nexstate_beta.alignscore, nexstate_beta.startHMMstate, nexstate_beta.endHMMstate, weight, verbose);
                            // cout << bestMulti_beta[s][item.first].alnobj.alignscore << endl;
                        }
                            
                        else if (m == 1){
                            State& nexstate_beta = nextstate3_beta.ins2obj;
                            update_if_better(i1, j1, i2, j2, state_beta, // beamPINS2[get_keys(j1, i1, i2)], 
                                             nexstate_beta.seq1foldscore + score1, 
                                             nexstate_beta.seq2foldscore + score2,
                                             nexstate_beta.manner, MANNER_P_eq_MULTI_INS2,
                                             nexstate_beta.alignscore, nexstate_beta.startHMMstate, nexstate_beta.endHMMstate, weight, verbose);
                            // cout << bestMulti_beta[s][item.first].alnobj.alignscore << endl;
                        }
                    } // 2. generate P (i, j)
                } // beam ALN/INS1/INS2
            } // end beamMulti
        } // end beam of Multi
    } // end of s

}
*/

void BeamSankoffParser::prepare(const vector<string> &seqs){
    // get number of sequences and sequence length
    num_seqs = seqs.size();
    sequences.clear();
    // sequences.resize(num_seqs);

    sum_len = 0;
    for (int i=0; i<num_seqs; i++){
        SeqObject seq("N" + seqs[i]); // N.B. add an element at the beginning for alignment
        seq.seq_len = seqs[i].length() + 1;
        seq.nucs.clear();
        seq.nucs.resize(seq.seq_len);
        for (int j=0; j<seq.seq_len; j++){
            if (j == 0) seq.nucs[j] = 4;
            else seq.nucs[j] = GET_ACGU_NUM(seqs[i][j-1]);
        }
        // for (int j=0; j<seq.seq_len; j++){
        //     cout << seq.nucs[j] << " " ;
        // }
        // cout << endl;
        sequences.push_back(seq);
        sum_len += seq.seq_len;
    }
    // cout << "sum of seq length: " << sum_len << endl;

    // inside
    bestH.clear();
    bestH.resize(sum_len);
    bestP.clear();
    bestP.resize(sum_len);
    bestMulti.clear();
    bestMulti.resize(sum_len);

    bestM2.clear();
    bestM2.resize(sum_len);
    bestM.clear();
    bestM.resize(sum_len);

    bestC.clear();
    bestC.resize(sum_len+1);

    // outside
    bestH_beta.clear();
    bestH_beta.resize(sum_len);
    bestP_beta.clear();
    bestP_beta.resize(sum_len);
    bestMulti_beta.clear();
    bestMulti_beta.resize(sum_len);

    bestM2_beta.clear();
    bestM2_beta.resize(sum_len);
    bestM_beta.clear();
    bestM_beta.resize(sum_len);

    bestC_beta.clear();
    bestC_beta.resize(sum_len);

#ifdef is_cube_pruning
    // cube pruning
    sorted_bestM.clear();
    sorted_bestM.resize(sum_len);
#endif

    scores.clear();

    // HMM align tool
    hmmalign.set(alnbeam, sequences[0].nucs, sequences[1].nucs);
    hmmalign.viterbi_path(false); // envelope
#ifndef dynalign
    cout << "recompute forward-backward viterbi with new parameters" << endl;
    aln_viterbi = hmmalign.viterbi_path(true);
    // save alignment backward score
    aln_backward_score.resize(seq1_len);
    ins1_backward_score.resize(seq1_len);
    ins2_backward_score.resize(seq1_len);
    for (int s=0; s<sum_len-1; s++) {
        for (auto &item : hmmalign.bestALN[s]) {
            AlignState state = item.second;
            int i = state.i;
            int k = state.k;
            if (k < hmmalign.low_bounds[i] || k > hmmalign.up_bounds[i]) continue;
            aln_backward_score[i][k] = state.beta;
        } 
        for (auto &item : hmmalign.bestINS1[s]) {
            AlignState state = item.second;
            int i = state.i;
            int k = state.k;
            if (k < hmmalign.low_bounds[i] || k > hmmalign.up_bounds[i]) continue;
            ins1_backward_score[i][k] = state.beta;
        }
        for (auto &item : hmmalign.bestINS2[s]) {
            AlignState state = item.second;
            int i = state.i;
            int k = state.k;
            if (k < hmmalign.low_bounds[i] || k > hmmalign.up_bounds[i]) continue;
            ins2_backward_score[i][k] = state.beta;
        }
    }
    delete[] hmmalign.bestALN;
    delete[] hmmalign.bestINS1;
    delete[] hmmalign.bestINS2;

    hmmalign.viterbi_path_all_locals();
#endif

    // single sequence folding
    BeamCKYParser* cky_parser = new BeamCKYParser();
    cky_parser->beam = lfbeam;

    seq1_out_H.clear();
    seq1_out_P.clear();
    seq1_out_M.clear();
    seq1_out_M2.clear();
    seq1_out_Multi.clear();
    seq1_out_C.clear();
    seq2_out_H.clear();
    seq2_out_P.clear();
    seq2_out_M.clear();
    seq2_out_M2.clear();
    seq2_out_Multi.clear();
    seq2_out_C.clear();

    seq1_out_H.resize(seq1_len);
    seq1_out_P.resize(seq1_len);
    seq1_out_M.resize(seq1_len); 
    seq1_out_M2.resize(seq1_len);
    seq1_out_Multi.resize(seq1_len);
    seq1_out_C.resize(seq1_len);
    seq2_out_H.resize(seq2_len);
    seq2_out_P.resize(seq2_len);
    seq2_out_M.resize(seq2_len); 
    seq2_out_M2.resize(seq2_len);
    seq2_out_Multi.resize(seq2_len);
    seq2_out_C.resize(seq2_len);

    for (int i_seq=0; i_seq<2; i_seq++) {
        // single sequence folding
        int seq_viterbi = cky_parser->parse(seqs[i_seq], NULL, NULL);
        
        // only keep suboptimal base pairs
        // the maximum % change in free energy from the lowest free energy structure
        float min_score = seq_viterbi * (1 - max_energy_diff); // max_score must be positive
        // int min_score = seq_viterbi - 1000;
        
        int seq_len = i_seq==0? seq1_len: seq2_len;
        for (int j=0; j<seq_len-1; j++) {
            vector<unordered_map<int, LFState>*> seq_in{&cky_parser->bestH[j], &cky_parser->bestP[j], &cky_parser->bestM[j], &cky_parser->bestM2[j], &cky_parser->bestMulti[j]};
            vector<unordered_map<int, LFState>*> seq_out{&cky_parser->bestH_beta[j], &cky_parser->bestP_beta[j], &cky_parser->bestM_beta[j], &cky_parser->bestM2_beta[j], &cky_parser->bestMulti_beta[j]};
            
            vector<unordered_map<int, int>*> seq_out_saved;
            if (i_seq==0)
                seq_out_saved = {&seq1_out_H[j+1], &seq1_out_P[j+1], &seq1_out_M[j+1], &seq1_out_M2[j+1], &seq1_out_Multi[j+1]};
            else
                seq_out_saved = {&seq2_out_H[j+1], &seq2_out_P[j+1], &seq2_out_M[j+1], &seq2_out_M2[j+1], &seq2_out_Multi[j+1]};

            for (int k=0; k < 5; k++){
                set<int> valid_pos; 
                unordered_map<int, LFState>& beamins = *seq_in[k];
                unordered_map<int, LFState>& beamout = *seq_out[k];
                for(auto& item : beamins) {
                    int i = item.first;
                    int in_score = item.second.score;
                    int out_score = beamout[i].score;

                    if (in_score == VALUE_MIN || out_score == VALUE_MIN)
                        continue;

                    if (in_score + out_score >= min_score) 
                        valid_pos.insert(i);
                }

                for (auto i : valid_pos)
                    (*seq_out_saved[k])[i+1] = beamout[i].score;
            }
        }
    }
    delete cky_parser;
}


#ifdef multilign
void BeamSankoffParser::parse(const vector<string> &seqs, bool limited, const set<pair<int, int>> &allowed_pairs, vector<pair<int, int>> &out_pairs, int num_pairs){
#else
void BeamSankoffParser::parse(const vector<string> &seqs){
#endif
    // test google dense hash map 
    // dense_hash_map<pair<int, int>, int, pair_hash, eqstr> test;
    // pair<int, int> test_key;
    // test.set_empty_key(test_key);
    // test[make_pair(1, 1)] = 2;
    // cout << "test google dense hash map : " << test[make_pair(1, 1)] << endl;

    struct timeval parse_starttime, parse_endtime;
    gettimeofday(&parse_starttime, NULL);
    
    // allocate space and pre-computation
    seq1_len = seqs[0].size() + 1;
    seq2_len = seqs[1].size() + 1;
    // cout << "sequence length in BeamSankoffParser parse: " << seq1_len << " " << seq2_len << endl;
    prepare(seqs);

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    printf("seqs %d %d only pre-computation time: %f seconds.\n", seq1_len, seq2_len, parse_elapsed_time);

    gettimeofday(&parse_starttime, NULL);

    SeqObject* seq1 = &sequences[0];
    SeqObject* seq2 = &sequences[1];

    // preprocess: get next posititon paired with i after j
    {    
        for (int i=0; i<num_seqs; i++){
            // cout << i << endl;
            SeqObject* seq = &sequences[i];
            // seq.next_pair.resize(NOTON);
            int seq_len = seq->seq_len;
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                seq->next_pair[nuci].resize(seq_len, -1);
                int next = -1;
                for (int nucj = seq_len-1; nucj>=1; --nucj) {
                    seq->next_pair[nuci][nucj] = next;
                    // cout << i << " " << nuci << " " << nucj << "  " << seq->next_pair[nuci][nucj] << endl;
                    if (_allowed_pairs[nuci][seq->nucs[nucj]]) next = nucj;
                }
            }
        }
    }
    // cout << "test next pair" << endl;

// #ifdef SPECIAL_HP
// #ifdef lv
    for (int i=0; i<num_seqs; i++){
        SeqObject* seq = &sequences[i];
        v_init_tetra_hex_tri(seq->raw_seq, seq->seq_len, seq->if_tetraloops, seq->if_hexaloops, seq->if_triloops);
    }
// #else
//     if (is_verbose)
//         v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
// #endif
// #endif

    // start CKY decoding
// #ifdef lv

#ifdef dynalign
    bestC[0][0].set(0, 0, 0, 0, 0, MANNER_NONE, MANNER_NONE, 0.0, 0.0, MANNER_ALN, MANNER_ALN);
    bestC[1][1].set(1, 0, 0, -v_score_external_unpaired(0, 0), 0, MANNER_NONE, MANNER_C_eq_C_plus_U_INS1, 1.0, weight, MANNER_ALN, MANNER_INS1);
    bestC[1][0].set(0, 0, 0, 0, -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_INS2, 1.0, weight, MANNER_ALN, MANNER_INS2);
    bestC[2][1].set(1, 0, 0, -v_score_external_unpaired(0, 0), -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_ALN, 0.0, 0.0, MANNER_ALN, MANNER_ALN);

#else
    bestC[0][0].alnobj.set(0, 0, 0, 0, 0, MANNER_NONE, MANNER_NONE, xlog(1.0), weight*xlog(1.0), MANNER_ALN, MANNER_ALN);

    if (0 >= hmmalign.low_bounds[1] && 0 <= hmmalign.up_bounds[1]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_INS1, 1, 0, true);
        bestC[1][1].ins1obj.set(1, 0, 0, -v_score_external_unpaired(0, 0), 0, MANNER_NONE, MANNER_C_eq_C_plus_U_INS1, alignscore, weight*alignscore, MANNER_ALN, MANNER_INS1);
    }
    if (1 >= hmmalign.low_bounds[0] && 1 <= hmmalign.up_bounds[0]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_INS2, 0, 1, true);
        bestC[1][0].ins2obj.set(0, 0, 0, 0, -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_INS2, alignscore, weight*alignscore, MANNER_ALN, MANNER_INS2);
    }
    if (1 >= hmmalign.low_bounds[1] && 1 <= hmmalign.up_bounds[1]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_ALN, 1, 1, true);
        bestC[2][1].alnobj.set(1, 0, 0, -v_score_external_unpaired(0, 0), -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_ALN, alignscore, weight*alignscore, MANNER_ALN, MANNER_ALN);
    }
#endif


    // from left to right
    // processMem_t mem = GetProcessMemory();
    // cout << "VmPeak: " << mem.VmPeak << endl;
    // float pseudo1=0.0, pseudo2=0.0;
    for(int s = 1; s < seq1_len + seq2_len - 1; ++s) {
        // if (s > 2700) verbose = true;

        // mem = GetProcessMemory();
        // cout << "s: " << s << " VmPeak: " << mem.VmPeak << endl;
        if (s % 100 == 1) cout << "s: " << s << endl;
        // if (s >= 199) verbose = true;

        unordered_map<int, State3>& beamH = bestH[s];
        unordered_map<int, State3>& beamP = bestP[s];
        unordered_map<int, State3>& beamMulti = bestMulti[s];
        unordered_map<int, State>& beamM2 = bestM2[s];
        unordered_map<int, State>& beamM = bestM[s];

        unordered_map<int, State3>& beamC = bestC[s];
      
        // hairpin candidates
        // for nucj push H(j, j_next)
        if (verbose) cout << "push H" << endl;
        for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
            int j2 = s - j1;
            if (j1 < 1 || j2 < 1 || j2 >= seq2->seq_len - 1) continue; // boundary case

            // hmm constraint, left bracket
            if (j2 > hmmalign.up_bounds[j1] || j2 < hmmalign.low_bounds[j1]) continue;

            int j1next = j1, j2next;
            pair<int, int> newscore1, newscore2;
            while (j1next != -1) {
                newscore1 = hairpinScore(j1, j1next, seq1);
                j1next = newscore1.first;
                if (j1next == -1) break;

                // single seq folding 
                if (seq1_out_H[j1next].find(j1) == seq1_out_H[j1next].end()) continue;

                j2next = j2;
                while (j2next != -1) {
                    newscore2 = hairpinScore(j2, j2next, seq2);
                    j2next = newscore2.first;
                    if (j2next == -1) break;

                    // single seq folding subopt
                    if (seq2_out_H[j2next].find(j2) == seq2_out_H[j2next].end()) continue;

                    // hmm constraint, right bracket   
                    if (j2next > hmmalign.up_bounds[j1next] || j2next < hmmalign.low_bounds[j1next]) continue;

                    if (verbose) cout << "hairpin candidates: " << j1 << " " << j1next << " " << j2 << " " << j2next << endl;

#ifdef dynalign
                    aln_value_type alignscore = abs(j1next - j1 - j2next + j2); // do not apply alignment constraint on right side
#else
                    float alignscore = get_hmm_score(j1, j1next, j2, j2next, 2);
#endif
                    if (alignscore > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1next+j2next][get_keys(j1next, j1, j2)].alnobj,
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_ALN, 
                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
#ifndef dynalign
                    float ins1score = get_hmm_score(j1, j1next, j2+1, j2next-1, 0);
                    if (ins1score > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1next+j2next][get_keys(j1next, j1, j2)].ins1obj,
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_INS1, 
                                        ins1score, MANNER_INS1, MANNER_INS1, weight, verbose);

                    float ins2score = get_hmm_score(j1+1, j1next-1, j2, j2next, 1);
                    if (ins2score > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1next+j2next][get_keys(j1next, j1, j2)].ins2obj,
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_INS2, 
                                        ins2score, MANNER_INS2, MANNER_INS2, weight, verbose);
#endif
                }
            }
        }

        // beam of H
        if (verbose) cout << "beam of H" << endl;
        {
            // for every state h in H[j]
            //   1. extend h(i, j) to h(i, jnext)
            //   2. generate p(i, j)
                
            // cout << "beam of H: " << beamH.size() << endl;
            if (beam > 0 && beamH.size() > beam) beam_prune(beamH, s, seq1_out_H, seq2_out_H);
            // cout << "beam of H after prune: " << beamH.size() << endl;

            for (auto &item : beamH) {
                for (int m=0; m<3; m++){
#ifdef dynalign
                    if (m != 2) continue; // hairpin must be aligned
#endif

                    State state;
                    switch (m)
                    {
                        case 0:
                            state = item.second.ins1obj;
                            break;
                        case 1:
                            state = item.second.ins2obj;
                            break;
                        case 2:
                            state = item.second.alnobj;
                            break;
                        
                        default:
                            break;
                    }
                    if (state.manner == MANNER_NONE) continue;

                    // State state = item.second;
                    int j1 = state.j1;
                    int i1 = state.i1;
                    int i2 = state.i2;
                    int j2 = s - j1;

                    // assert (j1 > 0); // lisiz TODO: DEBUG
                    // assert (i1 > 0); // lisiz TODO: DEBUG
                    // assert (i2 > 0); // lisiz TODO: DEBUG
                        
                    // 2. generate p(i, j)
                    // single seq folding 
                    {
                        if (seq1_out_P[j1].find(i1) != seq1_out_P[j1].end() && seq2_out_P[j2].find(i2) != seq2_out_P[j2].end()) {
#ifdef multilign
                            if (!limited || allowed_pairs.find(make_pair(i1, j1)) != allowed_pairs.end()) {
#else
                            {
#endif
                                if (verbose) cout << "H to P: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl; 

                                switch(m){
                                    case 0:
                                        update_if_better(i1, j1, i2, j2,
                                                        beamP[item.first].ins1obj,  // beamPINS1[get_keys(j1, i1, i2)], 
                                                        state.seq1foldscore, state.seq2foldscore, state.manner, MANNER_HAIRPIN_INS1, 
                                                        state.alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                        break;
                                    case 1:
                                        update_if_better(i1, j1, i2, j2,
                                                        beamP[item.first].ins2obj, // beamPINS2[get_keys(j1, i1, i2)],
                                                        state.seq1foldscore, state.seq2foldscore, state.manner, MANNER_HAIRPIN_INS2, 
                                                        state.alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                        break;
                                    case 2:
                                        update_if_better(i1, j1, i2, j2,
                                                        beamP[item.first].alnobj, // beamPALN[get_keys(j1, i1, i2)],
                                                        state.seq1foldscore, state.seq2foldscore, state.manner, MANNER_HAIRPIN_ALN, 
                                                        state.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                        break;
                                    default:
                                        cout << "error H->P" << endl;
                                        assert(false);
                                }
                            }
                        }
                    } // H->P(i, j)

                    //   1. extend h(i, j) to h(i, jnext)
//                     { 
//                         pair<int, int> newscore1 = hairpinScore(i1, j1, seq1);
//                         pair<int, int> newscore2 = hairpinScore(i2, j2, seq2);

//                         int j1next = newscore1.first;
//                         int j2next = newscore2.first;   

//                         // h(i1, j1; i2, j2; si, sj) -> h(i1, j1next; i2, j2next; si, sjnext)
//                         // sjnext = aln/ins1/ins2  
//                         // if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
//                         if (j1next != -1 && j2next != -1) {
//                             if (seq1_out_H[j1next].find(i1) != seq1_out_H[j1next].end() && seq2_out_H[j2next].find(i2) != seq2_out_H[j2next].end()) {

//                                 if (verbose) cout << "hairpin candidates jump: " << i1 << " " << j1 << " " << i2 << " " << j2  << " to "  << j1next << " " << j2next << " " << newscore1.second << " " << newscore2.second << endl; 

// #ifdef dynalign
//                                 aln_value_type alignscore = abs(j1next - i1 - j2next + i2);
// #else
//                                 float alignscore = get_hmm_score(i1, j1next, i2, j2next, 2, true);
// #endif
//                                 update_if_better(i1, j1next, i2, j2next,
//                                                 bestH[j1next+j2next][get_keys(j1next, i1, i2)].alnobj,
//                                                 newscore1.second, newscore2.second, 
//                                                 MANNER_NONE, MANNER_H_ALN, 
//                                                 alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
// #ifndef dynalign
//                                 alignscore = get_hmm_score(i1, j1next, i2+1, j2next-1, 0, true);
//                                 update_if_better(i1, j1next, i2, j2next,
//                                                 bestH[j1next+j2next][get_keys(j1next, i1, i2)].ins1obj,
//                                                 newscore1.second, newscore2.second, 
//                                                 MANNER_NONE, MANNER_H_INS1,
//                                                 alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                            
//                                 alignscore = get_hmm_score(i1+1, j1next-1, i2, j2next, 1, true);
//                                 update_if_better(i1, j1next, i2, j2next,
//                                                 bestH[j1next+j2next][get_keys(j1next, i1, i2)].ins2obj,
//                                                 newscore1.second, newscore2.second, 
//                                                 MANNER_NONE, MANNER_H_INS2, 
//                                                 alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
//                             }
// #endif
//                         }

//                         // h(i1, j1; i2, j2; si, sj) -> h(i1, j1next; i2, j2; si, sjnext)
//                         // sjnext = aln/ins1/ins2
//                         // if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) {
//                         if (j1next != -1) {
//                             if (seq1_out_H[j1next].find(i1) != seq1_out_H[j1next].end()) {
//                                 if (verbose) cout << "hairpin candidates jump: " << i1 << " " << j1 << " " << i2 << " " << j2  << " to "  << j1next << " " << j2 << " " << newscore1.second << " " << newscore2.second << endl; 
// #ifdef dynalign
//                                 aln_value_type alignscore = abs(j1next - i1 - j2 + i2);
//     #else
//                                 float alignscore = get_hmm_score(i1, j1next, i2, j2, 2, true);
//     #endif
//                                 update_if_better(i1, j1next, i2, j2,
//                                                 bestH[j1next+j2][get_keys(j1next, i1, i2)].alnobj,
//                                                 newscore1.second, state.seq2foldscore, 
//                                                 MANNER_NONE, MANNER_H_ALN, 
//                                                 alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
// #ifndef dynalign
//                                 alignscore = get_hmm_score(i1, j1next, i2+1, j2-1, 0, true);
//                                 update_if_better(i1, j1next, i2, j2,
//                                                 bestH[j1next+j2][get_keys(j1next, i1, i2)].ins1obj,
//                                                 newscore1.second, state.seq2foldscore, 
//                                                 MANNER_NONE, MANNER_H_INS1, 
//                                                 alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                
//                                 alignscore = get_hmm_score(i1+1, j1next-1, i2, j2, 1, true);
//                                 update_if_better(i1, j1next, i2, j2,
//                                                 bestH[j1next+j2][get_keys(j1next, i1, i2)].ins2obj,
//                                                 newscore1.second, state.seq2foldscore, 
//                                                 MANNER_NONE, MANNER_H_INS2, 
//                                                 alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
//                             }
// #endif
//                         }

//                         // h(i1, j1; i2, j2; si, sj) -> h(i1, j1; i2, j2next; si, sjnext)
//                         // sjnext = aln/ins1/ins2
//                         // if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) {
//                         if (j2next != -1) {
//                             if (seq2_out_H[j2next].find(i2) != seq2_out_H[j2next].end()) {
//                                 if (verbose) cout << "hairpin candidates jump: " << i1 << " " << j1 << " " << i2 << " " << j2  << " to "  << j1 << " " << j2next << " " << newscore1.second << " " << newscore2.second << endl; 
                            
// #ifdef dynalign
//                                 aln_value_type alignscore = abs(j1 - i1 - j2next + i2);
// #else
//                                 float alignscore = get_hmm_score(i1, j1, i2, j2next, 2, true);
// #endif
//                                 update_if_better(i1, j1, i2, j2next, 
//                                                 bestH[j1+j2next][get_keys(j1, i1, i2)].alnobj,
//                                                 state.seq1foldscore, newscore2.second, 
//                                                 MANNER_NONE, MANNER_H_ALN, 
//                                              alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
// #ifndef dynalign
//                                 alignscore = get_hmm_score(i1, j1, i2+1, j2next-1, 0, true);
//                                 update_if_better(i1, j1, i2, j2next, 
//                                                 bestH[j1+j2next][get_keys(j1, i1, i2)].ins1obj,
//                                                 state.seq1foldscore, newscore2.second, 
//                                                 MANNER_NONE, MANNER_H_INS1, 
//                                                 alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);

//                                 alignscore = get_hmm_score(i1+1, j1-1, i2, j2next, 1, true);
//                                 update_if_better(i1, j1, i2, j2next, 
//                                                 bestH[j1+j2next][get_keys(j1, i1, i2)].ins2obj,
//                                                 state.seq1foldscore, newscore2.second, 
//                                                 MANNER_NONE, MANNER_H_INS2, 
//                                                 alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);

//                             }
// #endif
//                         } 
//                     } // h(i, j) -> h(i, jnext)
                } // end loop ALN/INS1/INS2 H
            } // end beamH
        } // beam of H

        // beam of Multi
        // for every state in Multi[j]
        //   1. extend (i, j) to (i, jnext)
        //   2. generate P (i, j)
        if (verbose) cout << "beam of Multi" << endl;
        {
            // cout << "beam of Multi: " << beamMulti.size() << endl;
            float threhold;
            // cout << "beam of Multi: " << beamMulti.size() << endl;
            if (beam > 0 && beamMulti.size() > beam) threhold = beam_prune(beamMulti, s, seq1_out_Multi, seq2_out_Multi);
            // cout << "beam of Multi after prune: " << threhold << " " << beamMulti.size() << endl;

            for (auto &item : beamMulti) {
                for (int m=0; m<3; m++){
#ifdef dynalign
                    if (m != 2) continue; // closed pairs for multiloop must be aligned
#endif
                    State state;
                    switch (m)
                    {
                        case 0:
                            state = item.second.ins1obj;
                            break;
                        case 1:
                            state = item.second.ins2obj;
                            break;
                        case 2:
                            state = item.second.alnobj;
                            break;
                        
                        default:
                            break;
                    }
                    if (state.manner == MANNER_NONE) continue;

                    // State state = item.second;
                    int j1 = state.j1;
                    int i1 = state.i1;
                    int i2 = state.i2;
                    int j2 = s - j1; 
                    
                    // assert (j1 > 0); // lisiz TODO: DEBUG
                    // assert (i1 > 0); // lisiz TODO: DEBUG
                    // assert (i2 > 0); // lisiz TODO: DEBUG
                   
                    // 1. extend (i, j) to (i, jnext)
                    if (false) {
                        tuple<int, int, char, int> result1 = multiloopUnpairedScore(i1, j1, seq1, &state.trace1);
                        tuple<int, int, char, int> result2 = multiloopUnpairedScore(i2, j2, seq2, &state.trace2);

                        int j1next = get<0>(result1);
                        int j2next = get<0>(result2);

                        int new_seq1_l2 = get<3>(result1);
                        int new_seq2_l2 = get<3>(result2);

                        if (verbose) cout << "beamMulti jnext: " << m << " "  << i1  <<  " " << j1 << " " << i2 << " " << j2 << " "  << j1next  <<  " " << j2next << endl; 
                        
                        tuple<float, HMMManner, HMMManner, string, string> alnret;
                        // pair<float, HMMManner> alignscore;
                        float alignscore;
                        int key, news;
                       
                        if (m == 2) {

                            if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
#ifdef dynalign
                                // do not apply alignment constraint on right side
                                aln_value_type alignscore = state.alignscore - abs(state.trace1.paddings.l2 - state.trace2.paddings.l2) + abs(new_seq1_l2 - new_seq2_l2);
#else
                                alignscore = get_hmm_score_right(j1, j1next, j2, j2next, 2, 2, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
#endif
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1next, i2, j2next,  bestMulti[j1next+j2next][get_keys(j1next, i1, i2)].alnobj,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                    state.seq1foldscore + get<1>(result1),
                                                    state.seq2foldscore + get<1>(result2),
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                    state.trace1.paddings.l1, get<3>(result1),
                                                    state.trace2.paddings.l1, get<3>(result2),
                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                            }
                            if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) { // TODO: only the first sequence jump towards right

#ifdef dynalign
                                // do not apply alignment constraint on right side
                                aln_value_type alignscore = state.alignscore - abs(state.trace1.paddings.l2 - state.trace2.paddings.l2) + abs(new_seq1_l2-state.trace2.paddings.l2);
#else
                                alignscore = get_hmm_score_right(j1, j1next, j2, j2, 2, 2, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
#endif
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1next, i2, j2,  bestMulti[j1next+j2][get_keys(j1next, i1, i2)].alnobj,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                    state.seq1foldscore + get<1>(result1),
                                                    state.seq2foldscore,
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                    state.trace1.paddings.l1, new_seq1_l2,
                                                    state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                            }
                            if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) { // only the second sequence jump towards right
#ifdef dynalign
                                // do not apply alignment constraint on right side
                                aln_value_type alignscore = state.alignscore - abs(state.trace1.paddings.l2 - state.trace2.paddings.l2) + abs(state.trace1.paddings.l2-new_seq2_l2);
#else
                                alignscore = get_hmm_score_right(j1, j1, j2, j2next, 2, 2, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
#endif
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1, i2, j2next, bestMulti[j1+j2next][get_keys(j1, i1, i2)].alnobj,  // bestMultiALN[j1next+j2next][get_keys(j1next, i1, i2)], 
                                                    state.seq1foldscore,
                                                    state.seq2foldscore + get<1>(result2),
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                    state.trace1.paddings.l1, state.trace1.paddings.l2,
                                                    state.trace2.paddings.l1, new_seq2_l2,
                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                            }
                        }
                                
                        else if (m == 0){

                            if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {

                                alignscore = get_hmm_score_right(j1, j1next, j2, j2next-1, 0, 0, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1next, i2, j2next, bestMulti[j1next+j2next][get_keys(j1next, i1, i2)].ins1obj, // bestMultiINS1[j1next+j2][get_keys(j1next, i1, i2)],
                                                    state.seq1foldscore + get<1>(result1),
                                                    state.seq2foldscore + get<1>(result2),
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                    state.trace1.paddings.l1, new_seq1_l2,
                                                    state.trace2.paddings.l1, new_seq2_l2,
                                                    alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                            }
                            if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) { // TODO: only the first sequence jump towards right

                                alignscore = get_hmm_score_right(j1, j1next, j2, j2-1, 0, 0, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1next, i2, j2, bestMulti[j1next+j2][get_keys(j1next, i1, i2)].ins1obj, // bestMultiINS1[j1next+j2][get_keys(j1next, i1, i2)],
                                                    state.seq1foldscore + get<1>(result1),
                                                    state.seq2foldscore,
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                    state.trace1.paddings.l1, new_seq1_l2,
                                                    state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                    alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                            }
                            if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) { // only the second sequence jump towards right

                                alignscore = get_hmm_score_right(j1, j1, j2, j2next-1, 0, 0, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1, i2, j2next, bestMulti[j1+j2next][get_keys(j1, i1, i2)].ins1obj, // bestMultiINS1[j1next+j2][get_keys(j1next, i1, i2)],
                                                    state.seq1foldscore,
                                                    state.seq2foldscore + get<1>(result2),
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                    state.trace1.paddings.l1, state.trace1.paddings.l2,
                                                    state.trace2.paddings.l1, new_seq2_l2,
                                                    alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                            }
                        }
                                
                        else if (m == 1) {

                            if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {

                                alignscore = get_hmm_score_right(j1, j1next-1, j2, j2next, 1, 1, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1next, i2, j2next, bestMulti[j1next+j2next][get_keys(j1next, i1, i2)].ins2obj, // bestMultiINS1[j1next+j2][get_keys(j1next, i1, i2)],
                                                    state.seq1foldscore + get<1>(result1),
                                                    state.seq2foldscore + get<1>(result2),
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                    state.trace1.paddings.l1, new_seq1_l2,
                                                    state.trace2.paddings.l1, new_seq2_l2,
                                                    alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                            }

                            if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) { // TODO: only the first sequence jump towards right

                                alignscore = get_hmm_score_right(j1, j1next-1, j2, j2, 1, 1, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1next, i2, j2,  bestMulti[j1next+j2][get_keys(j1next, i1, i2)].ins2obj, // bestMultiINS1[j1next+j2][get_keys(j1next, i1, i2)],
                                                    state.seq1foldscore + get<1>(result1),
                                                    state.seq2foldscore,
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                    state.trace1.paddings.l1, new_seq1_l2,
                                                    state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                    alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                            }

                            if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) { // only the second sequence jump towards right

                                alignscore = get_hmm_score_right(j1, j1-1, j2, j2next, 1, 1, true); // true
                                alignscore = xlog_mul(state.alignscore, alignscore);
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(i1, j1, i2, j2next, bestMulti[j1+j2next][get_keys(j1, i1, i2)].ins2obj, // bestMultiINS1[j1next+j2][get_keys(j1next, i1, i2)],
                                                    state.seq1foldscore,
                                                    state.seq2foldscore + get<1>(result2),
                                                    state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                    state.trace1.paddings.l1, state.trace1.paddings.l2,
                                                    state.trace2.paddings.l1, new_seq2_l2,
                                                    alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                            }
                        }
                    } // 1. extend (i, j) to (i, jnext)

                    // 2. generate P (i, j)
                    // single seq folding 
                    if (seq1_out_P[j1].find(i1) != seq1_out_P[j1].end() && seq2_out_P[j2].find(i2) != seq2_out_P[j2].end()) {
#ifdef multilign
                        if (!limited || allowed_pairs.find(make_pair(i1, j1)) != allowed_pairs.end()) {
#else
                        {
#endif
                            if (verbose) cout << "Multi to P: " << m <<  " " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;

                            int score1 = multiloop2Pscore(i1, j1, seq1);
                            int score2 = multiloop2Pscore(i2, j2, seq2);
                            if (m == 2)
                                update_if_better(i1, j1, i2, j2, beamP[item.first].alnobj, // beamPALN[get_keys(j1, i1, i2)],
                                                state.seq1foldscore + score1, 
                                                state.seq2foldscore + score2,
                                                state.manner, MANNER_P_eq_MULTI_ALN,
                                                state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                            else if (m == 0)
                                update_if_better(i1, j1, i2, j2, beamP[item.first].ins1obj, // beamPINS1[get_keys(j1, i1, i2)],
                                                state.seq1foldscore + score1, 
                                                state.seq2foldscore + score2,
                                                state.manner, MANNER_P_eq_MULTI_INS1,
                                                state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                            else if (m == 1)
                                update_if_better(i1, j1, i2, j2, beamP[item.first].ins2obj, // beamPINS2[get_keys(j1, i1, i2)], 
                                                state.seq1foldscore + score1, 
                                                state.seq2foldscore + score2,
                                                state.manner, MANNER_P_eq_MULTI_INS2,
                                                state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);

                        }
                    } // 2. generate P (i, j)
                } // beam ALN/INS1/INS2
            } // end beamMulti
        } // end beam of Multi

        // beam of P
        // for every state in P[j]
        //   1. generate new helix/bulge
        //   2. M = P
        //   3. M2 = M + P
        //   4. C = C + P
        if (verbose) cout << "beam of P" << endl;
        {
            // cout << "beam of P: "  << beamP.size() << endl;
            // if (s == 383) verbose=true;
            float threhold;
            if (beam > 0 && beamP.size() > beam) threhold = beam_prune(beamP, s, seq1_out_P, seq2_out_P);

#ifdef is_cube_pruning
            bool use_cube_pruning = beam > MIN_CUBE_PRUNING_SIZE && beamP.size() > MIN_CUBE_PRUNING_SIZE;
            sort_keys(beamP, keys); // sort P states
#else
            bool use_cube_pruning = false;
#endif  

            for (auto &item : beamP) {
                State3& state3 = item.second;
                
                // ALN state
                {
                    State& alnstate = state3.alnobj;
                    if (alnstate.manner != MANNER_NONE){
                        int j1 = alnstate.j1;
                        int i1 = alnstate.i1;
                        int i2 = alnstate.i2;
                        int j2 = s - j1;
                        
                        // assert (j1 > 0); // lisiz TODO: DEBUG
                        // assert (i1 > 0); // lisiz TODO: DEBUG
                        // assert (i2 > 0); // lisiz TODO: DEBUG

                        int m = 2;

                        int nuci1 = seq1->nucs[i1];
                        int nuci1_1 = seq1->nucs[i1 - 1];
                        int nucj1 = seq1->nucs[j1];
                        int nucj1p1 = seq1->nucs[j1 + 1];

                        int nuci2 = seq2->nucs[i2];
                        int nuci2_1 = seq2->nucs[i2 - 1];
                        int nucj2 = seq2->nucs[j2];
                        int nucj2p1 = seq2->nucs[j2 + 1];

                        // 1. generate new helix / single_branch
                        // new state is of shape p..i..j..q
                        // Note: p >= 1
                        {
                            // for (int p1 = i1 - 1; p1 >= max(i1 - SINGLE_MAX_LEN, 1); --p1) {
                            int min_p1 = max(i1 - MAX_LOOP_LEN - 1, 1);
                            for (int p1 = i1-1; p1 >= min_p1; --p1) {
                                int nucp1 = seq1->nucs[p1];
                                int nucp1p1 = seq1->nucs[p1 + 1];
                                int q1 = seq1->next_pair[nucp1][j1];

                                // int i1_p1 = i1 - p1;
                                // while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                                while (q1 != -1 &&  (q1 - j1 - 1 <= MAX_LOOP_LEN)) {
                                    // single seq folding
                                    if (seq1_out_P[q1].find(p1) == seq1_out_P[q1].end()) {
                                        q1 = seq1->next_pair[nucp1][q1];
                                        continue;
                                    }

#ifdef multilign
                                    if (limited && allowed_pairs.find(make_pair(p1, q1)) == allowed_pairs.end()) {
                                        q1 = seq1->next_pair[nucp1][q1];
                                        continue;
                                    }
#endif
                                    int nucq1 = seq1->nucs[q1];
                                    int nucq1_1 = seq1->nucs[q1 - 1];

                                    int p2p1 = P2PScore(p1,q1,i1,j1,nucp1,nucp1p1,nucq1_1,nucq1,nuci1_1,nuci1,nucj1,nucj1p1); 

                                    int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                                    int min_p2 = max(hmmalign.low_bounds[p1], max(1, i2 - MAX_LOOP_LEN - 1));
                                    for (int p2 = max_p2; p2 >= min_p2; --p2) { // seq2 for loop
                                        int nucp2 = seq2->nucs[p2];
                                        int nucp2p1 = seq2->nucs[p2 + 1]; // hzhang: move here
                                        int q2 = seq2->next_pair[nucp2][j2];

                                        // speed up
                                        int i2_p2 = i2 - p2;
                                        if (q2 > hmmalign.up_bounds[q1] || q2 == -1 || (q2 - j2 - 1 > MAX_LOOP_LEN)) continue;
                                        if (q2 < hmmalign.low_bounds[q1])
                                            q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1] - 1]; 

                                        // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                        while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= MAX_LOOP_LEN)) {
                                            // single seq folding
                                            if (seq2_out_P[q2].find(p2) == seq2_out_P[q2].end()) {
                                                q2 = seq2->next_pair[nucp2][q2];
                                                continue;
                                            }

                                            // TODO: redundant calculation
                                            int nucq2 = seq2->nucs[q2];
                                            int nucq2_1 = seq2->nucs[q2 - 1];
                                            int p2p2 = P2PScore(p2,q2,i2,j2,nucp2,nucp2p1,nucq2_1,nucq2,nuci2_1,nuci2,nucj2,nucj2p1);

                                            if (verbose) cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;
                                            
                                            pair<float, HMMManner> pre_alignscore, post_alignscore;
                                            float pre_align_trans, post_align_trans, alignscore;
                                            
                                            // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj)
                                            //  si == sj == ALN/INS1/INS2
                                            {          
#ifdef dynalign                                     
                                                aln_value_type alignscore = ALN_VALUE_MIN;
                                                if (p2>=hmmalign.low_bounds[p1] && p2<=hmmalign.up_bounds[p1] && q2>=hmmalign.low_bounds[q1] && q2<=hmmalign.up_bounds[q1]){
                                                    alignscore = alnstate.alignscore + abs(i1 - p1 - 1) + abs(q1 - j1 - 1);
                                                }
#else
                                                // si == sj == ALN
                                                // align cost = (p, i) + (i, j) + (j, q)
                                                pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                                                
                                                alignscore = xlog_mul(pre_align_trans, alnstate.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);

#endif                                        
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].alnobj, // bestPALN[q1+q2][get_keys(q1, p1, p2)],  // helix, one branch using one state MANNER_SINGLE
                                                                        alnstate.seq1foldscore + p2p1,
                                                                        alnstate.seq2foldscore + p2p2,
                                                                        alnstate.manner, MANNER_SINGLE_ALN, 
                                                                        static_cast<char>(i1 - p1), q1 - j1,
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                            }
                                
#ifdef dynalign
                                            if (p1==(i1-1) && p2==(i2-1) && q1==(j1+1) && q2==(j2+1))
                                            { // ALN -> INS1/INS2, not allow continual inserted base pairs
                                                aln_value_type alignscore = ALN_VALUE_MIN;
                                                if (p2+1>=hmmalign.low_bounds[p1] && p2+1<=hmmalign.up_bounds[p1] && q2-1>=hmmalign.low_bounds[q1] && q2-1<=hmmalign.up_bounds[q1]){
                                                    alignscore = alnstate.alignscore + abs(i1 - p1 - i2 + p2) + abs(q1 - j1 - q2 + j2) + 2;
                                                }
#else
                                            {
                                                // si == sj == INS1
                                                // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                                pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0);
                                                
                                                alignscore = xlog_mul(pre_align_trans, alnstate.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans); 

#endif   
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins1obj, // bestPINS1[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                    alnstate.seq1foldscore + p2p1,
                                                                    alnstate.seq2foldscore + p2p2,
                                                                    alnstate.manner, MANNER_SINGLE_INS1, 
                                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                            }
                                            
                                            
#ifdef dynalign
                                            if (p1==(i1-1) && p2==(i2-1) && q1==(j1+1) && q2==(j2+1))
                                            { // ALN -> INS1/INS2, not allow continual inserted base pairs
                                                aln_value_type alignscore = ALN_VALUE_MIN;
                                                if (p2>=hmmalign.low_bounds[p1+1] && p2<=hmmalign.up_bounds[p1+1] && q2>=hmmalign.low_bounds[q1-1] && q2<=hmmalign.up_bounds[q1-1]){
                                                    alignscore = alnstate.alignscore + abs(i1 - p1 - i2 + p2) + abs(q1 - j1 - q2 + j2) + 2;
                                                }
#else
                                            {
                                                // si == sj == INS2
                                                // align cost = (p, i) + (i, j) + (j, q);
                                                pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, m);
                                                post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, m, 1);
                                                
                                                alignscore = xlog_mul(pre_align_trans, alnstate.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);

#endif 
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins2obj, // bestPINS2[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                    alnstate.seq1foldscore + p2p1,
                                                                    alnstate.seq2foldscore + p2p2,
                                                                    alnstate.manner, MANNER_SINGLE_INS2, 
                                                                    static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                            }
                                            q2 = seq2->next_pair[nucp2][q2];
                                        } // while loop enumerate q2
                                    } // for loop enumerate p2
                                    q1 = seq1->next_pair[nucp1][q1];
                                } // while loop enumerate q1
                            } // for loop enumerate p1
                        }

                        // 2. M = P
                        // accessible pairs in multiloop must be aligned
                        // singe seq folding
                        if (seq1_out_M[j1].find(i1) != seq1_out_M[j1].end() && seq2_out_M[j2].find(i2) != seq2_out_M[j2].end()) 
                        {   
                            if (verbose) cout << "P to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                            int newscore1 = branch_score(i1, j1, seq1);
                            int newscore2 = branch_score(i2, j2, seq2);
                            update_if_better(i1, j1, i2, j2, beamM[item.first], // bestM[key1][newi2][newj2], // beamM[j2][make_pair(i1, i2)], 
                                            alnstate.seq1foldscore + newscore1,
                                            alnstate.seq2foldscore + newscore2, 
                                            alnstate.manner, MANNER_M_eq_P,
                                            alnstate.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                        } // 2. M = P

                        // 3. M2 = M + P 
                        // accessible pairs in multiloop must be aligned
                        if (!use_cube_pruning){
                            int k1 = i1 - 1;
                            int k2 = i2 - 1;

                            if (k1 > 1 && k2 > 1 && !bestM[k1 + k2].empty()) {
                                int newscore1 = branch_score(i1, j1, seq1);
                                int newscore2 = branch_score(i2, j2, seq2); 
#ifndef is_candidate_list
                                for (auto &m : bestM[k1 + k2]) {
                                    State mstate = m.second; 
                                    int newj1 = mstate.j1;
                                    if (newj1 != k1) continue;
                                    int newi1 = mstate.i1;
                                    int newi2 = mstate.i2;

                                    // single seq folding 
                                    if (seq1_out_M2[j1].find(newi1) != seq1_out_M2[j1].end() && seq2_out_M2[j2].find(newi2) != seq2_out_M2[j2].end()) 
                                    { 

                                        if (verbose) cout << "M2=M+P: " << i1 << " " << j1 << " "  << i2 << " " << j2 <<  " " << newi1 << " " << j1 << " "  << newi2 << " " << j2  << endl;
                            
#ifdef dynalign
                                        aln_value_type alignscore = alnstate.alignscore + mstate.alignscore;
#else
                                        float alignscore = xlog_mul(mstate.alignscore, hmmalign.trans_probs[mstate.endHMMstate-1][MANNER_ALN-1]); 
                                        alignscore = xlog_mul(alignscore, alnstate.alignscore);
#endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(newi1, j1, newi2, j2, beamM2[get_keys(j1, newi1, newi2)], // bestM2[get_keys(newi1, j1)][mnewi2][newj2], // beamM2[j2][make_pair(newi1, newi2)], // Note: not i1 i2 but newi1 newi2new
                                                            mstate.seq1foldscore + newscore1 + alnstate.seq1foldscore,
                                                            mstate.seq2foldscore + newscore2 + alnstate.seq2foldscore,
                                                            mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);     
                                    }                     
                                }
#else 
                                // candidate list
                                auto bestM2_iter = beamM2.find(item.first);
                                int M1_score = newscore1 + newscore2 + state.score;
                                if (bestM2_iter==beamM2.end() || M1_score > bestM2_iter->second.score) {
                                    for (auto &m : bestM[k1+k2]) {
                                        State mstate = m.second; 
                                        int newj1 = mstate.j1;
                                        if (newj1 != k1) continue;
                                        int newi1 = mstate.i1;
                                        int newi2 = mstate.i2;
                                        
                                        // eq. to first convert P to M1, then M2/M = M + M1
                                        float alignscore = xlog_mul(mstate.alignscore, hmmalign.trans_probs[mstate.endHMMstate-1][MANNER_ALN-1]); 
                                        alignscore = xlog_mul(alignscore, alnstate.alignscore);

                                        update_if_better(newi1, j1, newi2, j2, beamM2[get_keys(j1, newi1, newi2)], // bestM2[get_keys(newi1, j1)][mnewi2][newj2], // beamM2[j2][make_pair(newi1, newi2)], // Note: not i1 i2 but newi1 newi2new
                                                        mstate.seq1foldscore + newscore1 + alnstate.seq1foldscore,
                                                        mstate.seq2foldscore + newscore2 + alnstate.seq2foldscore,
                                                        mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                                    }
                                }
#endif 
                            }
                        } // 3. M2 = M + P 

                        // 4. C = C + P
                        // external pairs must be aligned
                        {   
                            int k1 = i1 - 1;
                            int k2 = i2 - 1;

                            if (k1 >= 0 && k2 >= 0) {
                                if (bestC[k1+k2].find(k1) == bestC[k1+k2].end()) continue;

                                if (verbose) cout << "C+P: "<< i1 << " " << j1 << " "  << i2 << " " << j2 << endl;

                                int newscore1 = alnstate.seq1foldscore + external_paired_score(k1, j1, seq1);
                                int newscore2 = alnstate.seq2foldscore + external_paired_score(k2, j2, seq2);

                                {
                                    State& prefix_C = bestC[k1+k2][k1].alnobj; // bestC[k1][k2];
                                    if (prefix_C.endHMMstate != HMMMANNER_NONE) {
#ifdef dynalign
                                        int alignscore = alnstate.alignscore + prefix_C.alignscore;
#else
                                        float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[2][2]); 
                                        alignscore = xlog_mul(alignscore, alnstate.alignscore); 
#endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(0, j1, 0, j2, beamC[j1].alnobj,
                                                            prefix_C.seq1foldscore + newscore1,
                                                            prefix_C.seq2foldscore + newscore2, 
                                                            prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                            k1, k2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                    }
                                }
                                {
                                    State& prefix_C = bestC[k1+k2][k1].ins1obj;
                                    if (prefix_C.endHMMstate != HMMMANNER_NONE) {
#ifdef dynalign
                                        int alignscore = alnstate.alignscore + prefix_C.alignscore;
#else
                                        float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[0][2]); 
                                        alignscore = xlog_mul(alignscore, alnstate.alignscore); 
#endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(0, j1, 0, j2, beamC[j1].alnobj,
                                                            prefix_C.seq1foldscore + newscore1,
                                                            prefix_C.seq2foldscore + newscore2, 
                                                            prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                            k1, k2,
                                                            alignscore, MANNER_INS1, MANNER_ALN, weight, verbose);
                                    }
                                }
                                {
                                    State& prefix_C = bestC[k1+k2][k1].ins2obj;
                                    if (prefix_C.endHMMstate != HMMMANNER_NONE) {
#ifdef dynalign
                                        int alignscore = alnstate.alignscore + prefix_C.alignscore;
#else
                                        float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[1][2]); 
                                        alignscore = xlog_mul(alignscore, alnstate.alignscore); 
#endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(0, j1, 0, j2, beamC[j1].alnobj,
                                                            prefix_C.seq1foldscore + newscore1,
                                                            prefix_C.seq2foldscore + newscore2, 
                                                            prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                            k1, k2,
                                                            alignscore, MANNER_INS2, MANNER_ALN, weight, verbose);
                                    }
                                }   
                            }
                        } // 4. C = C + P
                    } 
                } // ALN state

                // INS2 state
                State& ins2state = state3.ins2obj; 
                if (ins2state.manner != MANNER_NONE) {
                    int j1 = ins2state.j1;
                    int i1 = ins2state.i1;
                    int i2 = ins2state.i2;
                    int j2 = s - j1;
                    
                    // assert (j1 > 0); // lisiz TODO: DEBUG
                    // assert (i1 > 0); // lisiz TODO: DEBUG
                    // assert (i2 > 0); // lisiz TODO: DEBUG

                    int m = 1;

                    int nuci1 = seq1->nucs[i1];
                    int nuci1_1 = seq1->nucs[i1 - 1];
                    int nucj1 = seq1->nucs[j1];
                    int nucj1p1 = seq1->nucs[j1 + 1];
                    int nuci2 = seq2->nucs[i2];
                    int nuci2_1 = seq2->nucs[i2 - 1];
                    int nucj2 = seq2->nucs[j2];
                    int nucj2p1 = seq2->nucs[j2 + 1];

                    // fixed p1==i1 q1==j1
                    int p1 = i1;
                    int q1 = j1;

                    int nucp1 = seq1->nucs[p1];
                    int nucp1p1 = seq1->nucs[p1 + 1];
                    int nucq1 = seq1->nucs[q1];
                    int nucq1_1 = seq1->nucs[q1 - 1];
                    int p2p1 = 0; 

                    int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                    int min_p2 = max(hmmalign.low_bounds[p1], max(1, i2 - MAX_LOOP_LEN - 1));

                    for (int p2 = max_p2; p2 >= min_p2; --p2) {
                        int nucp2 = seq2->nucs[p2];
                        int nucp2p1 = seq2->nucs[p2 + 1];
                        int q2 = seq2->next_pair[nucp2][j2];

                        // speed up
                        // int i2_p2 = i2 - p2;
                        if (q2 > hmmalign.up_bounds[q1] || q2 == -1 || (q2 - j2 - 1 > MAX_LOOP_LEN)) continue;
                        if (q2 < hmmalign.low_bounds[q1])
                            q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1] - 1]; 

                        // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                        while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= MAX_LOOP_LEN)) {
                            // single seq folding
                            if (seq2_out_P[q2].find(p2) == seq2_out_P[q2].end()) {
                                q2 = seq2->next_pair[nucp2][q2];
                                continue;
                            }

                            // TODO: redundant calculation
                            int nucq2 = seq2->nucs[q2];
                            int nucq2_1 = seq2->nucs[q2 - 1];
                            int p2p2 = P2PScore(p2,q2,i2,j2,nucp2,nucp2p1,nucq2_1,nucq2,nuci2_1,nuci2,nucj2,nucj2p1);

                            if (verbose) cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;
                            
                            pair<float, HMMManner> pre_alignscore, post_alignscore;
                            float pre_align_trans, post_align_trans, alignscore;

                            // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj)
                            // si == sj == ALN/INS1/INS2    
#ifdef dynalign                                     
                            aln_value_type alignscore = ALN_VALUE_MIN;
                            if (p2>=hmmalign.low_bounds[p1] && p2<=hmmalign.up_bounds[p1] && q2>=hmmalign.low_bounds[q1] && q2<=hmmalign.up_bounds[q1]){
                                alignscore = ins2state.alignscore + abs(i2 - p2 - 1) + abs(q2 - j2 - 1);
                            }
#else
                            // si == sj == ALN
                            // align cost = (p, i) + (i, j) + (j, q)
                            pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                            alignscore = xlog_mul(pre_align_trans, ins2state.alignscore);
                            alignscore = xlog_mul(alignscore, post_align_trans);
#endif                                        
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].alnobj, // bestPALN[q1+q2][get_keys(q1, p1, p2)],  // helix, one branch using one state MANNER_SINGLE
                                                    ins2state.seq1foldscore + p2p1,
                                                    ins2state.seq2foldscore + p2p2,
                                                    ins2state.manner, MANNER_SINGLE_ALN, 
                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                            // si == sj == INS1
                            // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                            pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0);
                            alignscore = xlog_mul(pre_align_trans, ins2state.alignscore);
                            alignscore = xlog_mul(alignscore, post_align_trans); 
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins1obj, // bestPINS1[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                ins2state.seq1foldscore + p2p1,
                                                ins2state.seq2foldscore + p2p2,
                                                ins2state.manner, MANNER_SINGLE_INS1, 
                                                static_cast<char>(i1 - p1), q1 - j1,
                                                static_cast<char>(i2 - p2), q2 - j2,
                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);

                            // si == sj == INS2
                            // align cost = (p, i) + (i, j) + (j, q);
                            pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, m);
                            post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, m, 1);
                            alignscore = xlog_mul(pre_align_trans, ins2state.alignscore);
                            alignscore = xlog_mul(alignscore, post_align_trans);
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins2obj, // bestPINS2[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                ins2state.seq1foldscore + p2p1,
                                                ins2state.seq2foldscore + p2p2,
                                                ins2state.manner, MANNER_SINGLE_INS2, 
                                                static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                static_cast<char>(i2 - p2), q2 - j2,
                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                        
                            q2 = seq2->next_pair[nucp2][q2];
                        } // while loop enumerate q2
                    } // for loop enumerate p2
                }
                    
                // INS1 state
                State& ins1state = state3.ins1obj;
                if (ins1state.manner != MANNER_NONE) {
                    int j1 = ins1state.j1;
                    int i1 = ins1state.i1;
                    int i2 = ins1state.i2;
                    int j2 = s - j1;
                    
                    // assert (j1 > 0); // lisiz TODO: DEBUG
                    // assert (i1 > 0); // lisiz TODO: DEBUG
                    // assert (i2 > 0); // lisiz TODO: DEBUG

                    int m = 0;

                    int nuci1 = seq1->nucs[i1];
                    int nuci1_1 = seq1->nucs[i1 - 1];
                    int nucj1 = seq1->nucs[j1];
                    int nucj1p1 = seq1->nucs[j1 + 1];
                    int nuci2 = seq2->nucs[i2];
                    int nuci2_1 = seq2->nucs[i2 - 1];
                    int nucj2 = seq2->nucs[j2];
                    int nucj2p1 = seq2->nucs[j2 + 1];

                    // fixed p2 q2
                    int p2 = i2;
                    int q2 = j2;

                    int nucp2 = seq2->nucs[p2];
                    int nucp2p1 = seq2->nucs[p2 + 1];
                    int nucq2 = seq2->nucs[q2];
                    int nucq2_1 = seq2->nucs[q2 - 1];
                    int p2p2 = 0;

                    int min_p1 = max(i1 - MAX_LOOP_LEN - 1, 1); // SINGLE_MAX_LEN
                    for (int p1 = i1-1; p1 >= min_p1; --p1) {
                        if (p2 < hmmalign.low_bounds[p1] || p2 > hmmalign.up_bounds[p1]) continue; // Note.

                        int nucp1 = seq1->nucs[p1];
                        int nucp1p1 = seq1->nucs[p1 + 1];
                        int q1 = seq1->next_pair[nucp1][j1];

                        // int i1_p1 = i1 - p1;
                        // while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                        while (q1 != -1 && (q1 - j1 - 1 <= MAX_LOOP_LEN)) {
                            if (q2 < hmmalign.low_bounds[q1] || q2 > hmmalign.up_bounds[q1]) { // Note.
                                q1 = seq1->next_pair[nucp1][q1];
                                continue;
                            }

                            // single seq folding
                            if (seq1_out_P[q1].find(p1) == seq1_out_P[q1].end()) {
                                q1 = seq1->next_pair[nucp1][q1];
                                continue;
                            }
#ifdef multilign
                            if (limited && allowed_pairs.find(make_pair(p1, q1)) == allowed_pairs.end()) {
                                q1 = seq1->next_pair[nucp1][q1];
                                continue;
                            }
#endif
                            int nucq1 = seq1->nucs[q1];
                            int nucq1_1 = seq1->nucs[q1 - 1];
                            int p2p1 = P2PScore(p1,q1,i1,j1,nucp1,nucp1p1,nucq1_1,nucq1,nuci1_1,nuci1,nucj1,nucj1p1); 
                            
                            if (verbose) cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;
                                    
                            pair<float, HMMManner> pre_alignscore, post_alignscore;
                            float pre_align_trans, post_align_trans, alignscore;
                            // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj), si == sj == ALN/INS1/INS2
                            {          
#ifdef dynalign                                     
                                aln_value_type alignscore = ALN_VALUE_MIN;
                                if (p2>=hmmalign.low_bounds[p1] && p2<=hmmalign.up_bounds[p1] && q2>=hmmalign.low_bounds[q1] && q2<=hmmalign.up_bounds[q1]){
                                    alignscore = ins1state.alignscore + abs(i1 - p1 - 1) + abs(q1 - j1 - 1);
                                }
#else
                                // si == sj == ALN
                                // align cost = (p, i) + (i, j) + (j, q)
                                pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                                
                                alignscore = xlog_mul(pre_align_trans, ins1state.alignscore);
                                alignscore = xlog_mul(alignscore, post_align_trans);

#endif                                        
                                if (alignscore > LOG_OF_ZERO)
                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].alnobj, // bestPALN[q1+q2][get_keys(q1, p1, p2)],  // helix, one branch using one state MANNER_SINGLE
                                                        ins1state.seq1foldscore + p2p1,
                                                        ins1state.seq2foldscore + p2p2,
                                                        ins1state.manner, MANNER_SINGLE_ALN, 
                                                        static_cast<char>(i1 - p1), q1 - j1,
                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                            }
                        

                            // si == sj == INS1
                            // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                            pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0);
                            alignscore = xlog_mul(pre_align_trans, ins1state.alignscore);
                            alignscore = xlog_mul(alignscore, post_align_trans);   
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins1obj, // bestPINS1[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                ins1state.seq1foldscore + p2p1,
                                                ins1state.seq2foldscore + p2p2,
                                                ins1state.manner, MANNER_SINGLE_INS1, 
                                                static_cast<char>(i1 - p1), q1 - j1,
                                                static_cast<char>(i2 - p2), q2 - j2,
                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                    
                            // si == sj == INS2
                            // align cost = (p, i) + (i, j) + (j, q);
                            pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, m);
                            post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, m, 1);
                            alignscore = xlog_mul(pre_align_trans, ins1state.alignscore);
                            alignscore = xlog_mul(alignscore, post_align_trans);
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins2obj, // bestPINS2[q1+q2][get_keys(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                ins1state.seq1foldscore + p2p1,
                                                ins1state.seq2foldscore + p2p2,
                                                ins1state.manner, MANNER_SINGLE_INS2, 
                                                static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                static_cast<char>(i2 - p2), q2 - j2,
                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                            
                            q1 = seq1->next_pair[nucp1][q1];
                        } // while loop enumerate q1
                    }
                } // INS1 state
            } // beamP
        }

#ifdef is_cube_pruning
                    // 3. M2 = M + P with cube pruning
                    if (use_cube_pruning) {
                        cout << "cube pruning" << endl; // debug
                        vector<int> valid_Ps;
                        vector<tuple<float, int, int>> M1_scores;
                        // sort_keys(beamP, keys);
                        for (auto &item : keys) { // sorted P state
                            int key = item.first;
                            State &state = item.second;
                            int i1 = state.i1;
                            int i2 = state.i2;
                            int j1 = state.j1;
                            int j2 = s - j1;

                            int k1 = i1 - 1;
                            int k2 = i2 - 1;
                            int k = i1 + i2 - 2;

                            // debug
                            // cout << "cube pruning: " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.score << endl;

                            // group candidate Ps
                            if (k1 > 1 && k2 > 1 && !sorted_bestM[k1+k2].empty()){ // sorted M state
                                // if (bestM[k1+k2].size() != sorted_bestM[k1+k2].size()) {
                                //     cout << "error: " << bestM[k1+k2].size() << " " << sorted_bestM[k1+k2].size() << endl;
                                // }
                                // assert(bestM[k1+k2].size() == sorted_bestM[k1+k2].size());

                                int newscore1 = branch_score(i1, j1, seq1);
                                int newscore2 = branch_score(i2, j2, seq2); 
                                float M1_score = newscore1 + newscore2 + state.score;

#ifndef is_candidate_list
                                valid_Ps.push_back(key);
                                M1_scores.push_back(make_tuple(M1_score, newscore1+state.seq1foldscore, newscore2+state.seq2foldscore));
#else 
                                auto bestM2_iter = beamM2.find(key); // TODO: why compare with M2
                                if (bestM2_iter == beamM2.end() || M1_score > bestM2_iter->second.score) {
                                    valid_Ps.push_back(key);
                                    M1_scores.push_back(make_tuple(M1_score, newscore1+state.seq1foldscore, newscore2+state.seq2foldscore));
                                }
#endif
                            }
                        }
                        // build max heap
                        // heap is of form (heuristic score, (index of i in valid_Ps, index of M in bestM[i-1]))
                        vector<pair<value_type, pair<int, int>>> heap;
                        for (int p = 0; p < valid_Ps.size(); ++p) {
                            int key = valid_Ps[p];
                            State &state = beamP[key].alnobj;
                            int i1 = state.i1;
                            int i2 = state.i2;
                            int j1 = state.j1;
                            int k = i1 + i2 - 2;

                            // cout << "build max heap0: " << i1 << " " << j1 << " " << i2 << " " << s-j1 << sorted_bestM[k].size() << " " <<M1_scores.size() << endl;

                            State &mstate = bestM[k][sorted_bestM[k][0].second]; // the best M state
                            float transit_prob = hmmalign.trans_probs[mstate.endHMMstate-1][MANNER_ALN-1];
                            float alignscore = xlog_mul(mstate.alignscore, transit_prob);
                            float sum_score = get<0>(M1_scores[p]) + mstate.seq1foldscore + mstate.seq2foldscore + weight * alignscore;

                            cout << "build max heap: " << i1 << " " << j1 << " " << i2 << " " << s-j1 << " " << sum_score << " " << p << " " <<valid_Ps.size() << endl;

                            heap.push_back(make_pair(sum_score, make_pair(p, 0)));
                            push_heap(heap.begin(), heap.end());
                        }
                        // start cube pruning
                        // stop after beam size M2 states being filled
                        int filled = 0;
                        // exit when filled >= beam and current score < prev score
                        float prev_score = VALUE_FMIN;
                        float current_score = VALUE_FMIN;
                        while ((filled < beam || current_score == prev_score) && !heap.empty()) {
                            auto &top = heap.front();
                            prev_score = current_score;
                            current_score = top.first;
                            int index_P = top.second.first;
                            int index_M = top.second.second;

                            int key = valid_Ps[index_P];
                            State &state = beamP[key].alnobj;
                            int i1 = state.i1;
                            int i2 = state.i2;
                            int j1 = state.j1;

                            int k = i1 + i2 - 2;
                            State &mstate = bestM[k][sorted_bestM[k][index_M].second];
                            int newi1 = mstate.i1;
                            int newi2 = mstate.i2;
                            float transit_prob = hmmalign.trans_probs[mstate.endHMMstate-1][MANNER_ALN-1];
                            float alignscore = xlog_mul(xlog_mul(mstate.alignscore, transit_prob), state.alignscore);
                            float newscore = get<0>(M1_scores[index_P]) + mstate.score + weight * transit_prob;

                            cout << "newscore: " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << newscore << " " << index_P << " " <<  index_M << endl;
                            
                            pop_heap(heap.begin(), heap.end());
                            heap.pop_back();
                            int newkey = get_keys(j1, newi1, newi2);
                            if (beamM2[newkey].manner == MANNER_NONE) {
                                ++filled;
                                // update_if_better(beamM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                // ++nos_M2;
                                update_if_better(newi1, j1, newi2, s-j1, beamM2[newkey], // bestM2[get_keys(newi1, j1)][mnewi2][newj2], // beamM2[j2][make_pair(newi1, newi2)], // Note: not i1 i2 but newi1 newi2new
                                                get<1>(M1_scores[index_P]) + mstate.seq1foldscore,
                                                get<2>(M1_scores[index_P]) + mstate.seq2foldscore,
                                                mstate.manner, MANNER_M2_eq_M_plus_P, i1-1, i2-1,
                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                            } else {
                                if (beamM2[newkey].score < newscore - 1e-8)
                                    cout << "old: " << beamM2[newkey].score << " new: " << newscore << " " << newscore - 1e-8 << endl;
                                // assert(beamM2[newkey].score > newscore - 1e-8); // why? 
                            }
                        
                            ++index_M;
                            while (index_M < sorted_bestM[k].size()) {
                                mstate = bestM[k][sorted_bestM[k][index_M].second]; // sorted_bestM includes outside score
                                transit_prob = hmmalign.trans_probs[mstate.endHMMstate-1][MANNER_ALN-1];
                                alignscore = xlog_mul(xlog_mul(mstate.alignscore, transit_prob), state.alignscore);
                                float candidate_score = get<0>(M1_scores[index_P]) + mstate.score + weight * transit_prob;
                                int candidate_newkey = get_keys(j1, mstate.i1, mstate.i2);
                                
                                // candidate_score is a heuristic score
                                // float candidate_score = M1_scores[index_P] + sorted_bestM[k][index_M].first;
                                // int candidate_newi = sorted_bestM[k][index_M].second;
                                if (beamM2.find(candidate_newkey) == beamM2.end()) {
                                    heap.push_back(make_pair(candidate_score, make_pair(index_P, index_M)));
                                    push_heap(heap.begin(), heap.end());
                                    break;
                                } else {
                                    // based on the property of cube pruning, the new score must be worse
                                    // than the state already inserted
                                    // so we keep iterate through the candidate list to find the next
                                    // candidate
                                    ++index_M;
                                    if (beamM2[candidate_newkey].score < candidate_score - 1e-8)
                                        cout << "old: " << beamM2[candidate_newkey].score << " new: " << candidate_score << " " << candidate_score - 1e-8 << endl;
                                    // assert(beamM2[candidate_newkey].score > candidate_score - 1e-8); // why? 
                                }
                            }
                        }
                    }              
#endif

        // beam of M2
        // for every state in M2[j]
        //   1. multi-loop  (by extending M2 on the left)
        //   2. M = M2
        if (verbose) cout << "beam of M2" << endl;
        {
            // sort_keys(beamM2, keys);
            // cout << "beam of M2: " << beamM2.size() << endl;
            float threhold;
            if (beam > 0 && beamM2.size() > beam) threhold = beam_prune(beamM2, s, seq1_out_M2, seq2_out_M2);
            // cout << "beam of M2 after prune threhold: " << threhold << " " << beamM2.size() << endl;
            for (auto &item : beamM2) {
                State state = item.second; 
                int j1 = state.j1;
                int i1 = state.i1;
                int i2 = state.i2;
                int j2 = s - j1;
                
                // assert (j1 > 0); // lisiz TODO: DEBUG
                // assert (i1 > 0); // lisiz TODO: DEBUG
                // assert (i2 > 0); // lisiz TODO: DEBUG
                
                // 1. multi-loop
                {
                    for (int p1 = i1-1; p1 >= max(i1 - SINGLE_MAX_LEN - 1, 1); --p1) {
                        int nucp1 = seq1->nucs[p1];
                        int q1 = seq1->next_pair[nucp1][j1];
                        
                        // if (q1 != -1 && ((i1 - p1 - 1) <= SINGLE_MAX_LEN)) {
                        int i1_p1 = i1 - p1;
                        while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                            // single seq folding
                            if (seq1_out_Multi[q1].find(p1) == seq1_out_Multi[q1].end()) {
                                q1 = seq1->next_pair[nucp1][q1];
                                continue;
                            }

                            // the current shape is p..i M2 j ..q
                            int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                            int min_p2 = max(max(1, i2 - SINGLE_MAX_LEN - 1), hmmalign.low_bounds[p1]);
                            for (int p2 = max_p2; p2 >= min_p2; --p2) {
                                // hmm constraint
                                // if (p2 < hmmalign.low_bounds[p1] || p2 > hmmalign.up_bounds[p1]) continue;

                                int nucp2 = seq2->nucs[p2];
                                int q2 = seq2->next_pair[nucp2][j2];

                                // speed up
                                int i2_p2 = i2 - p2;
                                if (q2 > hmmalign.up_bounds[q1] || q2 == -1 || (i2_p2 + (q2 - j2) - 2 > SINGLE_MAX_LEN)) continue;
                                if (q2 < hmmalign.low_bounds[q1]) 
                                    q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1] - 1];

                                // if (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && ((i2 - p2 - 1) <= SINGLE_MAX_LEN)) {
                                while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                    // single seq folding
                                    if (seq2_out_Multi[q2].find(p2) == seq2_out_Multi[q2].end()) {
                                        q2 = seq2->next_pair[nucp2][q2];
                                        continue;
                                    }

                                    int newscore1 = multi_unpaired_score2(i1, j1, p1, q1, seq1);
                                    int newscore2 = multi_unpaired_score2(i2, j2, p2, q2, seq2);

                                    if (verbose) cout << "M2 to Multi: " << i1 << " " << j1 << " " << i2 << " " << j2  << " " << p1 << " " << q1 << " " << p2 << " " << q2 << endl;

                                    float pre_align_trans, post_align_trans, alignscore;
                                    // update bestMultiALN
                                    {
#ifdef dynalign
                                        aln_value_type alignscore = ALN_VALUE_MIN;
                                        if (p2>=hmmalign.low_bounds[p1] && p2<=hmmalign.up_bounds[p1]){ // only apply alignment constraint on left side
                                            alignscore = state.alignscore + abs(i1 - p1 - i2 + p2) + abs(q1 - j1 - q2 + j2);
                                        }
#else
                                        pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, 2);
                                        post_align_trans = get_hmm_score_right(j1, q1, j2, q2, 2, 2, true); // true

                                        alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                        alignscore = xlog_mul(alignscore, post_align_trans);

#endif 

                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].alnobj, // bestMultiALN[q1+q2][get_keys(q1, p1, p2)], 
                                                            state.seq1foldscore + newscore1,
                                                            state.seq2foldscore + newscore2, 
                                                            state.manner, MANNER_MULTI_ALN,
                                                            static_cast<char>(i1 - p1), q1 - j1,
                                                            static_cast<char>(i2 - p2), q2 - j2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                    }
                                    // update bestMultiINS1
                                    // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                    // if (p2+1-hmmalign.low_bounds[p1] >= 0)
#ifndef dynalign
                                    {
                                        pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, 2);
                                        post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, 2, 0, true); // true
                                        
                                        alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].ins1obj, // bestMultiINS1[q1+q2][get_keys(q1, p1, p2)], 
                                                            state.seq1foldscore + newscore1,
                                                            state.seq2foldscore + newscore2, 
                                                            state.manner, MANNER_MULTI_INS1,
                                                            static_cast<char>(i1 - p1), q1 - j1,
                                                            static_cast<char>(i2 - p2), q2 - j2,
                                                            alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                    }
                                    // update bestMultiINS2
                                    // align cost = (p, i) + (i, j) + (j, q); p1+1, q1-1
                                    // if (p2-hmmalign.low_bounds[p1+1] >= 0)
                                    {
                                        pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, 2);
                                        post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, 2, 1, true); // true
                                        
                                        alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                        alignscore = xlog_mul(alignscore, post_align_trans);

                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].ins2obj, // bestMultiINS2[q1+q2][get_keys(q1, p1, p2)],
                                                            state.seq1foldscore + newscore1,
                                                            state.seq2foldscore + newscore2, 
                                                            state.manner, MANNER_MULTI_INS2,
                                                            static_cast<char>(i1 - p1), q1 - j1,
                                                            static_cast<char>(i2 - p2), q2 - j2,
                                                            alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                    }
#endif
                                    q2 = seq2->next_pair[nucp2][q2];
                                } // while loop enumerate q2
                            } // for loop enumerate p2
                            q1 = seq1->next_pair[nucp1][q1];
                        } // q1
                    } // p1
                }// 1. multi-loop

                // 2. M = M2
                // single seq folding
                if (seq1_out_M[j1].find(i1) != seq1_out_M[j1].end() && seq2_out_M[j2].find(i2) != seq2_out_M[j2].end()) 
                {
                    if (verbose) cout << "M2 to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                    update_if_better(i1, j1, i2, j2, beamM[item.first], // bestM[key1][newi2][newj2], //beamM[j2][make_pair(i1, i2)], 
                                     state.seq1foldscore,
                                     state.seq2foldscore,
                                     state.manner, MANNER_M_eq_M2,
                                     state.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                }// 2. M = M2
            }// end beamM2
        }// beam of M2
     
        // beam of M
        // for every state in M[j]
        //   1. M = M + unpaired
        if (verbose) cout << "beam of M" << endl;
        {
            // cout << "beam of M: "  << beamM.size() << endl;
            float threshold = VALUE_FMIN;
            if (beam > 0 && beamM.size() > beam) threshold = beam_prune(beamM, s, seq1_out_M, seq2_out_M);

#ifdef is_cube_pruning
            sortM(threshold, beamM, sorted_bestM[s], s, seq1_out_M2, seq2_out_M2); // sorted_bestM include outside score
#endif

            for (auto &item : beamM) {
                State state = item.second; 
                int j1 = state.j1;
                int i1 = state.i1;
                int i2 = state.i2;
                int j2 = s - j1;
                
                // assert (j1 > 0); // lisiz TODO: DEBUG
                // assert (i1 > 0); // lisiz TODO: DEBUG
                // assert (i2 > 0); // lisiz TODO: DEBUG

                float trans_emit_prob; //  alignscore;
                if (j1 < seq1->seq_len - 1 && j2 < seq2->seq_len - 1){
                    // hmm constraint
                    if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]){
                        // single seq folding
                        if (seq1_out_M[j1+1].find(i1) != seq1_out_M[j1+1].end() && seq2_out_M[j2+1].find(i2) != seq2_out_M[j2+1].end()) {
#ifdef dynalign
                            int alignscore = state.alignscore; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_ALN, j1+1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(i1, j1+1, i2, j2+1, bestM[s+2][get_keys(j1+1, i1, i2)], // bestM[key1p1][newi2][newj2p1], // [s+2][j2+1][make_pair(i1, i2)],
                                                state.seq1foldscore,
                                                state.seq2foldscore,
                                                state.manner, MANNER_M_eq_M_plus_U_ALN,
                                                alignscore, state.startHMMstate, MANNER_ALN, weight, verbose);
                        }
                    }
                }
                if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) {
                    // hmm constraint
                    if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) {
                        // single seq folding
                        if (seq1_out_M[j1+1].find(i1) != seq1_out_M[j1+1].end()) {

#ifdef dynalign
                            int alignscore = state.alignscore + 1; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_INS1, j1+1, j2, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
    #endif 
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(i1, j1+1, i2, j2, bestM[s+1][get_keys(j1+1, i1, i2)], // bestM[key1p1][newi2][newj2p1], // [s+1][j2][make_pair(i1, i2)],
                                                state.seq1foldscore,
                                                state.seq2foldscore,
                                                state.manner, MANNER_M_eq_M_plus_U_INS1,
                                                alignscore, state.startHMMstate, MANNER_INS1, weight, verbose);
                        }
                    }
                }
                if (j1 <= seq1->seq_len - 1 && j2 < seq2->seq_len - 1) {
                    // hmm constraint
                    if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) {
                        // single seq folding
                        if (seq2_out_M[j2+1].find(i2) != seq2_out_M[j2+1].end()) {

#ifdef dynalign
                            int alignscore = state.alignscore + 1; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_INS2, j1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(i1, j1, i2, j2+1, bestM[s+1][get_keys(j1, i1, i2)], // bestM[key1][newi2][newj2p1], //s+1][j2+1][make_pair(i1, i2)],
                                                state.seq1foldscore,
                                                state.seq2foldscore,
                                                state.manner, MANNER_M_eq_M_plus_U_INS2,
                                                alignscore, state.startHMMstate, MANNER_INS2, weight, verbose);
                        }
                    }
                }

            } // for loop j1
        }// beam of M

        // beam of C
        // C = C + U
        if (verbose) cout << "beam of C" << endl;
        {    
            for (auto &item : beamC) {
                int j1 = item.first;
                int j2 = s - j1;
                if (j2 < 0) continue;
                State3 state3 = item.second;

                int newscore1 = external_unpaired_score(j1 + 1, seq1);
                int newscore2 = external_unpaired_score(j2 + 1, seq2);

                for (int m=0; m<3; m++){
                    State state;
                    HMMManner pre_manner;
                    switch (m)
                    {
                        case 0:
                            state = state3.ins1obj;
                            pre_manner = MANNER_INS1;
                            break;
                        case 1:
                            state = state3.ins2obj;
                            pre_manner = MANNER_INS2;
                            break;
                        case 2:
                            state = state3.alnobj;
                            pre_manner = MANNER_ALN;
                            break;
                        
                        default:
                            break;
                    }
                    if (state.endHMMstate == HMMMANNER_NONE) continue;
                    assert (state.endHMMstate == pre_manner);

                    if (verbose) cout << "C+U: " << m << " " << j1 << " "  << j2  << " "<<  j1 << " "  << j2  << " " << state.endHMMstate << endl;

                    float trans_emit_prob;
                    if (j1 == seq1->seq_len - 1 && j2 == seq2->seq_len - 1) {
                        // cout << "C+U: " << m << " " << j1 << " "  << j2  << " "<<  j1 << " "  << j2  << " " << state.score << endl;

                        trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_ALN, j1+1, j2+1, true);
                        float alignscore = xlog_mul(state.alignscore, trans_emit_prob);

                        if (alignscore > LOG_OF_ZERO)
                            update_if_better(0, j1+1, 0, j2+1, bestC[s+2][j1+1].alnobj,
                                            state.seq1foldscore,
                                            state.seq2foldscore,
                                            state.manner, MANNER_C_eq_C_plus_U_ALN,
                                            alignscore, state.endHMMstate, MANNER_ALN, weight, verbose);
                        
                        continue;
                    }

                    if (j1 < seq1->seq_len - 1 && j2 < seq2->seq_len - 1){ // ALN
                        // hmm constraint
                        if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]) {

#ifdef dynalign
                            int alignscore = state.alignscore; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_ALN, j1+1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 

                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1+1, 0, j2+1, bestC[s+2][j1+1].alnobj,
                                                state.seq1foldscore + newscore1,
                                                state.seq2foldscore + newscore2,
                                                state.manner, MANNER_C_eq_C_plus_U_ALN,
                                                alignscore, state.endHMMstate, MANNER_ALN, weight, verbose);
                        }
                    } 
                    if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) { // INS1
                        // hmm constraint
                        if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) {

#ifdef dynalign
                            int alignscore = state.alignscore + 1; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_INS1, j1+1, j2, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1+1, 0, j2, bestC[s+1][j1+1].ins1obj,
                                                state.seq1foldscore + newscore1,
                                                state.seq2foldscore,
                                                state.manner, MANNER_C_eq_C_plus_U_INS1,
                                                alignscore, state.endHMMstate, MANNER_INS1, weight, verbose);
                        }
                    } 
                    if (j2 < seq2->seq_len - 1 && j1 <= seq1->seq_len - 1) { // INS2
                        // hmm constraint
                        if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) {

#ifdef dynalign
                            int alignscore = state.alignscore + 1; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_INS2, j1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 

                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1, 0, j2+1, bestC[s+1][j1].ins2obj,
                                                state.seq1foldscore,
                                                state.seq2foldscore + newscore2,
                                                state.manner, MANNER_C_eq_C_plus_U_INS2,
                                                alignscore, state.endHMMstate, MANNER_INS2, weight, verbose);
                        }
                    }
                }
            } // for loop j1
        } // beam of C

    } // for loop s
    cout << "end of loop s" << endl;

    gettimeofday(&parse_endtime, NULL);
    parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    printf("seqs %d %d only parse time: %f seconds.\n", seq1_len, seq2_len, parse_elapsed_time);  

    processMem_t mem = GetProcessMemory();
    cout << "VmPeak: " << mem.VmPeak / 1024.0 / 1024.0 << endl;

    auto &state  = bestC[seq1->seq_len + seq2->seq_len][seq1->seq_len].alnobj; // bestC[seq1->seq_len - 1][seq2->seq_len-1];
    float bestscore = state.score;
    cout << "inside: " << state.score << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << xlog_mul(state.alignscore, hmmalign.trans_probs[state.endHMMstate-1][2]) << endl;
    // cout << "inside: " << state.manner << " " << state.startHMMstate << " " << state.endHMMstate << endl;

    // backtrace
    tuple<string, string, string, string> ret = get_parentheses(*seq1, *seq2, hmmalign);
    cout << get<0>(ret) <<endl;
    cout << get<1>(ret) <<endl;
    cout << get<2>(ret) <<endl;
    cout << get<3>(ret) <<endl;

#ifndef multilign
    return;
#else
    // outside
    cout << "start outside computation" << endl;
    // verbose = false;
    outside(limited, allowed_pairs);
    state  = bestC_beta[0][0]; // bestC[seq1->seq_len - 1][seq2->seq_len-1];
    cout << "outside: " << state.score << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << endl;

    // list of base pairs
    // float best_score = bestC[seq1->seq_len - 1 + seq2->seq_len - 1][seq1->seq_len-1].score;
    float min_score = bestscore * 1.2;
    cout << "min_score: " << min_score << endl;
    vector<tuple<float, int, int>> out_pairs_candidates;
    for(int s = 1; s < sum_len - 1; ++s) {
        unordered_map<int, State3>& beamP = bestP[s];
        unordered_map<int, State3>& beamP_beta = bestP_beta[s];
        for (auto &item : beamP) {
            // for (int m=0; m<3; m++){
            //     State *stateptr, *state_betaptr;
            //     switch (m)
            //     {
            //         case 0:
            //             stateptr = &item.second.ins1obj;
            //             state_betaptr = &beamP_beta[item.first].ins1obj;
            //             break;
            //         case 1:
            //             stateptr = &item.second.ins2obj;
            //             state_betaptr = &beamP_beta[item.first].ins2obj;
            //             break;
            //         case 2:
            //             stateptr = &item.second.alnobj;
            //             state_betaptr = &beamP_beta[item.first].alnobj;
            //             break;
                    
            //         default:
            //             break;
            //     }

            State& state = item.second.alnobj;
            State& state_beta = beamP_beta[item.first].alnobj;
            if (state.i1 <= 0) continue;

            if (state.score == VALUE_FMIN || state_beta.score == VALUE_FMIN) continue;

            float score = state.score + state_beta.score;

            if (min_score > 0) {
                if (score > min_score) continue;
            } else {
                if (score < min_score) continue;
            }

            int i1 = state.i1;
            int j1 = state.j1;
            // cout << i1 << " " << j1 << endl;
            if (limited) {
                if (allowed_pairs.find(make_pair(i1, j1)) == allowed_pairs.end()) 
                    continue;
            }

            out_pairs_candidates.push_back(make_tuple(score, state.i1, state.j1));

            // cout << state.i1 << " " << state.j1 << " " << state.i2 << " " << s - state.j1 <<  " ";
            // cout << state.score  << " " << state_beta.score << " " << state.score + state_beta.score << endl;
            // cout << state.seq1foldscore + state_beta.seq1foldscore << " " <<  state.seq2foldscore + state_beta.seq2foldscore << " ";
            // cout << state.seq1foldscore + state_beta.seq1foldscore + state.seq2foldscore + state_beta.seq2foldscore << " ";
            // cout << state.alignscore  << " " << state_beta.alignscore << " " << state.alignscore + state_beta.alignscore << endl;
        }
    }
    cout << "out_pairs_candidates size: " << out_pairs_candidates.size() << " " << num_pairs << endl;

    // only save top num_pairs base pairs using quickselect
    int current_num = out_pairs_candidates.size();
    if (current_num > num_pairs) {
        float threshold = quickselect(out_pairs_candidates, 0, current_num-1, current_num - num_pairs);
        float threshold2 = bestscore * 1.01;
        cout << "threshold: " << threshold  << " threshold2: " << threshold2 << endl;
        threshold = min(threshold, threshold2);
        cout << "threshold: " << threshold  << endl;
        for (auto &p : out_pairs_candidates) {
            if (get<0>(p) > threshold)
                out_pairs.push_back(make_pair(get<1>(p), get<2>(p)));
        }
    } else {
        for (auto &p : out_pairs_candidates) {
            out_pairs.push_back(make_pair(get<1>(p), get<2>(p)));
        }
    }
    cout << "out_pairs size: " << out_pairs.size() << " " << num_pairs << endl;

#endif

    // clear allocated space in alignment
    // cout << "clear allocated space in alignment" << endl;
#ifdef dynalign
    hmmalign.clear(false);
#else
    hmmalign.clear(true);
#endif
}

BeamSankoffParser::BeamSankoffParser(float aln_weight, int beam_size, int LFbeam, int LAbeam, bool if_aster, float energy_diff, bool is_verbose)
    :weight(aln_weight),
     beam(beam_size),
     lfbeam(LFbeam),
     alnbeam(LAbeam),
     use_astar(if_aster),
     max_energy_diff(energy_diff),
     verbose(is_verbose){

    cout << "beam : " << beam << " lfbeam: " << lfbeam << " alnbeam: " << alnbeam << " use_astar: " << use_astar << " max_energy_diff: " << max_energy_diff << endl;  
    
    initialize();

}


// int main(int argc, char** argv){
//     vector<string> seqs;
//     int beam_size, beam_size2;
//     float aln_weight;
//     bool is_verbose=false;

//     if (argc > 1) {
//         beam_size = atoi(argv[1]);
//         beam_size2 = atoi(argv[2]);
//         aln_weight = atof(argv[3]);
//         is_verbose = atoi(argv[4]) == 1;
//     }

//     int seq1len, seq2len;
//     for (string seq; getline(cin, seq);) {
//         if (seq.length() == 0)
//             continue;

//         if (seq[0] == ';' || seq[0] == '>') {
//             printf("%s\n", seq.c_str());
//             continue;
//         }

//         // convert to uppercase
//         transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

//         // convert T to U
//         replace(seq.begin(), seq.end(), 'T', 'U');

//         printf("%s\n", seq.c_str());

//         seqs.push_back(seq);

//         if (seqs.size() == 1) {
//             seq1len = seq.size();
//         }  
//         else if (seqs.size() == 2) {
//             seq2len = seq.size();

//             struct timeval parse_starttime, parse_endtime;
//             gettimeofday(&parse_starttime, NULL);
            
//             if (seq1len < seq2len){ // TODO
//                 BeamSankoffParser parser(beam_size, beam_size2, seq2len, seq1len, aln_weight, is_verbose);
//                 swap(seqs[0], seqs[1]);
//                 parser.parse(seqs);
//             }
//             else {
//                 BeamSankoffParser parser(beam_size, beam_size2, seq1len, seq2len, aln_weight, is_verbose);
//                 parser.parse(seqs);
//             }
                
//             gettimeofday(&parse_endtime, NULL);
//             double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
//             printf("seqs %d %d time: %f seconds.\n", seq1len, seq2len, parse_elapsed_time);
//         }
//     }
// }