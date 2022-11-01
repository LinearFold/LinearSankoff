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


// #ifdef is_cube_pruning
// // cube pruning
// bool comparefunc(pair<int,State> a, pair<int,State> b) {
//     return a.first > b.first;
// }
// void BeamSankoffParser::sort_keys(unordered_map<int, State3> &map, vector<pair<int,State> > &sorted_keys) {
//     sorted_keys.clear();
//     for(auto &item : map) {
//         sorted_keys.push_back(make_pair(item.first, item.second.alnobj));
//     }
//     sort(sorted_keys.begin(), sorted_keys.end(), comparefunc);    
// }
// void BeamSankoffParser::sortM(float threshold,
//                               unordered_map<int, State> &beamstep,
//                               vector<pair<float, int>> &sorted_stepM,
//                               int s,
//                               vector<unordered_map<int, int> > seq1_outside, 
//                               vector<unordered_map<int, int> > seq2_outside) {
//     sorted_stepM.clear();
//     if (threshold == VALUE_FMIN) {
//         // no beam pruning before, so scores vector not usable
//         for (auto &item : beamstep) {
//             State &cand = item.second;
//             int i1 = cand.i1;
//             int i2 = cand.i2;
//             int j1 = cand.j1;
//             int j2 = s - j1;

//             int k1 = i1 - 1;
//             int k2 = i2 - 1;

//             // get new score
//             float newscore = VALUE_MIN;
//             // alignment score
//             float alignscore = xlog_mul(xlog_mul(bestC[k1+k2][k1].alignscore, cand.alignscore), hmmalign.trans_probs[bestC[k1+k2][k1].endHMMstate-1][2]);
//             float backward_score = hmmalign.bestALN[s][j2].beta;
//             alignscore = xlog_mul(alignscore, backward_score);
//             // alignscore = xlog_div(alignscore, aln_viterbi);

//             // folding score
//             int foldingscore = cand.seq1foldscore + cand.seq2foldscore;
//             // folding heuristic
//             int seq1_out = VALUE_MIN;
//             int seq2_out = VALUE_MIN;
//             if ((seq1_outside[j1].find(i1) != seq1_outside[j1].end())) {
//                 seq1_out = seq1_outside[j1][i1];
//             }
//             if ((seq2_outside[j2].find(i2) != seq2_outside[j2].end())) {
//                 seq2_out = seq2_outside[j2][i2];
//             }

//             if (seq1_out == VALUE_MIN || seq2_out == VALUE_MIN) {
//                 newscore = VALUE_FMIN;
//             } else {
//                 foldingscore += seq1_out + seq2_out;
//                 newscore = foldingscore + weight * alignscore;
//             }

//             // float newscore = (k >= 0 ? bestC[k][i1-1].score : 0) + cand.score; // old
//             // cout << "sorting M " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << newscore << endl;
//             if (newscore > VALUE_FMIN)
//                 sorted_stepM.push_back(make_pair(newscore, item.first));
//         }
//     } else {
//         for (auto &p : scores) {
//             if (p.first >= threshold) sorted_stepM.push_back(p);
//         }
//     }
//     sort(sorted_stepM.begin(), sorted_stepM.end(), greater<pair<float, int>>());
// }
// #endif

tuple<string, string, string, string> BeamSankoffParser::get_parentheses_H(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
    if ( !stk.empty() ) {
        tuple<int, int, int, int, State> top = stk.top();
        int i1 = get<0>(top), j1 = get<1>(top);
        int i2 = get<2>(top), j2 = get<3>(top);
        State& state = get<4>(top);
        stk.pop();

        if (verbose) printf("get_parentheses_H manner at %d, %d, %d, %d: manner %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.manner, state.alignscore, state.seq1foldscore, state.seq2foldscore);

        switch (state.endHMMstate) {
            case MANNER_ALN:
                {
                    string sub_struc1 = struc_hairpin_pair(j1 - i1 -1);
                    string sub_struc2 = struc_hairpin_pair(j2 - i2 -1);
                    pair<string, string> aln_ret = get_hmm_aln(i1, j1, i2, j2, MANNER_ALN, MANNER_ALN);

                    return make_tuple(sub_struc1, sub_struc2, get<0>(aln_ret), get<1>(aln_ret));
                }
                break;
            case MANNER_INS1:
                {
                    string sub_struc1 = struc_hairpin_pair(j1 - i1 -1);
                    string sub_struc2 = struc_hairpin_dots(j2 - i2 -1);
                    pair<string, string> aln_ret = get_hmm_aln(i1, j1, i2+1, j2-1, MANNER_INS1, MANNER_INS1);
                    return make_tuple(sub_struc1, sub_struc2, get<0>(aln_ret), get<1>(aln_ret));
                }
                break;
            case MANNER_INS2:
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

        // stk.push(make_tuple(p1, q1, p2, q2, bestM2[q1+q2][get_keys(q1, p1, p2)]));
        stk.push(make_tuple(p1, q1, p2, q2, bestM2[q1+q2][q1][make_pair(p1, p2)]));
        tuple<string, string, string, string> m2result = get_parentheses_M2(seq1, seq2, stk, hmmalign);

        pair<string, string> pre_aln, post_aln;
        string seq1_struc, seq2_struc, seq1_aln, seq2_aln;
        switch(state.endHMMstate) {
            case MANNER_ALN:
                pre_aln = get_hmm_aln_left(i1, p1, i2, p2, MANNER_ALN, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1, q2, j2, MANNER_ALN, MANNER_ALN);

                seq1_struc = struc_p2p_pair(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_pair(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            case MANNER_INS1:
                pre_aln = get_hmm_aln_left(i1, p1, i2+1, p2, MANNER_INS1, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1, q2, j2-1, MANNER_ALN, MANNER_INS1);

                seq1_struc = struc_p2p_pair(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_dots(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            case MANNER_INS2:
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

        if (verbose) {
            cout << seq1_aln << endl;
            cout << seq1_struc << endl;
            cout << seq2_aln << endl;
            cout << seq2_struc << endl;
        } 
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

        if (verbose) printf("get_parentheses_M2 manner at %d, %d, %d, %d: manner %d %d, alignscore: %f, foldscore: %d %d\n", i1, j1, i2, j2, state.premanner, state.manner, state.alignscore, state.seq1foldscore, state.seq2foldscore);

        int k1 = state.trace1.split;
        int k2 = state.trace2.split; // N.B.

        switch (state.premanner) {
            case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_P:
                stk.push(make_tuple(i1, k1, i2, k2, bestM[k1+k2][k1][2][make_pair(i1, i2)])); // push M1
                break;
            case MANNER_M_eq_M_plus_U_INS1: case MANNER_M_eq_P1:
                stk.push(make_tuple(i1, k1, i2, k2, bestM[k1+k2][k1][0][make_pair(i1, i2)])); // push M1
                break;
            case MANNER_M_eq_M_plus_U_INS2: case MANNER_M_eq_P2:
                stk.push(make_tuple(i1, k1, i2, k2, bestM[k1+k2][k1][1][make_pair(i1, i2)])); // push M1
                break;
            case MANNER_M_eq_M2: {
                HMMManner endHMMstate = bestM2[k1+k2][k1][make_pair(i1, i2)].endHMMstate;
                stk.push(make_tuple(i1, k1, i2, k2, bestM[k1+k2][k1][endHMMstate-1][make_pair(i1, i2)])); // push M1
            }
            break;
            default:
                printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                assert(false);
        }
        // stk.push(make_tuple(i1, k1, i2, k2, bestM[k1+k2][k1][make_pair(i1, i2)])); // push M1
        tuple<string, string, string, string> pre_result = get_parentheses_M1(seq1, seq2, stk, hmmalign);

        switch (state.manner) {
            case MANNER_M2_eq_M_plus_P1:
            {
                cout << "MANNER_M2_eq_M_plus_P1: " << i1 << " " << i2 << " " << k1+1 << " " << j1 << " " << k2+1 << " " << j2 << endl;

                tuple<string, string, string, string> branch_insertion = backtrace_branch_insertion(k1+1, j1, k2+1, j2);

                return make_tuple(get<0>(pre_result)+get<0>(branch_insertion), get<1>(pre_result), get<2>(pre_result)+get<2>(branch_insertion), get<3>(pre_result)+get<3>(branch_insertion));
                // return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                break;
            }
            case MANNER_M2_eq_M_plus_P2:
            {
                cout << "MANNER_M2_eq_M_plus_P2: " << i1 << " " << i2 << " " << k1+1 << " " << j1 << " " << k2+1 << " " << j2 << endl;

                tuple<string, string, string, string> branch_insertion = backtrace_branch_insertion(k1+1, j1, k2+1, j2);
                
                return make_tuple(get<0>(pre_result), get<1>(pre_result)+get<1>(branch_insertion), get<2>(pre_result)+get<2>(branch_insertion), get<3>(pre_result)+get<3>(branch_insertion));
                // return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                break;
            }
            default:
            {
                stack<tuple<int, int, int, int, State>> stk2; // push P
                // stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][get_keys(j1, k1+1, k2+1)].alnobj)); 
                stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][j1][2][make_pair(k1+1, k2+1)])); 
                tuple<string, string, string, string> result = get_parentheses_P(seq1, seq2, stk2, hmmalign);

                return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
            }
        }
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
            {
                if (verbose) cout << "MANNER_M_eq_M_plus_U_ALN" << endl;

                switch (state.premanner) {
                    case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_P:  
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[j1-1+j2-1][j1-1][2][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M_plus_U_INS1: case MANNER_M_eq_P1:
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[j1-1+j2-1][j1-1][0][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M_plus_U_INS2: case MANNER_M_eq_P2:
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[j1-1+j2-1][j1-1][1][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M2:
                        HMMManner endHMMstate = bestM2[j1-1+j2-1][j1-1][make_pair(i1, i2)].endHMMstate;
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[j1-1+j2-1][j1-1][endHMMstate-1][make_pair(i1, i2)]));
                        break;
                }   

                // stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[j1-1+j2-1][j1-1][make_pair(i1, i2)]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result)+".", get<1>(result)+".", get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+seq2.raw_seq.at(j2));
                break;
            }
            case MANNER_M_eq_M_plus_U_INS1:
            {
                if (verbose) cout << "MANNER_M_eq_M_plus_U_INS1" << endl;

                switch (state.premanner) {
                    case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_P:
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[j1-1+j2][j1-1][2][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M_plus_U_INS1: case MANNER_M_eq_P1:
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[j1-1+j2][j1-1][0][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M_plus_U_INS2: case MANNER_M_eq_P2:
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[j1-1+j2][j1-1][1][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M2:
                        HMMManner endHMMstate = bestM2[j1-1+j2][j1-1][make_pair(i1, i2)].endHMMstate;
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[j1-1+j2][j1-1][endHMMstate-1][make_pair(i1, i2)]));
                        break;
                }

                // stk.push(make_tuple(i1, j1-1, i2, j2, bestM[j1-1+j2][j1-1][make_pair(i1, i2)]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result)+".", get<1>(result), get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+"-");
                break;
            }
            case MANNER_M_eq_M_plus_U_INS2:
            {
                if (verbose) cout << "MANNER_M_eq_M_plus_U_INS2" << endl;
                
                switch (state.premanner) {
                    case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_P:
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[j1+j2-1][j1][2][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M_plus_U_INS1: case MANNER_M_eq_P1:
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[j1+j2-1][j1][0][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M_plus_U_INS2: case MANNER_M_eq_P2:
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[j1+j2-1][j1][1][make_pair(i1, i2)]));
                        break;
                    case MANNER_M_eq_M2:
                        HMMManner endHMMstate = bestM2[j1+j2-1][j1][make_pair(i1, i2)].endHMMstate;
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[j1+j2-1][j1][endHMMstate-1][make_pair(i1, i2)]));
                        break;
                }

                // stk.push(make_tuple(i1, j1, i2, j2-1, bestM[j1+j2-1][j1][make_pair(i1, i2)]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result), get<1>(result)+".", get<2>(result)+"-", get<3>(result)+seq2.raw_seq.at(j2));
                break;
            }
            case MANNER_M_eq_M2:
            {
                if (verbose) cout << "MANNER_M_eq_M2" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestM2[j1+j2][get_keys(j1, i1, i2)]));
                stk.push(make_tuple(i1, j1, i2, j2, bestM2[j1+j2][j1][make_pair(i1, i2)]));
                return get_parentheses_M2(seq1, seq2, stk, hmmalign);
                break;
            }
            case MANNER_M_eq_P:
            {
                if (verbose) cout << "MANNER_M_eq_P" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestP[j1+j2][get_keys(j1, i1, i2)].alnobj));
                stk.push(make_tuple(i1, j1, i2, j2, bestP[j1+j2][j1][2][make_pair(i1, i2)]));
                return get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;  
            }
            case MANNER_M_eq_P1: // case MANNER_M_eq_P2:
            {
                // if (verbose) 
                cout << "MANNER_M_eq_P1 " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                
                return backtrace_branch_insertion(i1, j1, i2, j2);
                break;  
            }
            case MANNER_M_eq_P2:
            {
                // if (verbose) 
                cout << "MANNER_M_eq_P2 " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                
                return backtrace_branch_insertion(i1, j1, i2, j2);
                break;  
            }
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
        
        switch(state.manner) {
            case MANNER_HAIRPIN_ALN: case MANNER_HAIRPIN_INS1: case MANNER_HAIRPIN_INS2:
                stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][j1][state.endHMMstate-1][make_pair(i1, i2)]));
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
            case MANNER_MULTI_ALN: case MANNER_MULTI_INS1: case MANNER_MULTI_INS2:
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][j1][state.endHMMstate-1][make_pair(i1, i2)]));
                return get_parentheses_Multi(seq1, seq2, stk, hmmalign);
        }

        // inner 
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
        HMMManner preHMMstate;
        switch(state.premanner)
        {
            case MANNER_SINGLE_ALN:
                if (verbose) cout << "MANNER_SINGLE_ALN" << endl;
                stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][q1][2][make_pair(p1, p2)]));
                preHMMstate = MANNER_ALN; // bestP[q1+q2][q1][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_SINGLE_INS1:
                if (verbose) cout << "MANNER_SINGLE_INS1" << endl; 
                stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][q1][0][make_pair(p1, p2)]));
                preHMMstate = MANNER_INS1; // bestP[q1+q2][q1][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_SINGLE_INS2:
                if (verbose) cout << "MANNER_SINGLE_INS1" << endl;
                stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][q1][1][make_pair(p1, p2)]));
                preHMMstate = MANNER_INS2; // bestP[q1+q2][q1][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_ALN: 
                if (verbose) cout << "MANNER_HAIRPIN_ALN" << endl;
                stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][q1][2][make_pair(p1, p2)]));
                preHMMstate = MANNER_ALN; //  bestH[q1+q2][q1][2][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_HAIRPIN_INS1: 
                if (verbose) cout << "MANNER_HAIRPIN_INS1" << endl;
                stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][q1][0][make_pair(p1, p2)]));
                preHMMstate = MANNER_INS1; //  bestH[q1+q2][q1][2][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
             case MANNER_HAIRPIN_INS2: 
                if (verbose) cout << "MANNER_HAIRPIN_INS2" << endl;
                stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][q1][1][make_pair(p1, p2)]));
                preHMMstate = MANNER_INS2; //  bestH[q1+q2][q1][2][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_ALN: 
                if (verbose) cout << "MANNER_P_eq_MULTI_ALN" << endl; 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][q1][2][make_pair(p1, p2)]));
                preHMMstate = MANNER_ALN; // bestMulti[q1+q2][q1][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_INS1: 
                if (verbose) cout << "MANNER_P_eq_MULTI_INS1" << endl; 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][q1][0][make_pair(p1, p2)]));
                preHMMstate = MANNER_INS1; // bestMulti[q1+q2][q1][make_pair(p1, p2)].endHMMstate;
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_INS2: 
                if (verbose) cout << "MANNER_P_eq_MULTI_INS2" << endl;
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][q1][1][make_pair(p1, p2)]));
                preHMMstate = MANNER_INS2; // bestMulti[q1+q2][q1][make_pair(p1, p2)].endHMMstate;
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

            switch (preHMMstate)
            {
                case MANNER_ALN:
                    {
                        hmm_pre_end = MANNER_ALN;
                        hmm_post_start = MANNER_ALN;
                    }
                    break;
                case MANNER_INS1:
                    {
                        // aln_p2++;
                        // aln_q2--; 

                        hmm_pre_end = MANNER_INS1;
                        hmm_post_start = MANNER_INS1;
                        break;
                    }
                case MANNER_INS2:
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

            switch (state.endHMMstate) {
                case MANNER_ALN: 
                    {
                        // pre_aln = get_hmm_aln(i1, aln_p1-1, i2, aln_p2-1, MANNER_ALN, hmm_pre_end);
                        // post_aln = get_hmm_aln(aln_q1+1, j1, aln_q2+1, j2, hmm_post_start, MANNER_ALN);
                        
                        pre_aln = get_hmm_aln_left(i1, aln_p1, i2, aln_p2, MANNER_ALN, hmm_pre_end);
                        post_aln = get_hmm_aln_right(aln_q1, j1, aln_q2, j2, hmm_post_start, MANNER_ALN);
                        
                        seq1_struc = struc_p2p_pair(get<0>(result), seq1_l1-1, seq1_l2-1);
                        seq2_struc = struc_p2p_pair(get<1>(result), seq2_l1-1, seq2_l2-1);
                        break;
                    }
                case MANNER_INS1:
                    {
                        // pre_aln = get_hmm_aln(i1, aln_p1-1, i2+1, aln_p2-1, MANNER_INS1, hmm_pre_end);
                        // post_aln = get_hmm_aln(aln_q1+1, j1, aln_q2+1, j2-1, hmm_post_start, MANNER_INS1);

                        pre_aln = get_hmm_aln_left(i1, aln_p1, i2+1, aln_p2, MANNER_INS1, hmm_pre_end);
                        post_aln = get_hmm_aln_right(aln_q1, j1, aln_q2, j2-1, hmm_post_start, MANNER_INS1);

                        seq1_struc = struc_p2p_pair(get<0>(result), seq1_l1-1, seq1_l2-1);
                        seq2_struc = struc_p2p_dots(get<1>(result), seq2_l1-1, seq2_l2-1);
                        break;
                    }
                case MANNER_INS2:
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

        if (verbose) {
            cout << seq1_aln << endl;
            cout << seq1_struc << endl;
            cout << seq2_aln << endl;
            cout << seq2_struc << endl;
        }

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

                    Manner premanner = state.premanner;
                    switch (premanner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].alnobj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins1obj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1])); // push C
                    tuple<string, string, string, string> pre_result = get_parentheses_C(seq1, seq2, stk, hmmalign);

                    stack<tuple<int, int, int, int, State>> stk2; // push P
                    // stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][get_keys(j1, k1+1, k2+1)].alnobj)); 
                    stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][j1][2][make_pair(k1+1, k2+1)])); 
                    tuple<string, string, string, string> result = get_parentheses_P(seq1, seq2, stk2, hmmalign);

                    return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                    break;
                }

            case MANNER_C_eq_C_plus_P1:
                {   
                    int k1 = state.trace1.split;
                    int k2 = state.trace2.split;

                    // if (verbose) 
                    cout <<  "MANNER_C_eq_C_plus_P1 " << k1 << " " << k2 << endl;

                    Manner premanner = state.premanner;
                    switch (premanner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].alnobj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins1obj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1])); // push C
                    tuple<string, string, string, string> pre_result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    tuple<string, string, string, string> branch_insertion = backtrace_branch_insertion(k1+1, j1, j2, j2);

                    // return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                    return make_tuple(get<0>(pre_result)+get<0>(branch_insertion), get<1>(pre_result), get<2>(pre_result)+get<2>(branch_insertion), get<3>(pre_result)+get<3>(branch_insertion));
                    break;
                }
            
            case MANNER_C_eq_C_plus_P2:
                {   
                    int k1 = state.trace1.split;
                    int k2 = state.trace2.split;

                    // if (verbose) 
                    cout <<  "MANNER_C_eq_C_plus_P2 " << k1 << " " << k2 << endl;

                    Manner premanner = state.premanner;
                    switch (premanner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].alnobj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins1obj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1].ins2obj));
                            break;
                        
                        default:
                            break;
                    }
                    // stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1])); // push C
                    tuple<string, string, string, string> pre_result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    tuple<string, string, string, string> branch_insertion = backtrace_branch_insertion(j1, j1, k2+1, j2);

                    // return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                    return make_tuple(get<0>(pre_result), get<1>(pre_result)+get<1>(branch_insertion), get<2>(pre_result)+get<2>(branch_insertion), get<3>(pre_result)+get<3>(branch_insertion));
                    break;
                }
            case MANNER_C_eq_C_plus_U_ALN:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_ALN" << endl;
                    Manner premanner = state.premanner;
                    switch (premanner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1+j2-1][j1-1].alnobj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1+j2-1][j1-1].ins1obj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
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
                    Manner premanner = state.premanner;
                    switch (premanner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1+j2][j1-1].alnobj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1+j2][j1-1].ins1obj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
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
                    Manner premanner = state.premanner;
                    switch (premanner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1+j2-1][j1].alnobj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1+j2-1][j1].ins1obj));
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
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
    Manner premanner = bestC[seq1.seq_len + seq2.seq_len][seq1.seq_len].alnobj.premanner;
    switch (premanner)
    {
        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1.seq_len - 1 + seq2.seq_len-1][seq1.seq_len-1].alnobj));
            break;
        case MANNER_C_eq_C_plus_U_INS1:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1.seq_len - 1 + seq2.seq_len-1][seq1.seq_len-1].ins1obj));
            break;
        case MANNER_C_eq_C_plus_U_INS2:
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

    return hmmalign.all_local_scores[s1][i1][i2][j1][j2];
}

float BeamSankoffParser::get_hmm_score_left(int i1, int j1, int i2, int j2, int s1, int s2){
    // (-->(
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return xlog(0.0);
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return xlog(0.0);

    if ((j1-i1+1 > 30) || (j2-i2+1 > 30))
        return hmmalign.viterbi_path_local_left(i1, j1, i2, j2, static_cast<HMMManner>(s1+1), static_cast<HMMManner>(s2+1));

    return hmmalign.left_local_scores[s1][s2][i1][i2][j1][j2];
}

float BeamSankoffParser::get_hmm_score_right(int i1, int j1, int i2, int j2, int s1, int s2, bool allowout){
    // )-->)
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return xlog(0.0);
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return xlog(0.0);

    if ((j1-i1+1 > 30) || (j2-i2+1 > 30)) {
        return hmmalign.viterbi_path_local_right(i1, j1, i2, j2, static_cast<HMMManner>(s1+1), static_cast<HMMManner>(s2+1));
    }

    return hmmalign.right_local_scores[s1][s2][i1][i2][j1][j2];
}

string BeamSankoffParser::backtrace_single_seq(SeqObject* seq, int i, int j, vector<unordered_map<int, int> > &seq_in_P) {
    int nuci = seq->nucs[i];
    int nuci1 = seq->nucs[i+1];
    int nucj = seq->nucs[j];
    int nucj_1 = seq->nucs[j-1];

    // hairpin
    int tetra_hex_tri = -1;
    if (j-i-1 == 4) // 6:tetra
        tetra_hex_tri = seq->if_tetraloops[i];
    else if (j-i-1 == 6) // 8:hexa
        tetra_hex_tri = seq->if_hexaloops[i];
    else if (j-i-1 == 3) // 5:tri
        tetra_hex_tri = seq->if_triloops[i];

    int hairpin_score = -v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri); // hairpinScore(i, j, seq);
    if (hairpin_score == seq_in_P[j][i]) 
        return struc_hairpin_pair(j-i-1);

    int maxp=-1;
    int maxq=-1;
    string max_sub_struc;
    for (int p=i+1; p<=min(i+MAX_LOOP_LEN, j-1); p++) {
        int nucp = seq->nucs[p];
        int nucp_1 = seq->nucs[p-1]; 

        int i_p = i - p;
        int q = seq->next_pair[nucp][max(p+3, j-MAX_LOOP_LEN)];
        while (q != -1 && q < j && j - q + 1 > MAX_LOOP_LEN) {
            q = seq->next_pair[nucp][q];
        }

        while (q != -1 && q < j) {
            if (seq_in_P[q].find(p) != seq_in_P[q].end()) {
                int nucq = seq->nucs[q];
                int nucq1 = seq->nucs[q+1];
                
                // stacking, internal loop, bulges
                int newscore = seq_in_P[q][p] - v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                
                if (newscore == seq_in_P[j][i]) {
                    maxp = p;
                    maxq = q;
                    max_sub_struc = backtrace_single_seq(seq, p, q, seq_in_P);
                    break;
                }
            }

            q = seq->next_pair[nucp][q]; // next q
        }
    }
    
    assert (maxp > 0);

    return struc_p2p_pair(max_sub_struc, maxp-i-1, j-maxq-1);
}

tuple<string, string, string, string> BeamSankoffParser::backtrace_branch_insertion(int i1, int j1, int i2, int j2) {
    if (i2 == (j2+1)) { // branch insertion in seq1
        string seq1_struc = backtrace_single_seq(&sequences[0], i1, j1, seq1_in_P);
        pair<string, string> alignment = get_hmm_aln(i1, j1, i2, j2, MANNER_INS1, MANNER_INS1);
        cout << seq1_struc << " " << alignment.first << " " << alignment.second << endl;
        return make_tuple(seq1_struc, "", alignment.first, alignment.second);
    } 
    else if (i1 == (j1+1)) { // branch insertion in seq2
        string seq2_struc = backtrace_single_seq(&sequences[1], i2, j2, seq2_in_P);
        pair<string, string> alignment = get_hmm_aln(i1, j1, i2, j2, MANNER_INS2, MANNER_INS2);
        cout << seq2_struc << " " << alignment.first << " " << alignment.second << endl;
        return make_tuple("", seq2_struc, alignment.first, alignment.second);
    } else {
        cout << "wrong parameters for backtrace_branch_insertion: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
    }

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

// int BeamSankoffParser::get_keys(int j1, int i1, int i2){
//     return (j1 * seq1_len + i1) * seq2_len + i2;
// }

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

unsigned long quickselect_partition(vector<tuple<float, int, int, pair<int, int>> >& scores, unsigned long lower, unsigned long upper) {
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
float quickselect(vector<tuple<float, int, int, pair<int, int>> >& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return get<0>(scores[lower]);
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return get<0>(scores[split]);
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

float BeamSankoffParser::beam_prune(unordered_map<pair<int, int>, State, pair_hash> **beamstep, int s, vector<unordered_map<int, int> > seq1_outside, vector<unordered_map<int, int> > seq2_outside) {
    vector<tuple<float, int, int, pair<int, int>> > candidates;
    for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
        int j2 = s - j1;
        if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

        for (int m=0; m<3; m++) {
            for (auto &item: beamstep[j1][m]) {
                State &cand = item.second;

                int i1 = cand.i1;
                int i2 = cand.i2;

                int k1, k2;
                float forward_score, backward_score, alignscore, newscore; // alignment backward/heuristic score
                switch (m) 
                {
                    case 0:
                    {
                        k1 = i1 - 1;
                        k2 = i2 + 1;
                        backward_score = ins1_backward_score[j1][j2-1];
                        break;
                    }
                    case 1:
                    {
                        k1 = i1 + 1;
                        k2 = i2 - 1;
                        backward_score = ins2_backward_score[j1-1][j2];
                        break;
                    }
                    default:
                    {
                        k1 = i1 - 1;
                        k2 = i2 - 1;
                        backward_score = aln_backward_score[j1][j2];
                        break;
                    }
                }

                // alignment forward score
                forward_score = xlog_mul(bestC[k1+k2][k1].alnobj.alignscore, hmmalign.trans_probs[2][cand.startHMMstate-1]);
                float tmp_score = xlog_mul(bestC[k1+k2][k1].ins1obj.alignscore, hmmalign.trans_probs[0][cand.startHMMstate-1]);
                forward_score = (tmp_score > forward_score) ? tmp_score : forward_score;
                tmp_score = xlog_mul(bestC[k1+k2][k1].ins2obj.alignscore, hmmalign.trans_probs[1][cand.startHMMstate-1]);
                forward_score = (tmp_score > forward_score) ? tmp_score : forward_score;
                // joint alignment score
                alignscore = xlog_mul(xlog_mul(forward_score, cand.alignscore), backward_score);

                // + folding score
                // joint folding score
                int foldingscore = cand.seq1foldscore + cand.seq2foldscore;
                // branch insertion P1/P2 to M
                int seq1_out, seq2_out; 
                if (i1 == (j1+1)) seq1_out = seq1_mfe;
                else seq1_out = seq1_outside[j1][i1];
                if (i2 == (j2+1)) seq2_out = seq2_mfe;
                else seq2_out = seq2_outside[j2][i2];

                foldingscore += seq1_out + seq2_out;

                // newscore
                if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN)
                    newscore = LOG_OF_ZERO;
                else
                    newscore = foldingscore + weight * alignscore;

                candidates.push_back(make_tuple(newscore, j1, m, item.first));
            }
        }
    }

    if (candidates.size() <= beam) return VALUE_FMIN;
    
    float threshold = quickselect(candidates, 0, candidates.size() - 1, candidates.size() - beam);
    for (auto &p : candidates) {
        float score = get<0>(p);
        if (score < threshold) {
            int j1 = get<1>(p);
            int m = get<2>(p);
            pair<int, int> key = get<3>(p);
            beamstep[j1][m].erase(key);
        }
    }

    return threshold;
}

float BeamSankoffParser::beam_prune(unordered_map<pair<int, int>, State, pair_hash> *beamstep, int s, vector<unordered_map<int, int> > seq1_outside, vector<unordered_map<int, int> > seq2_outside) {
    int count = 0;
    for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
        int j2 = s - j1;
        if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;
        count += beamstep[j1].size();
    }
    if (count <= beam) return VALUE_FMIN;

    vector<tuple<float, int, int, pair<int, int>> > candidates;
    for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
        int j2 = s - j1;
        if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

        unordered_map<pair<int, int>, State, pair_hash>& items = beamstep[j1];
        for (auto &item: items) {
            State &cand = item.second;

            int i1 = cand.i1;
            int i2 = cand.i2;

            int k1 = i1 - 1;
            int k2 = i2 - 1;

            State& prefix_C = bestC[k1+k2][k1].alnobj; // todo
            // // HMMManner forward_endmanner = prefix_C.endHMMstate;
            // if (prefix_C.manner == MANNER_NONE) {
            //     candidates.push_back(make_tuple(LOG_OF_ZERO, j1, make_pair(i1, i2)));
            //     continue;
            // }

            float newscore; // final score
            if (use_astar) {
                // + alignment sore 
                float alignscore, forward_score, backward_score;
#ifdef dynalign
                int alignscore = prefix_C.alignscore + cand.alignscore + abs(seq2_len-j2-seq1_len+j1);
#else       
                // alignment forward score
                forward_score = xlog_mul(bestC[k1+k2][k1].alnobj.alignscore, hmmalign.trans_probs[2][cand.startHMMstate-1]);
                float tmp_score = xlog_mul(bestC[k1+k2][k1].ins1obj.alignscore, hmmalign.trans_probs[0][cand.startHMMstate-1]);
                forward_score = (tmp_score > forward_score) ? tmp_score : forward_score;
                tmp_score = xlog_mul(bestC[k1+k2][k1].ins2obj.alignscore, hmmalign.trans_probs[1][cand.startHMMstate-1]);
                forward_score = (tmp_score > forward_score) ? tmp_score : forward_score;
            
                // alignment backward/heuristic score
                switch (cand.endHMMstate) {
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

                // joint alignment score
                alignscore = xlog_mul(xlog_mul(forward_score, cand.alignscore), backward_score);
#endif
                // + folding score
                // joint folding score
                int foldingscore = cand.seq1foldscore + cand.seq2foldscore;
                // folding heuristic
                if (!use_suffix) {
                    // foldingscore += seq1_outside[j1][i1] + seq2_outside[j2][i2];

                    // branch insertion P1/P2 to M
                    int seq1_out, seq2_out; 
                    if (i1 == (j1+1)) seq1_out = seq1_mfe;
                    else seq1_out = seq1_outside[j1][i1];
                    if (i2 == (j2+1)) seq2_out = seq2_mfe;
                    else seq2_out = seq2_outside[j2][i2];

                    foldingscore += seq1_out + seq2_out;
                } else {
                    foldingscore += prefix_C.seq1foldscore + prefix_C.seq2foldscore;
                    foldingscore += seq1_out_C[j1] + seq2_out_C[j2];
                }

                // final score
                if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN)
                    newscore = LOG_OF_ZERO;
                else
                    newscore = foldingscore + weight * alignscore;
            } else {
                // original
                int foldingscore = (prefix_C.seq1foldscore + prefix_C.seq2foldscore) + (cand.seq1foldscore + cand.seq2foldscore);
                float alignscore = xlog_mul(xlog_mul(prefix_C.alignscore, cand.alignscore), hmmalign.trans_probs[2][2]);
                
                if (alignscore <= LOG_OF_ZERO || foldingscore == VALUE_MIN) // DEBUG
                    newscore = LOG_OF_ZERO;
                else
                    newscore = foldingscore + weight * alignscore;
            }

            candidates.push_back(make_tuple(newscore, j1, 0, make_pair(i1, i2)));
        }
    }

    if (candidates.size() <= beam) return VALUE_FMIN;

    float threshold = quickselect(candidates, 0, candidates.size() - 1, candidates.size() - beam);
    for (auto &p : candidates) {
        float score = get<0>(p);
        if (score < threshold) {
            int j1 = get<1>(p);
            pair<int, int> key = get<3>(p);
            beamstep[j1].erase(key);
        }
    }

    return threshold;
}

void BeamSankoffParser::prepare(const vector<string> &seqs){
    num_seqs = seqs.size();
    sequences.clear();

    // preprocess: map nucleotides to integers
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
        sequences.push_back(seq);
        sum_len += seq.seq_len;
    }
    cout << "sum of seq length: " << sum_len << " " << seq1_len << " " << seq2_len << endl;

    // preprocess: get next posititon paired with i after j
    {    
        for (int i=0; i<num_seqs; i++){
            SeqObject* seq = &sequences[i];
            int seq_len = seq->seq_len;
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                seq->next_pair[nuci].resize(seq_len, -1);
                int next = -1;
                for (int nucj = seq_len-1; nucj>=1; --nucj) {
                    seq->next_pair[nuci][nucj] = next;
                    if (_allowed_pairs[nuci][seq->nucs[nucj]]) next = nucj;
                }
            }
        }
    }

    // preprocess: HMM align tool
    cout << endl;
    cout << "********** HMM alignment preprocessing **********" << endl;
    hmmalign.set(alnbeam, sequences[0].nucs, sequences[1].nucs);
    float similarity = hmmalign.viterbi_path(false);
    hmmalign.set_parameters_by_sim(similarity); // load new parameters
    
    // cout << endl;
    // float newsimilarity = hmmalign.viterbi_path(true); // get new similarity 
    // similarity = newsimilarity;
    // hmmalign.set_parameters_by_sim(similarity); // load new parameters

    float avg_width = hmmalign.envelope(similarity); // alignment envelope
    cout << endl;
    
    // update alignment weight
    // cout << "weight old: " << weight << endl;
    // weight *= (1 - 2 * avg_width / (seq1_len + seq2_len)); // avg_width / avg_seq_len / 2
    // // weight = 100 * similarity * (1 - 2 * avg_width / (seq1_len + seq2_len)); // avg_width / avg_seq_len / 2
    // cout << "weight new: " << weight << endl << endl;

    // save alignment backward score
    hmmalign.viterbi_path(true);
    hmmalign.viterbi_backward(true);
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
    cout << endl;

    // local alignment scores
    hmmalign.viterbi_path_all_locals();

    // preprocess: allocate space
    /********** allocate space preprocessing **********/
    bestH = new unordered_map<pair<int, int>, State, pair_hash>**[sum_len+1];
    bestP = new unordered_map<pair<int, int>, State, pair_hash>**[sum_len+1];
    bestMulti = new unordered_map<pair<int, int>, State, pair_hash>**[sum_len+1];
    bestM = new unordered_map<pair<int, int>, State, pair_hash>**[sum_len+1];

    bestM2 = new unordered_map<pair<int, int>, State, pair_hash>*[sum_len+1];
    for (int s = 1; s <= sum_len; ++s) {
        int minj1 = hmmalign.min_j1[s];
        int maxj1 = hmmalign.max_j1[s];
        if (maxj1 == -1) continue; // j1+j1==s not exist

        // cout << s << " " << minj1 << " " << maxj1 << endl;
        int interval = maxj1 - minj1 + 1;

        bestH[s] = new unordered_map<pair<int, int>, State, pair_hash>*[interval];
        bestH[s] = bestH[s] - minj1;

        bestP[s] = new unordered_map<pair<int, int>, State, pair_hash>*[interval];
        bestP[s] = bestP[s] - minj1;

        bestMulti[s] = new unordered_map<pair<int, int>, State, pair_hash>*[interval];
        bestMulti[s] = bestMulti[s] - minj1;

        bestM[s] = new unordered_map<pair<int, int>, State, pair_hash>*[interval];
        bestM[s] = bestM[s] - minj1;

        for (int j=minj1; j<=maxj1; j++) {
            bestH[s][j] = new unordered_map<pair<int, int>, State, pair_hash>[3];
            bestP[s][j] = new unordered_map<pair<int, int>, State, pair_hash>[3];
            bestMulti[s][j] = new unordered_map<pair<int, int>, State, pair_hash>[3];
            bestM[s][j] = new unordered_map<pair<int, int>, State, pair_hash>[3];
        }

        // bestM[s] = new unordered_map<pair<int, int>, State, pair_hash>[interval];
        // bestM[s] = bestM[s] - minj1;

        bestM2[s] = new unordered_map<pair<int, int>, State, pair_hash>[interval];
        bestM2[s] = bestM2[s] - minj1;
    }
    bestC.clear();
    bestC.resize(sum_len+1);

    // single sequence folding
    cout << endl;
    cout << "********** single sequence folding preprocessing **********" << endl;
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

    seq1_in_H.clear();
    seq1_in_P.clear();
    seq1_in_M.clear();
    seq1_in_M2.clear();
    seq1_in_Multi.clear();
    seq1_in_C.clear();
    seq2_in_H.clear();
    seq2_in_P.clear();
    seq2_in_M.clear();
    seq2_in_M2.clear();
    seq2_in_Multi.clear();
    seq2_in_C.clear();

    seq1_in_H.resize(seq1_len);
    seq1_in_P.resize(seq1_len);
    seq1_in_M.resize(seq1_len); 
    seq1_in_M2.resize(seq1_len);
    seq1_in_Multi.resize(seq1_len);
    seq1_in_C.resize(seq1_len);
    seq2_in_H.resize(seq2_len);
    seq2_in_P.resize(seq2_len);
    seq2_in_M.resize(seq2_len); 
    seq2_in_M2.resize(seq2_len);
    seq2_in_Multi.resize(seq2_len);
    seq2_in_C.resize(seq2_len);

    for (int i_seq=0; i_seq<2; i_seq++) {
        int seq_viterbi = cky_parser->parse(seqs[i_seq], NULL, NULL);
        if (i_seq == 0) seq1_mfe = seq_viterbi;
        else seq2_mfe = seq_viterbi;
        
        // only keep suboptimal base pairs
        // the maximum % change in free energy from the lowest free energy structure
        float min_score = seq_viterbi * (1 - max_energy_diff); // max_score must be positive
        cout << "min_score: " << min_score << endl;
        
        int seq_len = i_seq==0? seq1_len: seq2_len;
        for (int j=0; j<seq_len-1; j++) {
            vector<unordered_map<int, LFState>*> seq_in{&cky_parser->bestH[j], &cky_parser->bestP[j], &cky_parser->bestM[j], &cky_parser->bestM2[j], &cky_parser->bestMulti[j]};
            vector<unordered_map<int, LFState>*> seq_out{&cky_parser->bestH_beta[j], &cky_parser->bestP_beta[j], &cky_parser->bestM_beta[j], &cky_parser->bestM2_beta[j], &cky_parser->bestMulti_beta[j]};
            
            vector<unordered_map<int, int>*> seq_out_saved;
            if (i_seq==0)
                seq_out_saved = {&seq1_out_H[j+1], &seq1_out_P[j+1], &seq1_out_M[j+1], &seq1_out_M2[j+1], &seq1_out_Multi[j+1]};
            else
                seq_out_saved = {&seq2_out_H[j+1], &seq2_out_P[j+1], &seq2_out_M[j+1], &seq2_out_M2[j+1], &seq2_out_Multi[j+1]};

            // debug
            vector<unordered_map<int, int>*> seq_in_saved;
            if (i_seq==0)
                seq_in_saved = {&seq1_in_H[j+1], &seq1_in_P[j+1], &seq1_in_M[j+1], &seq1_in_M2[j+1], &seq1_in_Multi[j+1]};
            else
                seq_in_saved = {&seq2_in_H[j+1], &seq2_in_P[j+1], &seq2_in_M[j+1], &seq2_in_M2[j+1], &seq2_in_Multi[j+1]};

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

                for (auto i : valid_pos) {
                    (*seq_in_saved[k])[i+1] = beamins[i].score; // debug
                    (*seq_out_saved[k])[i+1] = beamout[i].score;
                }
            }

            // suffix C
            if (i_seq == 0) seq1_out_C[j+1] = cky_parser->bestC_beta[j].score;
            else seq2_out_C[j+1] = cky_parser->bestC_beta[j].score;
        }
    }
    delete cky_parser;
}


#ifdef multilign
void BeamSankoffParser::parse(const vector<string> &seqs, bool limited, const set<pair<int, int>> &allowed_pairs, vector<pair<int, int>> &out_pairs, int num_pairs){
#else
void BeamSankoffParser::parse(const vector<string> &seqs){
#endif

    /********** preprocessing **********/
    struct timeval parse_starttime, parse_endtime;
    gettimeofday(&parse_starttime, NULL);
    
    // allocate space and pre-computation
    seq1_len = seqs[0].size() + 1;
    seq2_len = seqs[1].size() + 1;
    prepare(seqs);

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    printf("seqs %d %d pre-processing time: %f seconds.\n", seq1_len, seq2_len, parse_elapsed_time);

    /********** core code **********/
    cout << endl;
    cout << "********** start computing **********" << endl;
    gettimeofday(&parse_starttime, NULL);

    SeqObject* seq1 = &sequences[0];
    SeqObject* seq2 = &sequences[1];

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
    bestC[0][0].set(0, 0, 0, 0, 0, MANNER_NONE, MANNER_NONE, 0.0, MANNER_ALN, 0.0);
    bestC[1][1].set(1, 0, 0, -v_score_external_unpaired(0, 0), 0, MANNER_NONE, MANNER_C_eq_C_plus_U_INS1, 1.0, MANNER_ALN, weight);
    bestC[1][0].set(0, 0, 0, 0, -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_INS2, 1.0, MANNER_ALN, weight);
    bestC[2][1].set(1, 0, 0, -v_score_external_unpaired(0, 0), -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_ALN, 0.0, MANNER_ALN, 0.0);

#else
    bestC[0][0].alnobj.set(0, 0, 0, 0, 0, MANNER_START, MANNER_START, xlog(1.0), MANNER_ALN, MANNER_ALN, weight*xlog(1.0));

    if (0 >= hmmalign.low_bounds[1] && 0 <= hmmalign.up_bounds[1]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_INS1, 1, 0, true);
        bestC[1][1].ins1obj.set(1, 0, 0, -v_score_external_unpaired(0, 0), 0, MANNER_START, MANNER_C_eq_C_plus_U_INS1, alignscore, MANNER_ALN, MANNER_INS1, weight*alignscore);
    }
    if (1 >= hmmalign.low_bounds[0] && 1 <= hmmalign.up_bounds[0]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_INS2, 0, 1, true);
        bestC[1][0].ins2obj.set(0, 0, 0, 0, -v_score_external_unpaired(0, 0), MANNER_START, MANNER_C_eq_C_plus_U_INS2, alignscore, MANNER_ALN, MANNER_INS2, weight*alignscore);
    }
    if (1 >= hmmalign.low_bounds[1] && 1 <= hmmalign.up_bounds[1]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_ALN, 1, 1, true);
        bestC[2][1].alnobj.set(1, 0, 0, -v_score_external_unpaired(0, 0), -v_score_external_unpaired(0, 0), MANNER_START, MANNER_C_eq_C_plus_U_ALN, alignscore, MANNER_ALN, MANNER_ALN, weight*alignscore);
    }
#endif


    // from left to right
    processMem_t mem;
    float threshold;
    for(int s = 1; s < seq1_len + seq2_len - 1; ++s) {
        // mem = GetProcessMemory();
        // cout << "s: " << s << " VmPeak: " << mem.VmPeak << endl;
        if (s % 10 == 1) cout << "s: " << s << endl;

        unordered_map<pair<int, int>, State, pair_hash>** beamH = bestH[s];
        unordered_map<pair<int, int>, State, pair_hash>** beamP = bestP[s];
        unordered_map<pair<int, int>, State, pair_hash>** beamMulti = bestMulti[s];
        unordered_map<pair<int, int>, State, pair_hash>** beamM = bestM[s];
        unordered_map<pair<int, int>, State, pair_hash>* beamM2 = bestM2[s];
        unordered_map<int, State3>& beamC = bestC[s];
      
        // hairpin candidates
        // for nucj push H(j, j_next)
        if (verbose) cout << "push H" << endl;
        for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
            int j2 = s - j1;
            if (j1 < 1 || j2 < 1 || j2 >= seq2->seq_len - 1) continue; // boundary case

            // hmm constraint, left bracket
            if (j2 > hmmalign.up_bounds[j1] || j2 < hmmalign.low_bounds[j1]) continue;

            int j1next = j1;
            pair<int, int> newscore1, newscore2;
            while (j1next != -1) {
                newscore1 = hairpinScore(j1, j1next, seq1);
                j1next = newscore1.first;
                if (j1next == -1) break;

                // single seq folding 
                if (seq1_out_H[j1next].find(j1) == seq1_out_H[j1next].end()) continue;

                int j2next = j2;
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
                                        bestH[j1next+j2next][j1next][2][make_pair(j1, j2)],
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_ALN, 
                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
#ifndef dynalign
                    float ins1score = get_hmm_score(j1, j1next, j2+1, j2next-1, 0);
                    if (ins1score > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1next+j2next][j1next][0][make_pair(j1, j2)],
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_INS1, 
                                        ins1score, MANNER_INS1, MANNER_INS1, weight, verbose);

                    float ins2score = get_hmm_score(j1+1, j1next-1, j2, j2next, 1);
                    if (ins2score > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1next+j2next][j1next][1][make_pair(j1, j2)],
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
            if (beam > 0) beam_prune(beamH, s, seq1_out_H, seq2_out_H);

            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (int m=0; m<3; m++) {
                    Manner next_manner;
                    switch (m) {
                        case 0:
                            next_manner = MANNER_HAIRPIN_INS1;
                            break;
                        case 1:
                            next_manner = MANNER_HAIRPIN_INS2;
                            break;
                        case 2:
                            next_manner = MANNER_HAIRPIN_ALN;
                            break;
                    }
                    for (auto &item: beamH[j1][m]) {
                        State &state = item.second;
                        if (state.manner == MANNER_NONE) continue;

                        int i1 = state.i1;
                        int i2 = state.i2;

                        assert (i1 > 0);
                        assert (i2 > 0);

                        // single seq folding
                        // H->P(i, j)
                        if (seq1_out_P[j1].find(i1) != seq1_out_P[j1].end() && seq2_out_P[j2].find(i2) != seq2_out_P[j2].end()) {
                            if (verbose) cout << "H to P: " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.endHMMstate << endl; 
                            
                            update_if_better(i1, j1, i2, j2,
                                            beamP[j1][m][make_pair(i1, i2)],  // beamPINS1[get_keys(j1, i1, i2)], 
                                            state.seq1foldscore, state.seq2foldscore,
                                            state.manner, next_manner, 
                                            state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                    }
                } // end of beamH[j1]
                }
            } // end of s
        } // end of beam H

        // beam of Multi
        // for every state in Multi[j]
        //   1. extend (i, j) to (i, jnext)
        //   2. generate P (i, j)
        if (verbose) cout << "beam of Multi" << endl;
        {
            if (beam > 0) threshold = beam_prune(beamMulti, s, seq1_out_Multi, seq2_out_Multi);

            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (int m=0; m<3; m++) {
                    Manner next_manner;
                    switch (m) {
                        case 0:
                            next_manner = MANNER_P_eq_MULTI_INS1;
                            break;
                        case 1:
                            next_manner = MANNER_P_eq_MULTI_INS2;
                            break;
                        case 2:
                            next_manner = MANNER_P_eq_MULTI_ALN;
                            break;
                    }
                    for (auto &item : beamMulti[j1][m]) {
                        State &state = item.second;
                        if (state.manner == MANNER_NONE) continue;

                        int i1 = state.i1;
                        int i2 = state.i2;
                        assert (i1 > 0); // DEBUG commet it later
                        assert (i2 > 0);

                        // 1. extend (i, j) to (i, jnext)
                        int j1next = j1;
                        int j2next = j2;
                        // while (j1next != -1 && j2next != -1 && j1next <= min(seq1_len, j1+10) && j2next <= min(seq2_len, j2+10))
                        // {
                        pair<int, int> result1 = multiloopUnpairedScore(i1, j1next, seq1);
                        pair<int, int> result2 = multiloopUnpairedScore(i2, j2next, seq2);

                        j1next = get<0>(result1);
                        j2next = get<0>(result2);

                        int new_seq1_l2 = state.trace1.paddings.l2 + j1next - j1;
                        int new_seq2_l2 = state.trace2.paddings.l2 + j2next - j2;

                        HMMManner premanner = state.endHMMstate;

                        // if (new_seq1_l2 <= SINGLE_MAX_LEN && new_seq2_l2 <= SINGLE_MAX_LEN)
                        {
                            // lengths of aligned internal loops are close
                            // int loop_diff = new_seq2_l2 - new_seq1_l2;
                            // if (loop_diff < -10) {
                            //     q2 = seq2->next_pair[nucp2][q2];
                            //     continue;
                            // }
                            // if (loop_diff > 10) break;

                            if (verbose) cout << "beamMulti jnext: " << i1  <<  " " << j1 << " " << i2 << " " << j2 << " "  << j1next  <<  " " << j2next << endl; 
                            
                            float alignscore;
                            switch (premanner)
                            {
                                case MANNER_ALN:
                                {
                                    if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
#ifdef dynalign
                                        // do not apply alignment constraint on right side
                                        aln_value_type alignscore = state.alignscore - abs(state.trace1.paddings.l2 - state.trace2.paddings.l2) + abs(new_seq1_l2 - new_seq2_l2);
#else
                                        alignscore = get_hmm_score_right(j1, j1next, j2, j2next, 2, 2, true); // true
                                        alignscore = xlog_mul(state.alignscore, alignscore);

#endif
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1next, i2, j2next, bestMulti[j1next+j2next][j1next][2][make_pair(i1, i2)], // bestMulti[j1next+j2next][get_keys(j1next, i1, i2)].alnobj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                                    }
                                    break;
                                }
                                case MANNER_INS1:
                                {
                                    if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
                                        alignscore = get_hmm_score_right(j1, j1next, j2, j2next-1, 0, 0, true); // true
                                        alignscore = xlog_mul(state.alignscore, alignscore);
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1next, i2, j2next, bestMulti[j1next+j2next][j1next][0][make_pair(i1, i2)], // bestMulti[j1next+j2next][get_keys(j1next, i1, i2)].ins1obj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                                    }
                                    if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) {
                                        alignscore = get_hmm_score_right(j1, j1next, j2, j2-1, 0, 0, true); // true
                                        alignscore = xlog_mul(state.alignscore, alignscore);

                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1next, i2, j2, bestMulti[j1next+j2][j1next][0][make_pair(i1, i2)], // bestMulti[j1next+j2][get_keys(j1next, i1, i2)].ins1obj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                            alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                                    }
                                    break;
                                }  
                                case MANNER_INS2:   
                                {
                                    if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
                                        alignscore = get_hmm_score_right(j1, j1next-1, j2, j2next, 1, 1, true); // true
                                        alignscore = xlog_mul(state.alignscore, alignscore);
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1next, i2, j2next, bestMulti[j1next+j2next][j1next][1][make_pair(i1, i2)], // bestMulti[j1next+j2next][get_keys(j1next, i1, i2)].ins2obj
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                                    }
                                    if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) {
                                        alignscore = get_hmm_score_right(j1, j1-1, j2, j2next, 1, 1, true); // true
                                        alignscore = xlog_mul(state.alignscore, alignscore);

                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1, i2, j2next, bestMulti[j1+j2next][j1][1][make_pair(i1, i2)], // bestMulti[j1+j2next][get_keys(j1, i1, i2)].ins2obj,
                                                            state.seq1foldscore,
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                            state.trace1.paddings.l1, state.trace1.paddings.l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                                    }
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
                                if (verbose) cout << "Multi to P: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;

                                int score1 = multiloop2Pscore(i1, j1, seq1);
                                int score2 = multiloop2Pscore(i2, j2, seq2);
                            
                                update_if_better(i1, j1, i2, j2, beamP[j1][m][make_pair(i1,i2)], //beamP[item.first].alnobj,
                                                state.seq1foldscore + score1, 
                                                state.seq2foldscore + score2,
                                                state.manner, next_manner,
                                                state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                            }
                        } // 2. generate P (i, j)
                    } // end of items
                }
            } // end of j1
        } // end beam of Multi

        // beam of P
        // for every state in P[j]
        //   1. generate new helix/bulge
        //   2. M = P
        //   3. M2 = M + P
        //   4. C = C + P
        if (verbose) cout << "beam of P" << endl;
        {
            if (beam > 0) threshold = beam_prune(beamP, s, seq1_out_P, seq2_out_P);

            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (int m=0; m<3; m++) {
                    for (auto &item: beamP[j1][m]) {
                        State& state = item.second;
                        if (state.manner == MANNER_NONE) continue;

                        int i1 = state.i1;
                        int i2 = state.i2;
                        assert (i1 > 0);
                        assert (i2 > 0);

                        int nuci1 = seq1->nucs[i1];
                        int nuci1_1 = seq1->nucs[i1 - 1];
                        int nucj1 = seq1->nucs[j1];
                        int nucj1p1 = seq1->nucs[j1 + 1];

                        int nuci2 = seq2->nucs[i2];
                        int nuci2_1 = seq2->nucs[i2 - 1];
                        int nucj2 = seq2->nucs[j2];
                        int nucj2p1 = seq2->nucs[j2 + 1];

                        HMMManner premanner = state.endHMMstate;

                        // 1. generate new helix / single_branch
                        // new state is of shape p..i..j..q
                        // Note: p >= 1
                        {
                            // for (int p1 = i1 - 1; p1 >= max(i1 - SINGLE_MAX_LEN, 1); --p1) {
                            int min_p1 = max(i1 - MAX_LOOP_LEN - 1, 1);
                            for (int p1 = i1; p1 >= min_p1; --p1) {
                                if (m == 1 && p1 != i1) continue; // Noted. 

                                int nucp1 = seq1->nucs[p1];
                                int nucp1p1 = seq1->nucs[p1 + 1];
                                int q1 = (p1==i1) ? j1 : seq1->next_pair[nucp1][j1];

                                int i1_p1 = i1 - p1;
                                // while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                                while (q1 != -1 &&  (q1 - j1 - 1 <= MAX_LOOP_LEN)) {
                                    if ((p1==i1) && (q1!=j1)) break;
                                    // single seq folding
                                    if (seq1_out_P[q1].find(p1) == seq1_out_P[q1].end()) {
                                        q1 = seq1->next_pair[nucp1][q1];
                                        continue;
                                    }
                                    int q1_j1 = q1 - j1;

    // #ifdef multilign
    //                                         if (limited && allowed_pairs.find(make_pair(p1, q1)) == allowed_pairs.end()) {
    //                                             q1 = seq1->next_pair[nucp1][q1];
    //                                             continue;
    //                                         }
    // #endif
                                    int nucq1 = seq1->nucs[q1];
                                    int nucq1_1 = seq1->nucs[q1 - 1];

                                    int p2p1 = (p1 == i1) ? 0 : P2PScore(p1,q1,i1,j1,nucp1,nucp1p1,nucq1_1,nucq1,nuci1_1,nuci1,nucj1,nucj1p1); 

                                    int max_p2 = min(i2, hmmalign.up_bounds[p1]);
                                    int min_p2 = max(hmmalign.low_bounds[p1], max(1, i2 - MAX_LOOP_LEN - 1));

                                    // max_p2 = min(max_p2, i2-i1_p1+10); // lengths of aligned internal loops are close. p2 = i2+i1-p1 +- 10
                                    // min_p2 = max(min_p2, i2-i1_p1-10); // lengths of aligned internal loops are close 

                                    for (int p2 = max_p2; p2 >= min_p2; --p2) { // seq2 for loop
                                        if ((p1 == i1) && (p2 == i2)) continue;
                                        if (m == 0 && p2 != i2) continue; // Noted. 

                                        int nucp2 = seq2->nucs[p2];
                                        int nucp2p1 = seq2->nucs[p2 + 1];
                                        int q2 = (p2 == i2) ? j2 : seq2->next_pair[nucp2][j2];

                                        // speed up
                                        int i2_p2 = i2 - p2;
                                        if (q2 > hmmalign.up_bounds[q1] || q2 == -1 || (q2 - j2 - 1 > MAX_LOOP_LEN)) continue;
                                        if (q2 < hmmalign.low_bounds[q1])
                                            q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1] - 1]; 

                                        // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                        while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= MAX_LOOP_LEN)) {
                                            if ((p2 == i2) && (q2 != j2)) break;
                                            // single seq folding
                                            if (seq2_out_P[q2].find(p2) == seq2_out_P[q2].end()) {
                                                q2 = seq2->next_pair[nucp2][q2];
                                                continue;
                                            }
                                            
                                            // lengths of aligned internal loops are close
                                            // int loop_diff = q2 - j2 - q1_j1;
                                            // if (loop_diff < -10) {
                                            //     q2 = seq2->next_pair[nucp2][q2];
                                            //     continue;
                                            // }
                                            // if (loop_diff > 10) break;

                                            // TODO: redundant calculation
                                            int nucq2 = seq2->nucs[q2];
                                            int nucq2_1 = seq2->nucs[q2 - 1];
                                            int p2p2 = (p2 == i2) ? 0 : P2PScore(p2,q2,i2,j2,nucp2,nucp2p1,nucq2_1,nucq2,nuci2_1,nuci2,nucj2,nucj2p1);

                                            if (verbose) cout << "P2P: " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.endHMMstate << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;
                                            
                                            pair<float, HMMManner> pre_alignscore, post_alignscore;
                                            float pre_align_trans, post_align_trans, alignscore;

                                            
                                            // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj)
                                            //  si == sj == ALN/INS1/INS2
                                            // if (p1 != i1 && p2 != i2) 
                                            {
#ifdef dynalign                                     
                                                aln_value_type alignscore = ALN_VALUE_MIN;
                                                if (p2>=hmmalign.low_bounds[p1] && p2<=hmmalign.up_bounds[p1] && q2>=hmmalign.low_bounds[q1] && q2<=hmmalign.up_bounds[q1]){
                                                    alignscore = state.alignscore + abs(i1 - p1 - 1) + abs(q1 - j1 - 1);
                                                }
#else
                                                // si == sj == ALN
                                                // align cost = (p, i) + (i, j) + (j, q)
                                                pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, premanner-1);
                                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2, premanner-1, 2);
                                                
                                                alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);
#endif                                        
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][q1][2][make_pair(p1, p2)],//bestP[q1+q2][get_keys(q1, p1, p2)].alnobj, 
                                                                    state.seq1foldscore + p2p1,
                                                                    state.seq2foldscore + p2p2,
                                                                    state.manner, MANNER_SINGLE_ALN, 
                                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                            }
                                            // si == sj == INS1
                                            // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                            // if (p1 < i1) 
                                            {
#ifdef dynalign
                                                if (p1==(i1-1) && p2==(i2-1) && q1==(j1+1) && q2==(j2+1))
                                                { // ALN -> INS1/INS2, not allow continual inserted base pairs
                                                    aln_value_type alignscore = ALN_VALUE_MIN;
                                                    if (p2+1>=hmmalign.low_bounds[p1] && p2+1<=hmmalign.up_bounds[p1] && q2-1>=hmmalign.low_bounds[q1] && q2-1<=hmmalign.up_bounds[q1]){
                                                        alignscore = alnstate.alignscore + abs(i1 - p1 - i2 + p2) + abs(q1 - j1 - q2 + j2) + 2;
                                                    }
                                                }
#else
                                                pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, premanner-1);
                                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, premanner-1, 0);
                                                
                                                alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans); 
#endif   
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][q1][0][make_pair(p1, p2)], // bestP[q1+q2][get_keys(q1, p1, p2)].ins1obj, 
                                                                    state.seq1foldscore + p2p1,
                                                                    state.seq2foldscore + p2p2,
                                                                    state.manner, MANNER_SINGLE_INS1, 
                                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                                            }
                                            // si == sj == INS2
                                            // align cost = (p, i) + (i, j) + (j, q);
                                            // if (p2 < i2) 
                                            {
#ifdef dynalign
                                                if (p1==(i1-1) && p2==(i2-1) && q1==(j1+1) && q2==(j2+1))
                                                { // ALN -> INS1/INS2, not allow continual inserted base pairs
                                                    aln_value_type alignscore = ALN_VALUE_MIN;
                                                    if (p2>=hmmalign.low_bounds[p1+1] && p2<=hmmalign.up_bounds[p1+1] && q2>=hmmalign.low_bounds[q1-1] && q2<=hmmalign.up_bounds[q1-1]){
                                                        alignscore = state.alignscore + abs(i1 - p1 - i2 + p2) + abs(q1 - j1 - q2 + j2) + 2;
                                                    }
                                                }
#else
                                                pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, premanner-1);
                                                post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, premanner-1, 1);
                                                
                                                alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);
#endif 
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestP[q1+q2][q1][1][make_pair(p1, p2)], // bestP[q1+q2][get_keys(q1, p1, p2)].ins2obj,
                                                                    state.seq1foldscore + p2p1,
                                                                    state.seq2foldscore + p2p2,
                                                                    state.manner, MANNER_SINGLE_INS2, 
                                                                    static_cast<char>(i1 - p1), q1 - j1,
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

                        if (m != 2) continue;

                        // 2. M = P
                        // accessible pairs in multiloop must be aligned
                        // singe seq folding
                        if (seq1_out_M[j1].find(i1) != seq1_out_M[j1].end() && seq2_out_M[j2].find(i2) != seq2_out_M[j2].end()) 
                        {
                            if (verbose) cout << "P to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                            int newscore1 = branch_score(i1, j1, seq1);
                            int newscore2 = branch_score(i2, j2, seq2);
                            update_if_better(i1, j1, i2, j2, beamM[j1][2][make_pair(i1, i2)],
                                            state.seq1foldscore + newscore1,
                                            state.seq2foldscore + newscore2, 
                                            state.manner, MANNER_M_eq_P,
                                            state.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                        } // 2. M = P

                        // 3. M2 = M + P 
                        // accessible pairs in multiloop must be aligned
                        {
                            int k1 = i1 - 1;
                            int k2 = i2 - 1;

                            if (k1 > 1 && k2 > 1 && k2 >= hmmalign.low_bounds[k1] && k2 <= hmmalign.up_bounds[k1]) {
                                int newscore1 = branch_score(i1, j1, seq1);
                                int newscore2 = branch_score(i2, j2, seq2); 

                                for (int mm=0; mm<3; mm++) {
                                    if (!bestM[k1+k2][k1][mm].empty()) {
                                        for (auto &m : bestM[k1+k2][k1][mm]) {
                                            State mstate = m.second; 
                                            int newi1 = mstate.i1;
                                            int newi2 = mstate.i2;

                                            // single seq folding 
                                            if (newi1 > 0 && newi2 > 0 && seq1_out_M2[j1].find(newi1) != seq1_out_M2[j1].end() && seq2_out_M2[j2].find(newi2) != seq2_out_M2[j2].end()) 
                                            {
                                                if (verbose) cout << "M2=M+P: " << i1 << " " << j1 << " "  << i2 << " " << j2 <<  " " << newi1 << " " << j1 << " "  << newi2 << " " << j2  << endl;
#ifdef dynalign
                                                aln_value_type alignscore = state.alignscore + mstate.alignscore;
#else
                                                float alignscore = xlog_mul(mstate.alignscore, hmmalign.trans_probs[mstate.endHMMstate-1][2]); 
                                                alignscore = xlog_mul(alignscore, state.alignscore);
#endif 
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(newi1, j1, newi2, j2, beamM2[j1][make_pair(newi1, newi2)], // beamM2[get_keys(j1, newi1, newi2)],
                                                                    mstate.seq1foldscore + newscore1 + state.seq1foldscore,
                                                                    mstate.seq2foldscore + newscore2 + state.seq2foldscore,
                                                                    mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                                    alignscore, mstate.startHMMstate, MANNER_ALN, weight, verbose);
                                            }
                                        }
                                    }
                                }
                            }
                        } // 3. M2 = M + P 

                        // 4. C = C + P
                        // external pairs must be aligned
                        {   
                            int k1 = i1 - 1;
                            int k2 = i2 - 1;

                            if (k1 >= 0 && k2 >= 0 && bestC[k1+k2].find(k1) != bestC[k1+k2].end()) {
                                // if (bestC[k1+k2].find(k1) == bestC[k1+k2].end()) continue;

                                if (verbose) cout << "C+P: "<< i1 << " " << j1 << " "  << i2 << " " << j2 << endl;

                                int newscore1 = state.seq1foldscore + external_paired_score(k1, j1, seq1);
                                int newscore2 = state.seq2foldscore + external_paired_score(k2, j2, seq2);

                                {
                                    State& prefix_C = bestC[k1+k2][k1].alnobj; // bestC[k1][k2];
                                    if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                        int alignscore = state.alignscore + prefix_C.alignscore;
    #else
                                        float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[2][2]); 
                                        alignscore = xlog_mul(alignscore, state.alignscore); 
    #endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(0, j1, 0, j2, beamC[j1].alnobj,
                                                            prefix_C.seq1foldscore + newscore1,
                                                            prefix_C.seq2foldscore + newscore2, 
                                                            prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                            k1, k2,
                                                            alignscore, prefix_C.startHMMstate, MANNER_ALN, weight, verbose);
                                    }
                                }
                                {
                                    State& prefix_C = bestC[k1+k2][k1].ins1obj;
                                    if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                        int alignscore = state.alignscore + prefix_C.alignscore;
    #else
                                        float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[0][2]); 
                                        alignscore = xlog_mul(alignscore, state.alignscore); 
    #endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(0, j1, 0, j2, beamC[j1].alnobj,
                                                            prefix_C.seq1foldscore + newscore1,
                                                            prefix_C.seq2foldscore + newscore2, 
                                                            prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                            k1, k2,
                                                            alignscore, prefix_C.startHMMstate, MANNER_ALN, weight, verbose);
                                    }
                                }
                                {
                                    State& prefix_C = bestC[k1+k2][k1].ins2obj;
                                    if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                        int alignscore = state.alignscore + prefix_C.alignscore;
    #else
                                        float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[1][2]); 
                                        alignscore = xlog_mul(alignscore, state.alignscore); 
    #endif 
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(0, j1, 0, j2, beamC[j1].alnobj,
                                                            prefix_C.seq1foldscore + newscore1,
                                                            prefix_C.seq2foldscore + newscore2, 
                                                            prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                            k1, k2,
                                                            alignscore, prefix_C.startHMMstate, MANNER_ALN, weight, verbose);
                                    }
                                }   
                            }
                        } // 4. C = C + P
                    } // beamP
                }
            } 
        }

        // branch insertion 
        // fully insertion P in seq1 
        // if (false) 
        {
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j1 < 1 || j2 < 1 || j2 >= seq2->seq_len - 1) continue; // boundary case

                // hmm constraint, TODO
                if (j2 > hmmalign.up_bounds[j1] || j2 < hmmalign.low_bounds[j1]) continue;

                for(auto& item : seq1_in_P[j1]) {
                    int i1 = item.first;
                    int seq1foldscore = item.second;
        
                    if (j2 >= hmmalign.low_bounds[i1] && j2 <= hmmalign.up_bounds[i1]) {
                        // P1 to M
                        // single sequnce folding
                        if (seq1_out_M[j1].find(i1) != seq1_out_M[j1].end()) 
                        {
                            if (verbose) cout << "P1 to M: " << i1 << " " << j1 << " " << j2 << endl;

                            int newscore1 = branch_score(i1, j1, seq1);
                            float alignscore = get_hmm_score_right(i1, j1, j2, j2, 0, 0); // fully inserted i1-j1
                            
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(i1, j1, j2+1, j2, beamM[j1][0][make_pair(i1, j2+1)],
                                                seq1foldscore + newscore1,
                                                0, 
                                                MANNER_NONE, MANNER_M_eq_P1,
                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                        }
                        
                        // M + P1
                        int k1 = i1 - 1;
                        int k2 = j2;
                        if (k1 > 1 && k2 > 1 && k2 >= hmmalign.low_bounds[k1] && k2 <= hmmalign.up_bounds[k1]) {
                            int newscore1 = branch_score(i1, j1, seq1);

                            for (int mm=0; mm<3; mm++) {
                                if (!bestM[k1+k2][k1][mm].empty()) {
                                    for (auto &m : bestM[k1+k2][k1][mm]) {
                                        // ALN state
                                        {                                     
                                            State mstate = m.second; 
                                            int newi1 = mstate.i1;
                                            int newi2 = mstate.i2;

                                            // single sequnce folding
                                            if (newi1>0 && newi2>0 && seq1_out_M2[j1].find(newi1) != seq1_out_M2[j1].end() && seq2_out_M2[j2].find(newi2) != seq2_out_M2[j2].end()) 
                                            {
                                                if (verbose) cout << "M2=M+P1: " << i1 << " " << j1 << " "  << j2 << " " << j2 <<  " " << newi1 << " " << j1 << " "  << newi2 << " " << j2  << endl;
// #ifdef dynalign
//                                              aln_value_type alignscore = alnstate.alignscore + mstate.alignscore;
// #else
                                                float alignscore = get_hmm_score_right(i1, j1, j2, j2, mstate.endHMMstate-1, 0); // fully inserted i1-j1
                                                alignscore = xlog_mul(mstate.alignscore, alignscore);
// #endif 
                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(newi1, j1, newi2, j2, beamM2[j1][make_pair(newi1, newi2)],
                                                                    mstate.seq1foldscore + newscore1 + seq1foldscore,
                                                                    mstate.seq2foldscore,
                                                                    mstate.manner, MANNER_M2_eq_M_plus_P1, k1, k2,
                                                                    alignscore, mstate.startHMMstate, MANNER_INS1, weight, verbose);     
                                            }
                                        }
                                    }  
                                }
                            }
                        } // M + P1

                        // C + P1
                        if (false && k1 >= 0 && k2 >= 0 && bestC[k1+k2].find(k1) != bestC[k1+k2].end()) {
                            if (verbose) cout << "C+P1: "<< i1 << " " << j1 << " "  << j2 << " " << j2 << endl;

                            int newscore1 = external_paired_score(k1, j1, seq1);

                            {
                                State& prefix_C = bestC[k1+k2][k1].alnobj; // bestC[k1][k2];
                                if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                    int alignscore = alnstate.alignscore + prefix_C.alignscore;
    #else
                                    float alignscore = get_hmm_score_right(i1, j1, j2, j2, 2, 0); // fully inserted i1--j1
                                    alignscore = xlog_mul(prefix_C.alignscore, alignscore);
    #endif 
                                    if (alignscore > LOG_OF_ZERO)
                                        update_if_better(0, j1, 0, j2, beamC[j1].ins1obj,
                                                        prefix_C.seq1foldscore + newscore1 + seq1foldscore,
                                                        prefix_C.seq2foldscore, 
                                                        prefix_C.manner, MANNER_C_eq_C_plus_P1,
                                                        k1, k2,
                                                        alignscore, prefix_C.startHMMstate, MANNER_INS1, weight, verbose);
                                }
                            }
                            {
                                State& prefix_C = bestC[k1+k2][k1].ins1obj;
                                if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                    int alignscore = alnstate.alignscore + prefix_C.alignscore;
    #else
                                    float alignscore = get_hmm_score_right(i1, j1, j2, j2, 0, 0); // fully inserted i1--j1
                                    alignscore = xlog_mul(prefix_C.alignscore, alignscore);
    #endif 
                                    if (alignscore > LOG_OF_ZERO)
                                        update_if_better(0, j1, 0, j2, beamC[j1].ins1obj,
                                                        prefix_C.seq1foldscore + newscore1 + seq1foldscore,
                                                        prefix_C.seq2foldscore,
                                                        prefix_C.manner, MANNER_C_eq_C_plus_P1,
                                                        k1, k2,
                                                        alignscore, prefix_C.startHMMstate, MANNER_INS1, weight, verbose);
                                }
                            }
                            {
                                State& prefix_C = bestC[k1+k2][k1].ins2obj;
                                if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                    int alignscore = alnstate.alignscore + prefix_C.alignscore;
    #else
                                    float alignscore = get_hmm_score_right(i1, j1, j2, j2, 1, 0); // fully inserted i1--j1
                                    alignscore = xlog_mul(prefix_C.alignscore, alignscore);
    #endif 
                                    if (alignscore > LOG_OF_ZERO)
                                        update_if_better(0, j1, 0, j2, beamC[j1].ins1obj,
                                                        prefix_C.seq1foldscore + newscore1 + seq1foldscore,
                                                        prefix_C.seq2foldscore, 
                                                        prefix_C.manner, MANNER_C_eq_C_plus_P1,
                                                        k1, k2,
                                                        alignscore, prefix_C.startHMMstate, MANNER_INS1, weight, verbose);
                                }
                            }
                        }  // C + P1
                    } // valid j2
                } // traverse i1

                // fully inertion P in seq2
                for(auto& item : seq2_in_P[j2]) {
                    int i2 = item.first;
                    int seq2foldscore = item.second;

                    if (i2 >= hmmalign.low_bounds[j1] && i2 <= hmmalign.up_bounds[j1]) {
                        // P2 to M
                        // single sequnce folding
                        if (seq2_out_M[j2].find(i2) != seq2_out_M[j2].end()) 
                        {
                            int newscore2 = branch_score(i2, j2, seq2);
                            float alignscore = get_hmm_score_right(j1, j1, i2, j2, 1, 1); // fully inserted i2--j2

                            if (verbose) cout << "P2 to M: " << j1 << " " << i2 << " " << j2 << endl;

                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(j1+1, j1, i2, j2, beamM[j1][1][make_pair(j1+1, i2)],
                                                0,
                                                seq2foldscore + newscore2,
                                                MANNER_NONE, MANNER_M_eq_P2,
                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                        }
                
                        // M + P2
                        int k1 = j1;
                        int k2 = i2 - 1;
                        if (k1 > 1 && k2 > 1 && k2 >= hmmalign.low_bounds[k1] && k2 <= hmmalign.up_bounds[k1]) {
                            int newscore2 = branch_score(i2, j2, seq2);

                            for (int mm=0; mm<3; mm++) {
                                if (!bestM[k1+k2][k1][mm].empty()) {
                                    for (auto &m : bestM[k1+k2][k1][mm]) {
                                        State mstate = m.second;
                                        int newi1 = mstate.i1;
                                        int newi2 = mstate.i2;

                                        // single sequnce folding
                                        if (seq1_out_M2[j1].find(newi1) != seq1_out_M2[j1].end() && seq2_out_M2[j2].find(newi2) != seq2_out_M2[j2].end()) 
                                        {
                                            if (verbose) cout << "M2=M+P2: " << j1 << " " << j1 << " "  << i2 << " " << j2 <<  " " << newi1 << " " << j1 << " "  << newi2 << " " << j2  << endl;
// #ifdef dynalign
//                                          aln_value_type alignscore = alnstate.alignscore + mstate.alignscore;
// #else
                                            float alignscore = get_hmm_score_right(j1, j1, i2, j2, mstate.endHMMstate-1, 1); // fully inserted i2--j2
                                            alignscore = xlog_mul(mstate.alignscore, alignscore); 
// #endif 
                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(newi1, j1, newi2, j2, beamM2[j1][make_pair(newi1, newi2)],
                                                                mstate.seq1foldscore, 
                                                                mstate.seq2foldscore + newscore2 + seq2foldscore,
                                                                mstate.manner, MANNER_M2_eq_M_plus_P2, k1, k2,
                                                                alignscore, mstate.startHMMstate, MANNER_INS2, weight, verbose);
                                        }
                                    }
                                }
                            }
                        } // M + P2

                        // C + P2
                        if (false && k1 >= 0 && k2 >= 0 && bestC[k1+k2].find(k1) != bestC[k1+k2].end()) {
                            if (verbose) cout << "C+P2: "<< j1 << " " << j1 << " "  << i2 << " " << j2 << endl;
                            int newscore2 = external_paired_score(k2, j2, seq2);
                            {
                                State& prefix_C = bestC[k1+k2][k1].alnobj; // bestC[k1][k2];
                                if (prefix_C.manner != MANNER_NONE) {
#ifdef dynalign
                                    int alignscore = alnstate.alignscore + prefix_C.alignscore;
#else
                                    float alignscore = get_hmm_score_right(j1, j1, i2, j2, 2, 1); // fully inserted i2--j2
                                    alignscore = xlog_mul(prefix_C.alignscore, alignscore);
#endif 
                                    if (alignscore > LOG_OF_ZERO)
                                        update_if_better(0, j1, 0, j2, beamC[j1].ins2obj,
                                                        prefix_C.seq1foldscore, 
                                                        prefix_C.seq2foldscore + newscore2 + seq2foldscore,
                                                        prefix_C.manner, MANNER_C_eq_C_plus_P2,
                                                        k1, k2,
                                                        alignscore, prefix_C.startHMMstate, MANNER_INS2, weight, verbose);
                                }
                            }
                            {
                                State& prefix_C = bestC[k1+k2][k1].ins1obj;
                                if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                    int alignscore = alnstate.alignscore + prefix_C.alignscore;
    #else
                                    float alignscore = get_hmm_score_right(j1, j1, i2, j2, 0, 1); // fully inserted i2--j2
                                    alignscore = xlog_mul(prefix_C.alignscore, alignscore);
    #endif 
                                    if (alignscore > LOG_OF_ZERO)
                                        update_if_better(0, j1, 0, j2, beamC[j1].ins2obj,
                                                        prefix_C.seq1foldscore, 
                                                        prefix_C.seq2foldscore + newscore2 + seq2foldscore,
                                                        prefix_C.manner, MANNER_C_eq_C_plus_P2,
                                                        k1, k2,
                                                        alignscore, prefix_C.startHMMstate, MANNER_INS2, weight, verbose);
                                }
                            }
                            {
                                State& prefix_C = bestC[k1+k2][k1].ins2obj;
                                if (prefix_C.manner != MANNER_NONE) {
    #ifdef dynalign
                                    int alignscore = alnstate.alignscore + prefix_C.alignscore;
    #else
                                    float alignscore = get_hmm_score_right(j1, j1, i2, j2, 1, 1); // fully inserted i2--j2
                                    alignscore = xlog_mul(prefix_C.alignscore, alignscore);
    #endif 
                                    if (alignscore > LOG_OF_ZERO)
                                        update_if_better(0, j1, 0, j2, beamC[j1].ins2obj,
                                                        prefix_C.seq1foldscore, 
                                                        prefix_C.seq2foldscore + newscore2 + seq2foldscore,
                                                        prefix_C.manner, MANNER_C_eq_C_plus_P2,
                                                        k1, k2,
                                                        alignscore, prefix_C.startHMMstate, MANNER_INS2, weight, verbose);
                                }
                            }
                        }  // C + P2
                    } // valid j1
                } // traverse i2
            } // traverse j1
        }

        // beam of M2
        // for every state in M2[j]
        //   1. multi-loop  (by extending M2 on the left)
        //   2. M = M2
        if (verbose) cout << "beam of M2" << endl;
        {
            // sort_keys(beamM2, keys); // cube pruning? 
            if (beam > 0) threshold = beam_prune(beamM2, s, seq1_out_M2, seq2_out_M2);
            
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (auto &item : beamM2[j1]) {
                    State state = item.second; 
                    int i1 = state.i1;
                    int i2 = state.i2;

                    assert (i1 > 0);
                    assert (i2 > 0);

                    // 2. M = M2
                    // single seq folding
                    if (seq1_out_M[j1].find(i1) != seq1_out_M[j1].end() && seq2_out_M[j2].find(i2) != seq2_out_M[j2].end()) 
                    {
                        if (verbose) cout << "M2 to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                        int nextm = state.endHMMstate - 1;
                        update_if_better(i1, j1, i2, j2, beamM[j1][nextm][make_pair(i1, i2)],
                                        state.seq1foldscore,
                                        state.seq2foldscore,
                                        state.manner, MANNER_M_eq_M2,
                                        state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                    } // 2. M = M2

                    // 1. multi-loop
                    {
                        for (int p1 = i1-1; p1 >= max(i1 - SINGLE_MAX_LEN - 1, 1); --p1) {
                            int nucp1 = seq1->nucs[p1];
                            int q1 = seq1->next_pair[nucp1][j1];
                            
                            int i1_p1 = i1 - p1;
                            if (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN) && seq1_out_Multi[q1].find(p1) != seq1_out_Multi[q1].end()) {
                            // while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                            // while (q1 != -1 && (q1 - j1 - 1 <= 5)) {
                                // single seq folding
                                // if (seq1_out_Multi[q1].find(p1) == seq1_out_Multi[q1].end()) {
                                //     q1 = seq1->next_pair[nucp1][q1];
                                //     continue; // TODO under if-expression not while-expression
                                // }
                                int q1_j1 = q1 - j1;

                                // the current shape is p..i M2 j ..q
                                int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                                int min_p2 = max(hmmalign.low_bounds[p1], max(1, i2 - SINGLE_MAX_LEN));

                                // max_p2 = min(max_p2, i2-i1_p1+10); // lengths of aligned internal loops are close. p2 = i2+i1-p1 +- 10
                                // min_p2 = max(min_p2, i2-i1_p1-10); // lengths of aligned internal loops are close 

                                for (int p2 = max_p2; p2 >= min_p2; --p2) {
                                    int nucp2 = seq2->nucs[p2];
                                    int q2 = seq2->next_pair[nucp2][j2];

                                    // speed up
                                    int i2_p2 = i2 - p2;
                                    if (q2 > hmmalign.up_bounds[q1] || q2 == -1 || (i2_p2 + (q2 - j2) - 2 > SINGLE_MAX_LEN)) continue;
                                    // if (q2 > hmmalign.up_bounds[q1] || q2 == -1 || (q2 - j2 - 1 > 5)) continue;
                                    if (q2 < hmmalign.low_bounds[q1]) 
                                        q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1] - 1];

                                    if (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN) && seq2_out_Multi[q2].find(p2) != seq2_out_Multi[q2].end()) {
                                    // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                    // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= 5)) {
                                        // single seq folding
                                        // if (seq2_out_Multi[q2].find(p2) == seq2_out_Multi[q2].end()) {
                                        //     q2 = seq2->next_pair[nucp2][q2];
                                        //     continue;
                                        // }

                                        int newscore1 = multi_unpaired_score2(i1, j1, p1, q1, seq1);
                                        int newscore2 = multi_unpaired_score2(i2, j2, p2, q2, seq2);

                                        if (verbose) cout << "M2 to Multi: " << i1 << " " << j1 << " " << i2 << " " << j2  << " " << p1 << " " << q1 << " " << p2 << " " << q2 << endl;
                                        // cout << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << " " << bestMulti[q1+q2][q1][make_pair(p1, p2)].seq1foldscore << " " << bestMulti[q1+q2][q1][make_pair(p1, p2)].seq2foldscore << " " << bestMulti[q1+q2][q1][make_pair(p1, p2)].alignscore << endl;

                                        float pre_align_trans, post_align_trans, alignscore;
                                        // update bestMultiALN
                                        {
#ifdef dynalign
                                            aln_value_type alignscore = ALN_VALUE_MIN;
                                            if (p2>=hmmalign.low_bounds[p1] && p2<=hmmalign.up_bounds[p1]){ // only apply alignment constraint on left side
                                                alignscore = state.alignscore + abs(i1 - p1 - i2 + p2) + abs(q1 - j1 - q2 + j2);
                                            }
#else
                                            pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, state.startHMMstate-1);
                                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2, state.endHMMstate-1, 2, true); // true

                                            alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans);
#endif 
                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestMulti[q1+q2][q1][2][make_pair(p1, p2)], // bestMulti[q1+q2][get_keys(q1, p1, p2)].alnobj,
                                                                state.seq1foldscore + newscore1,
                                                                state.seq2foldscore + newscore2, 
                                                                state.manner, MANNER_MULTI_ALN,
                                                                static_cast<char>(i1 - p1), q1 - j1,
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                        }
#ifndef dynalign
                                        // update bestMultiINS1
                                        // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                        // if (p2+1-hmmalign.low_bounds[p1] >= 0)
                                        {
                                            pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, state.startHMMstate-1);
                                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, state.endHMMstate-1, 0, true); // true
                                            
                                            alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans);

                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestMulti[q1+q2][q1][0][make_pair(p1, p2)], // bestMulti[q1+q2][get_keys(q1, p1, p2)].ins1obj,
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
                                            pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, state.startHMMstate-1);
                                            post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, state.endHMMstate-1, 1, true); // true
                                            
                                            alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans);

                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestMulti[q1+q2][q1][1][make_pair(p1, p2)], // bestMulti[q1+q2][get_keys(q1, p1, p2)].ins2obj,
                                                                state.seq1foldscore + newscore1,
                                                                state.seq2foldscore + newscore2, 
                                                                state.manner, MANNER_MULTI_INS2,
                                                                static_cast<char>(i1 - p1), q1 - j1,
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                        }
#endif
                                        // q2 = seq2->next_pair[nucp2][q2];
                                    } // while loop enumerate q2
                                } // for loop enumerate p2
                                // q1 = seq1->next_pair[nucp1][q1];
                            } // q1
                        } // p1
                    } // 1. multi-loop
                }
            }// end beamM2
        }// beam of M2
     
        // beam of M
        // for every state in M[j]
        //   1. M = M + unpaired
        if (verbose) cout << "beam of M" << endl;
        {
            if (beam > 0) threshold = beam_prune(beamM, s, seq1_out_M, seq2_out_M);
            // if (beam > 0 && (s > 291)) cout << "s: " << s << " threshold: " << threshold << endl;

// #ifdef is_cube_pruning
//             sortM(threshold, beamM, sorted_bestM[s], s, seq1_out_M2, seq2_out_M2); // sorted_bestM include outside score
// #endif

            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (int m=0; m<3; m++) {
                    for (auto &item : beamM[j1][m]) {
                        State state = item.second;

                        int i1 = state.i1;
                        int i2 = state.i2;
                        assert (i1 > 0);
                        assert (i2 > 0);

                        if (verbose) cout << "M = M + U " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.seq1foldscore << " " << state.seq2foldscore << " " <<state.alignscore<< endl;

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
                                        update_if_better(i1, j1+1, i2, j2+1, bestM[s+2][j1+1][2][make_pair(i1, i2)], // bestM[s+2][get_keys(j1+1, i1, i2)]
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
                                        update_if_better(i1, j1+1, i2, j2, bestM[s+1][j1+1][0][make_pair(i1, i2)], // bestM[s+1][get_keys(j1+1, i1, i2)]
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
                                        update_if_better(i1, j1, i2, j2+1, bestM[s+1][j1][1][make_pair(i1, i2)], // bestM[s+1][get_keys(j1, i1, i2)]
                                                        state.seq1foldscore,
                                                        state.seq2foldscore,
                                                        state.manner, MANNER_M_eq_M_plus_U_INS2,
                                                        alignscore, state.startHMMstate, MANNER_INS2, weight, verbose);
                                }
                            }
                        }
                    } // for loop j1
                }
            }
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
                    if (state.manner == MANNER_NONE) continue;
                    // assert (state.endHMMstate == pre_manner);

                    if (verbose) cout << "C+U: " << j1 << " "  << j2  << " "<<  j1 << " "  << j2  << " " << pre_manner << endl;

                    float trans_emit_prob;
                    if (j1 == seq1->seq_len - 1 && j2 == seq2->seq_len - 1) {
                        cout << j1 << " "  << j2  << " " << state.score << endl;

                        trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_ALN, j1+1, j2+1, true);
                        float alignscore = xlog_mul(state.alignscore, trans_emit_prob);

                        if (alignscore > LOG_OF_ZERO)
                            update_if_better(0, j1+1, 0, j2+1, bestC[s+2][j1+1].alnobj,
                                            state.seq1foldscore,
                                            state.seq2foldscore,
                                            state.manner, MANNER_C_eq_C_plus_U_ALN,
                                            alignscore, state.startHMMstate, MANNER_ALN, weight, verbose);
                        
                        continue;
                    }

                    if (j1 < seq1->seq_len - 1 && j2 < seq2->seq_len - 1){ // ALN
                        // hmm constraint
                        if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]) {

#ifdef dynalign
                            int alignscore = state.alignscore; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_ALN, j1+1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 

                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1+1, 0, j2+1, bestC[s+2][j1+1].alnobj,
                                                state.seq1foldscore + newscore1,
                                                state.seq2foldscore + newscore2,
                                                state.manner, MANNER_C_eq_C_plus_U_ALN,
                                                alignscore, state.startHMMstate, MANNER_ALN, weight, verbose);
                        }
                    } 
                    if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) { // INS1
                        // hmm constraint
                        if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) {

#ifdef dynalign
                            int alignscore = state.alignscore + 1; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_INS1, j1+1, j2, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1+1, 0, j2, bestC[s+1][j1+1].ins1obj,
                                                state.seq1foldscore + newscore1,
                                                state.seq2foldscore,
                                                state.manner, MANNER_C_eq_C_plus_U_INS1,
                                                alignscore, state.startHMMstate, MANNER_INS1, weight, verbose);
                        }
                    } 
                    if (j2 < seq2->seq_len - 1 && j1 <= seq1->seq_len - 1) { // INS2
                        // hmm constraint
                        if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) {

#ifdef dynalign
                            int alignscore = state.alignscore + 1; 
#else
                            trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_INS2, j1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
#endif 

                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1, 0, j2+1, bestC[s+1][j1].ins2obj,
                                                state.seq1foldscore,
                                                state.seq2foldscore + newscore2,
                                                state.manner, MANNER_C_eq_C_plus_U_INS2,
                                                alignscore, state.startHMMstate, MANNER_INS2, weight, verbose);
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

    mem = GetProcessMemory();
    cout << "VmPeak: " << mem.VmPeak / 1024.0 / 1024.0 << endl;

    auto &state  = bestC[seq1->seq_len + seq2->seq_len][seq1->seq_len].alnobj; // bestC[seq1->seq_len - 1][seq2->seq_len-1];
    float bestscore = state.score;
    cout << "inside: " << state.score << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << xlog_mul(state.alignscore, hmmalign.trans_probs[2][2]) << endl;

    // backtrace
    // verbose = true;
    tuple<string, string, string, string> ret = get_parentheses(*seq1, *seq2, hmmalign);
    cout << get<0>(ret) <<endl;
    cout << get<1>(ret) <<endl;
    cout << get<2>(ret) <<endl;
    cout << get<3>(ret) <<endl;

    // release space
    for (int s = 1; s <= sum_len; ++s) {
        int minj1 = hmmalign.min_j1[s];
        int maxj1 = hmmalign.max_j1[s];
        if (maxj1 == -1) continue; // j1+j1==s not exist

        for (int j=minj1; j<maxj1; j++) {
            delete[] bestH[s][j];
            delete[] bestP[s][j];
            delete[] bestMulti[s][j];
        }

        bestH[s] = bestH[s] + minj1;
        bestP[s] = bestP[s] + minj1;
        bestMulti[s] = bestMulti[s] + minj1;
        bestM[s] = bestM[s] + minj1;
        bestM2[s] = bestM2[s] + minj1;

        delete[] bestH[s];
        delete[] bestP[s];
        delete[] bestMulti[s];
        delete[] bestM[s];
        delete[] bestM2[s];
    }
    delete[] bestH;
    delete[] bestP;
    delete[] bestMulti;
    delete[] bestM;
    delete[] bestM2;

// #ifndef multilign
//     return;
// #else
//     // outside
//     cout << "start outside computation" << endl;
//     // verbose = false;
//     outside(limited, allowed_pairs);
//     state  = bestC_beta[0][0]; // bestC[seq1->seq_len - 1][seq2->seq_len-1];
//     cout << "outside: " << state.score << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << endl;

//     // list of base pairs
//     // float best_score = bestC[seq1->seq_len - 1 + seq2->seq_len - 1][seq1->seq_len-1].score;
//     float min_score = bestscore * 1.2;
//     cout << "min_score: " << min_score << endl;
//     vector<tuple<float, int, int>> out_pairs_candidates;
//     for(int s = 1; s < sum_len - 1; ++s) {
//         unordered_map<int, State3>& beamP = bestP[s];
//         unordered_map<int, State3>& beamP_beta = bestP_beta[s];
//         for (auto &item : beamP) {
//             // for (int m=0; m<3; m++){
//             //     State *stateptr, *state_betaptr;
//             //     switch (m)
//             //     {
//             //         case 0:
//             //             stateptr = &item.second.ins1obj;
//             //             state_betaptr = &beamP_beta[item.first].ins1obj;
//             //             break;
//             //         case 1:
//             //             stateptr = &item.second.ins2obj;
//             //             state_betaptr = &beamP_beta[item.first].ins2obj;
//             //             break;
//             //         case 2:
//             //             stateptr = &item.second.alnobj;
//             //             state_betaptr = &beamP_beta[item.first].alnobj;
//             //             break;
                    
//             //         default:
//             //             break;
//             //     }

//             State& state = item.second.alnobj;
//             State& state_beta = beamP_beta[item.first].alnobj;
//             if (state.i1 <= 0) continue;

//             if (state.score == VALUE_FMIN || state_beta.score == VALUE_FMIN) continue;

//             float score = state.score + state_beta.score;

//             if (min_score > 0) {
//                 if (score > min_score) continue;
//             } else {
//                 if (score < min_score) continue;
//             }

//             int i1 = state.i1;
//             int j1 = state.j1;
//             // cout << i1 << " " << j1 << endl;
//             if (limited) {
//                 if (allowed_pairs.find(make_pair(i1, j1)) == allowed_pairs.end()) 
//                     continue;
//             }

//             out_pairs_candidates.push_back(make_tuple(score, state.i1, state.j1));

//             // cout << state.i1 << " " << state.j1 << " " << state.i2 << " " << s - state.j1 <<  " ";
//             // cout << state.score  << " " << state_beta.score << " " << state.score + state_beta.score << endl;
//             // cout << state.seq1foldscore + state_beta.seq1foldscore << " " <<  state.seq2foldscore + state_beta.seq2foldscore << " ";
//             // cout << state.seq1foldscore + state_beta.seq1foldscore + state.seq2foldscore + state_beta.seq2foldscore << " ";
//             // cout << state.alignscore  << " " << state_beta.alignscore << " " << state.alignscore + state_beta.alignscore << endl;
//         }
//     }
//     cout << "out_pairs_candidates size: " << out_pairs_candidates.size() << " " << num_pairs << endl;

//     // only save top num_pairs base pairs using quickselect
//     int current_num = out_pairs_candidates.size();
//     if (current_num > num_pairs) {
//         float threshold = quickselect(out_pairs_candidates, 0, current_num-1, current_num - num_pairs);
//         float threshold2 = bestscore * 1.01;
//         cout << "threshold: " << threshold  << " threshold2: " << threshold2 << endl;
//         threshold = min(threshold, threshold2);
//         cout << "threshold: " << threshold  << endl;
//         for (auto &p : out_pairs_candidates) {
//             if (get<0>(p) > threshold)
//                 out_pairs.push_back(make_pair(get<1>(p), get<2>(p)));
//         }
//     } else {
//         for (auto &p : out_pairs_candidates) {
//             out_pairs.push_back(make_pair(get<1>(p), get<2>(p)));
//         }
//     }
//     cout << "out_pairs size: " << out_pairs.size() << " " << num_pairs << endl;
// #endif

    // clear allocated space in alignment
    // cout << "clear allocated space in alignment" << endl;
#ifdef dynalign
    hmmalign.clear(false);
#else
    hmmalign.clear(true);
#endif
}

BeamSankoffParser::BeamSankoffParser(float aln_weight, int beam_size, int LFbeam, int LAbeam, bool if_aster, bool if_suffix, float energy_diff, bool is_verbose)
    :weight(aln_weight),
     beam(beam_size),
     lfbeam(LFbeam),
     alnbeam(LAbeam),
     use_astar(if_aster),
     use_suffix(if_suffix),
     max_energy_diff(energy_diff),
     verbose(is_verbose){

    cout << "beam : " << beam << " lfbeam: " << lfbeam << " alnbeam: " << alnbeam << " use_astar: " << use_astar << " use_suffix: " << use_suffix <<  " max_energy_diff: " << max_energy_diff << endl;  
    
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