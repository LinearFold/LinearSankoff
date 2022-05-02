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

#include "LinearSankoff_array.h"
// #include "LinearSankoff.cpp"
// #include "HMMAlign.h"

// #include <sparsehash/dense_hash_map>
// using google::dense_hash_map;      // namespace where class lives by default

// #define dynalign

using namespace std;

tuple<string, string, string, string> SankoffParser::get_parentheses_H(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
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

tuple<string, string, string, string> SankoffParser::get_parentheses_Multi(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
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
        stk.push(make_tuple(p1, q1, p2, q2, bestM2[p1][q1][p2][q2]));
        tuple<string, string, string, string> m2result = get_parentheses_M2(seq1, seq2, stk, hmmalign);

        pair<string, string> pre_aln, post_aln;
        string seq1_struc, seq2_struc, seq1_aln, seq2_aln;
        switch(state.manner) {
            case MANNER_MULTI_eq_MULTI_plus_U_ALN: case MANNER_MULTI_ALN:
                pre_aln = get_hmm_aln_left(i1, p1, i2, p2, MANNER_ALN, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1, q2, j2, MANNER_ALN, MANNER_ALN);

                seq1_struc = struc_p2p_pair(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_pair(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            case MANNER_MULTI_eq_MULTI_plus_U_INS1: case MANNER_MULTI_INS1:
                pre_aln = get_hmm_aln_left(i1, p1, i2+1, p2, MANNER_INS1, MANNER_ALN);
                post_aln = get_hmm_aln_right(q1, j1, q2, j2-1, MANNER_ALN, MANNER_INS1);

                seq1_struc = struc_p2p_pair(get<0>(m2result), seq1_l1-1, seq1_l2-1);
                seq2_struc = struc_p2p_dots(get<1>(m2result), seq2_l1-1, seq2_l2-1);

                break;
            case MANNER_MULTI_eq_MULTI_plus_U_INS2: case MANNER_MULTI_INS2:
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

tuple<string, string, string, string> SankoffParser::get_parentheses_M2(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
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

        Manner pre_manner = state.premanner;
        switch (pre_manner)
        {
            case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_M2: case MANNER_M_eq_P:
                stk.push(make_tuple(i1, k1, i2, k2, bestM[i1][k1][i2][k2].alnobj)); // push M1
                break;
            case MANNER_M_eq_M_plus_U_INS1:
                stk.push(make_tuple(i1, k1, i2, k2, bestM[i1][k1][i2][k2].ins1obj)); // push M1
                break;
            case MANNER_M_eq_M_plus_U_INS2:
                stk.push(make_tuple(i1, k1, i2, k2, bestM[i1][k1][i2][k2].ins2obj)); // push M1
                break;
            
            default:
                break;
        }
        // stk.push(make_tuple(i1, k1, i2, k2, bestM[i1][k1][i2][k2])); // push M1
        tuple<string, string, string, string> pre_result = get_parentheses_M1(seq1, seq2, stk, hmmalign);

        stack<tuple<int, int, int, int, State>> stk2; // push P
        // stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[j1+j2][get_keys(j1, k1+1, k2+1)].alnobj)); 
        stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[k1+1][j1][k2+1][j2].alnobj)); 
        tuple<string, string, string, string> result = get_parentheses_P(seq1, seq2, stk2, hmmalign);

        return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
    }
    
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> SankoffParser::get_parentheses_M1(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
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
                Manner pre_manner = state.premanner;
                switch (pre_manner)
                {
                    case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_M2: case MANNER_M_eq_P:
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[i1][j1-1][i2][j2-1].alnobj)); // push M1
                        break;
                    case MANNER_M_eq_M_plus_U_INS1:
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[i1][j1-1][i2][j2-1].ins1obj)); // push M1
                        break;
                    case MANNER_M_eq_M_plus_U_INS2:
                        stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[i1][j1-1][i2][j2-1].ins2obj)); // push M1
                        break;
                    
                    default:
                        break;
                }
                // stk.push(make_tuple(i1, j1-1, i2, j2-1, bestM[i1][j1-1][i2][j2-1]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result)+".", get<1>(result)+".", get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+seq2.raw_seq.at(j2));
                break;
            }
            case MANNER_M_eq_M_plus_U_INS1:
            {
                if (verbose) cout << "MANNER_M_eq_M_plus_U_INS1" << endl;
                Manner pre_manner = state.premanner;
                switch (pre_manner)
                {
                    case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_M2: case MANNER_M_eq_P:
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[i1][j1-1][i2][j2].alnobj)); // push M1
                        break;
                    case MANNER_M_eq_M_plus_U_INS1:
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[i1][j1-1][i2][j2].ins1obj)); // push M1
                        break;
                    case MANNER_M_eq_M_plus_U_INS2:
                        stk.push(make_tuple(i1, j1-1, i2, j2, bestM[i1][j1-1][i2][j2].ins2obj)); // push M1
                        break;
                    
                    default:
                        break;
                }
                // stk.push(make_tuple(i1, j1-1, i2, j2, bestM[i1][j1-1][i2][j2])); 
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result)+".", get<1>(result), get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+"-");
                break;
            }
            case MANNER_M_eq_M_plus_U_INS2:
            {
                if (verbose) cout << "MANNER_M_eq_M_plus_U_INS2" << endl;
                Manner pre_manner = state.premanner;
                switch (pre_manner)
                {
                    case MANNER_M_eq_M_plus_U_ALN: case MANNER_M_eq_M2: case MANNER_M_eq_P:
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[i1][j1][i2][j2-1].alnobj)); // push M1
                        break;
                    case MANNER_M_eq_M_plus_U_INS1:
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[i1][j1][i2][j2-1].ins1obj)); // push M1
                        break;
                    case MANNER_M_eq_M_plus_U_INS2:
                        stk.push(make_tuple(i1, j1, i2, j2-1, bestM[i1][j1][i2][j2-1].ins2obj)); // push M1
                        break;
                    
                    default:
                        break;
                }
                // stk.push(make_tuple(i1, j1, i2, j2-1, bestM[i1][j1][i2][j2-1]));
                result = get_parentheses_M1(seq1, seq2, stk, hmmalign);
                return make_tuple(get<0>(result), get<1>(result)+".", get<2>(result)+"-", get<3>(result)+seq2.raw_seq.at(j2));
                break;
            }
            case MANNER_M_eq_M2:
                if (verbose) cout << "MANNER_M_eq_M2" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestM2[j1+j2][get_keys(j1, i1, i2)])); 
                stk.push(make_tuple(i1, j1, i2, j2, bestM2[i1][j1][i2][j2])); 
                return get_parentheses_M2(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_M_eq_P:
                if (verbose) cout << "MANNER_M_eq_P" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestP[j1+j2][get_keys(j1, i1, i2)].alnobj));
                stk.push(make_tuple(i1, j1, i2, j2, bestP[i1][j1][i2][j2].alnobj)); 
                return get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;  
            default:
                printf("wrong manner at %d, %d, %d, %d: manner %d %d\n", i1, j1, i2, j2, state.premanner, state.manner); fflush(stdout);
                assert(false);
        }
    }
    
    return make_tuple("", "", "", "");
}

tuple<string, string, string, string> SankoffParser::get_parentheses_P(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
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
                // stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][get_keys(j1, i1, i2)].alnobj)); 
                stk.push(make_tuple(i1, j1, i2, j2, bestH[i1][j1][i2][j2].alnobj)); 
                
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_HAIRPIN_INS1:
                if (verbose) cout << "MANNER_HAIRPIN_INS1" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][get_keys(j1, i1, i2)].ins1obj)); 
                stk.push(make_tuple(i1, j1, i2, j2, bestH[i1][j1][i2][j2].ins1obj)); 
                
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_HAIRPIN_INS2:
                if (verbose) cout << "MANNER_HAIRPIN_INS2" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestH[j1+j2][get_keys(j1, i1, i2)].ins2obj)); 
                stk.push(make_tuple(i1, j1, i2, j2, bestH[i1][j1][i2][j2].ins2obj)); 
                
                return get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_ALN:
                if (verbose) cout << "MANNER_P_eq_MULTI_ALN" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][get_keys(j1, i1, i2)].alnobj)); 
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[i1][j1][i2][j2].alnobj)); 
                
                return get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_INS1:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS1" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][get_keys(j1, i1, i2)].ins1obj));
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[i1][j1][i2][j2].ins1obj)); 
                
                return get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;
            case MANNER_P_eq_MULTI_INS2:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS2" << endl;
                // stk.push(make_tuple(i1, j1, i2, j2, bestMulti[j1+j2][get_keys(j1, i1, i2)].ins2obj)); 
                stk.push(make_tuple(i1, j1, i2, j2, bestMulti[i1][j1][i2][j2].ins2obj)); 
                
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

        if (verbose) printf("get_parentheses_P manner at %d, %d, %d, %d: manner %d %d\n", p1, q1, p2, q2, state.premanner, state.manner);

        // get inner result 
        tuple<string, string, string, string> result;
        switch(state.premanner)
        {
            case MANNER_SINGLE_ALN: 
                if (verbose) cout << "MANNER_SINGLE_ALN" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].alnobj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestP[p1][q1][p2][q2].alnobj)); 
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_SINGLE_INS1:
                if (verbose) cout << "MANNER_SINGLE_INS1" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins1obj));
                stk.push(make_tuple(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins1obj));
                result = get_parentheses_P(seq1, seq2, stk, hmmalign); 
                break;

            case MANNER_SINGLE_INS2:
                if (verbose) cout << "MANNER_SINGLE_INS2" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestP[q1+q2][get_keys(q1, p1, p2)].ins2obj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins2obj)); 
                result = get_parentheses_P(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_ALN: 
                if (verbose) cout << "MANNER_HAIRPIN_ALN" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][get_keys(q1, p1, p2)].alnobj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestH[p1][q1][p2][q2].alnobj)); 
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_INS1:
                if (verbose) cout << "MANNER_HAIRPIN_INS1" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][get_keys(q1, p1, p2)].ins1obj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestH[p1][q1][p2][q2].ins1obj)); 
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_HAIRPIN_INS2:
                if (verbose) cout << "MANNER_HAIRPIN_INS2" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestH[q1+q2][get_keys(q1, p1, p2)].ins2obj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestH[p1][q1][p2][q2].ins2obj)); 
                result = get_parentheses_H(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_ALN: 
                if (verbose) cout << "MANNER_P_eq_MULTI_ALN" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].alnobj));  
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[p1][q1][p2][q2].alnobj));  
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_INS1:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS1" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].ins1obj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[p1][q1][p2][q2].ins1obj)); 
                result = get_parentheses_Multi(seq1, seq2, stk, hmmalign);
                break;

            case MANNER_P_eq_MULTI_INS2:
                if (verbose) cout << "MANNER_P_eq_MULTI_INS2" << endl;
                // stk.push(make_tuple(p1, q1, p2, q2, bestMulti[q1+q2][get_keys(q1, p1, p2)].ins2obj)); 
                stk.push(make_tuple(p1, q1, p2, q2, bestMulti[p1][q1][p2][q2].ins2obj)); 
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

tuple<string, string, string, string> SankoffParser::get_parentheses_C(SeqObject& seq1, SeqObject& seq2, stack<tuple<int, int, int, int, State>> &stk, BeamAlign &hmmalign) {
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

                    // stk.push(make_tuple(i1, k1, i2, k2, bestC[k1+k2][k1])); // push C
                    Manner pre_manner = state.premanner;
                    switch (pre_manner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1][k2].alnobj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1][k2].ins1obj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, k1, i2, k2, bestC[k1][k2].ins2obj)); // push C
                            break;
                        
                        default:
                            break;
                    }
                    tuple<string, string, string, string> pre_result = get_parentheses_C(seq1, seq2, stk, hmmalign);

                    stack<tuple<int, int, int, int, State>> stk2; // push P
                    stk2.push(make_tuple(k1+1, j1, k2+1, j2, bestP[k1+1][j1][k2+1][j2].alnobj)); 
                    tuple<string, string, string, string> result = get_parentheses_P(seq1, seq2, stk2, hmmalign);

                    return make_tuple(get<0>(pre_result)+get<0>(result), get<1>(pre_result)+get<1>(result), get<2>(pre_result)+get<2>(result), get<3>(pre_result)+get<3>(result));
                    break;
                }
                
            case MANNER_C_eq_P:
                {
                    if (verbose) cout <<  "MANNER_C_eq_P" << endl;
                    stk.push(make_tuple(i1+1, j1, i2+1, j2, bestP[i1+1][j1][i2+1][j2].alnobj)); // N.B.
                    
                    return get_parentheses_P(seq1, seq2, stk, hmmalign);
                    break;
                }
            case MANNER_C_eq_C_plus_U_ALN:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_ALN" << endl;
                    // stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1][j2-1].alnobj)); // push C
                    Manner pre_manner = state.premanner;
                    switch (pre_manner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1][j2-1].alnobj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1][j2-1].ins1obj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, j1-1, i2, j2-1, bestC[j1-1][j2-1].ins2obj)); // push C
                            break;
                        
                        default:
                            break;
                    }
                    tuple<string, string, string, string> result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    return make_tuple(get<0>(result)+".", get<1>(result)+".", get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+seq2.raw_seq.at(j2));
                    break;
                }
            case MANNER_C_eq_C_plus_U_INS1:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_INS1" << endl;
                    // stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1][j2].ins1obj)); // push C
                    Manner pre_manner = state.premanner;
                    switch (pre_manner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1][j2].alnobj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1][j2].ins1obj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, j1-1, i2, j2, bestC[j1-1][j2].ins2obj)); // push C
                            break;
                        
                        default:
                            break;
                    }
                    tuple<string, string, string, string> result = get_parentheses_C(seq1, seq2, stk, hmmalign);
                    return make_tuple(get<0>(result)+".", get<1>(result), get<2>(result)+seq1.raw_seq.at(j1), get<3>(result)+"-");
                    break;
                }
            case MANNER_C_eq_C_plus_U_INS2:
                {
                    if (verbose) cout <<  "MANNER_C_eq_C_plus_U_INS2" << endl;
                    // stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1][j2-1].ins2obj)); // push C
                    Manner pre_manner = state.premanner;
                    switch (pre_manner)
                    {
                        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1][j2-1].alnobj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS1:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1][j2-1].ins1obj)); // push C
                            break;
                        case MANNER_C_eq_C_plus_U_INS2:
                            stk.push(make_tuple(i1, j1, i2, j2-1, bestC[j1][j2-1].ins2obj)); // push C
                            break;
                        
                        default:
                            break;
                    }
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

tuple<string, string, string, string> SankoffParser::get_parentheses(SeqObject& seq1, SeqObject& seq2, BeamAlign &hmmalign) {
    stack<tuple<int, int, int, int, State>> stk;
    
    Manner pre_manner = bestC[seq1_len][seq2_len].alnobj.premanner;
    switch (pre_manner)
    {
        case MANNER_C_eq_C_plus_U_ALN: case MANNER_C_eq_C_plus_P:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1_len-1][seq2_len-1].alnobj));
            break;
        case MANNER_C_eq_C_plus_U_INS1:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1_len-1][seq2_len-1].ins1obj));
            break;
        case MANNER_C_eq_C_plus_U_INS2:
            stk.push(make_tuple(0, seq1.seq_len-1, 0, seq2.seq_len-1, bestC[seq1_len-1][seq2_len-1].ins2obj));
            break;
        
        default:
            break;
    }
    return get_parentheses_C(seq1, seq2, stk, hmmalign);
}

float SankoffParser::get_hmm_score(int i1, int j1, int i2, int j2, int s1, bool allowout){
    // (...)
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return LOG_OF_ZERO;
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return LOG_OF_ZERO;

    // return hmmalign.all_local_scores[s1][i1][i2-hmmalign.low_bounds[i1]][j1-i1+1][j2-i2+1];
    return hmmalign.all_local_scores[s1][i1][i2][j1][j2];
}

float SankoffParser::get_hmm_score_left(int i1, int j1, int i2, int j2, int s1, int s2){
    // (-->(
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return LOG_OF_ZERO;
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return LOG_OF_ZERO;

    // return hmmalign.left_local_scores[s1][s2][i1][i2-hmmalign.low_bounds[i1]][j1-i1+1][j2-i2+1];
    return hmmalign.left_local_scores[s1][s2][i1][i2][j1][j2];
}

float SankoffParser::get_hmm_score_right(int i1, int j1, int i2, int j2, int s1, int s2, bool allowout){
    // )-->)
    if (i2 < hmmalign.low_bounds[i1] || i2 > hmmalign.up_bounds[i1]) return LOG_OF_ZERO;
    if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) return LOG_OF_ZERO;

    if ((j1-i1+1 > 30) || (j2-i2+1 > 30))
        return hmmalign.viterbi_path_local_right(i1, j1, i2, j2, static_cast<HMMManner>(s1+1), static_cast<HMMManner>(s2+1));

    // return hmmalign.right_local_scores[s1][s2][i1][i2-hmmalign.low_bounds[i1]][j1-i1+1][j2-i2+1];
    return hmmalign.right_local_scores[s1][s2][i1][i2][j1][j2];
}

pair<string, string> SankoffParser::get_hmm_aln(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2){
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

pair<string, string> SankoffParser::get_hmm_aln_left(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2){
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

pair<string, string> SankoffParser::get_hmm_aln_right(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2){
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

/*
void SankoffParser::outside(bool limited, const set<pair<int, int>> &allowed_pairs){
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

void SankoffParser::prepare(const vector<string> &seqs){
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
    cout << "sum of seq length: " << sum_len << endl;

    // HMM align tool
    hmmalign.set(alnbeam, sequences[0].nucs, sequences[1].nucs);
    hmmalign.viterbi_path(false); // alignment envelope
    hmmalign.viterbi_path_all_locals(); // compute and save local alignments

    // single sequence folding
    BeamCKYParser* cky_parser = new BeamCKYParser();
    cky_parser->beam = lfbeam;

    // seq1_H_pairs = new bool*[seq1_len+1];
    // seq1_P_pairs = new bool*[seq1_len+1];
    // seq1_Multi_pairs = new bool*[seq1_len+1];
    // seq1_M_pairs = new bool*[seq1_len+1];
    // seq1_M2_pairs = new bool*[seq1_len+1];
    // for (int j=1; j<=seq1_len; j++) {
    //     seq1_H_pairs[j] = new bool[j];
    //     seq1_P_pairs[j] = new bool[j];
    //     seq1_Multi_pairs[j] = new bool[j];
    //     seq1_M_pairs[j] = new bool[j];
    //     seq1_M2_pairs[j] = new bool[j];
    //     for (int i=0; i<j; i++) {
    //         seq1_H_pairs[j][i] = false;
    //         seq1_P_pairs[j][i] = false;
    //         seq1_Multi_pairs[j][i] = false;
    //         seq1_M_pairs[j][i] = false;
    //         seq1_M2_pairs[j][i] = false;
    //     }
    // }
    // seq2_H_pairs = new bool*[seq2_len+1];
    // seq2_P_pairs = new bool*[seq2_len+1];
    // seq2_Multi_pairs = new bool*[seq2_len+1];
    // seq2_M_pairs = new bool*[seq2_len+1];
    // seq2_M2_pairs = new bool*[seq2_len+1];
    // for (int j=1; j<=seq2_len; j++) {
    //     seq2_H_pairs[j] = new bool[j];
    //     seq2_P_pairs[j] = new bool[j];
    //     seq2_Multi_pairs[j] = new bool[j];
    //     seq2_M_pairs[j] = new bool[j];
    //     seq2_M2_pairs[j] = new bool[j];
    //     for (int i=0; i<j; i++) {
    //         seq2_H_pairs[j][i] = false;
    //         seq2_P_pairs[j][i] = false;
    //         seq2_Multi_pairs[j][i] = false;
    //         seq2_M_pairs[j][i] = false;
    //         seq2_M2_pairs[j][i] = false;
    //     }
    // }


    seq1_H_pairs.clear();
    seq1_H_pairs.resize(seq1_len);
    seq1_P_pairs.clear();
    seq1_P_pairs.resize(seq1_len);
    seq1_M_pairs.clear();
    seq1_M_pairs.resize(seq1_len);
    seq1_M2_pairs.clear();
    seq1_M2_pairs.resize(seq1_len);
    seq1_Multi_pairs.clear();
    seq1_Multi_pairs.resize(seq1_len);

    seq2_H_pairs.clear();
    seq2_H_pairs.resize(seq2_len);
    seq2_P_pairs.clear();
    seq2_P_pairs.resize(seq2_len);
    seq2_M_pairs.clear();
    seq2_M_pairs.resize(seq2_len);
    seq2_M2_pairs.clear();
    seq2_M2_pairs.resize(seq2_len);
    seq2_Multi_pairs.clear();
    seq2_Multi_pairs.resize(seq2_len);

    // internal loop score
    seq1_internal = new short***[seq1_len];
    seq2_internal = new short***[seq2_len];
    for (int i1=2; i1<seq1_len-3; i1++) {
        seq1_internal[i1] = new short**[seq1_len - i1 - 3];
        seq1_internal[i1] = seq1_internal[i1] - (i1+4);
        int nuci = sequences[0].nucs[i1];

        short mini = max(1, i1-MAX_LOOP_LEN-1);
        for (int j1 = i1+4; j1<seq1_len-1; j1++) {
            int nucj = sequences[0].nucs[j1];
            if (_allowed_pairs[nuci][nucj]) {
                seq1_internal[i1][j1] = new short*[i1-mini+1];
                seq1_internal[i1][j1] = seq1_internal[i1][j1] - mini;

                short minj = j1+1;
                short maxj = min(minj+MAX_LOOP_LEN+1, seq1_len);
                for (int k=mini; k<i1; k++) {
                    seq1_internal[i1][j1][k] = new short[maxj-minj+1];
                    seq1_internal[i1][j1][k] = seq1_internal[i1][j1][k] - minj; 

                    // for (int l=minj; l<maxj; l++) {
                    //     seq1_internal[i1][j1][k][l] = -1;
                    // }

                }
            }            
        }
    }
    for (int i2=2; i2<seq2_len-3; i2++) {
        seq2_internal[i2] = new short**[seq2_len - i2 - 3];
        seq2_internal[i2] = seq2_internal[i2] - (i2+4); 
        int nuci = sequences[1].nucs[i2];

        short mini = max(1, i2-MAX_LOOP_LEN-1);
        for (int j2 = i2+4; j2<seq2_len-1; j2++) {
            int nucj = sequences[1].nucs[j2];
            if (_allowed_pairs[nuci][nucj]) {
                seq2_internal[i2][j2] = new short*[i2-mini+1];
                seq2_internal[i2][j2] = seq2_internal[i2][j2] - mini;
                
                short minj = j2+1;
                short maxj = min(minj+MAX_LOOP_LEN+1, seq2_len);
                for (int k=mini; k<i2; k++) {
                    seq2_internal[i2][j2][k] = new short[maxj-minj+1];
                    seq2_internal[i2][j2][k] = seq2_internal[i2][j2][k] - minj;

                    // for (int l=minj; l<maxj; l++) {
                    //     seq2_internal[i2][j2][k][l] = -1;
                    // }
                }
            }
        }
    }

    // int max_range = -1;
    for (int i_seq=0; i_seq<2; i_seq++) {
        // single sequence folding
        int seq_viterbi;
        if (i_seq == 0)
            seq_viterbi = cky_parser->parse(seqs[i_seq], seq1_internal, NULL);
        else
            seq_viterbi = cky_parser->parse(seqs[i_seq], seq2_internal, NULL);
        
        // the maximum % change in free energy from the lowest free energy structure
        float min_score = seq_viterbi * (1 - max_energy_diff); // max_score must be positive
        cout << "min_score: " << min_score << endl;
        // int min_score = seq_viterbi - 1000;
        
        // only keep suboptimal base pairs
        int seq_len = i_seq==0? seq1_len: seq2_len;
        for (int j=0; j<seq_len-1; j++) {
            vector<unordered_map<int, LFState>*> seq_in{&cky_parser->bestH[j], &cky_parser->bestP[j], &cky_parser->bestM[j], &cky_parser->bestM2[j], &cky_parser->bestMulti[j]};
            vector<unordered_map<int, LFState>*> seq_out{&cky_parser->bestH_beta[j], &cky_parser->bestP_beta[j], &cky_parser->bestM_beta[j], &cky_parser->bestM2_beta[j], &cky_parser->bestMulti_beta[j]};
            
            vector<unordered_map<int, BoolState>*> valid_pairs;
            if (i_seq==0)
                valid_pairs = {&seq1_H_pairs[j+1], &seq1_P_pairs[j+1], &seq1_M_pairs[j+1], &seq1_M2_pairs[j+1], &seq1_Multi_pairs[j+1]};
            else
                valid_pairs = {&seq2_H_pairs[j+1], &seq2_P_pairs[j+1], &seq2_M_pairs[j+1], &seq2_M2_pairs[j+1], &seq2_Multi_pairs[j+1]};

            for (int k=0; k < 5; k++){
                unordered_map<int, LFState>& beamins = *seq_in[k];
                unordered_map<int, LFState>& beamout = *seq_out[k];
                for(auto& item : beamins) {
                    int i = item.first;
                    int in_score = item.second.score;
                    int out_score = beamout[i].score;

                    if (in_score == VALUE_MIN || out_score == VALUE_MIN) {
                        // (*valid_pairs[k])[i+1] = false;
                        continue;
                    }
                        
                    if (in_score + out_score >= min_score) {
                        (*valid_pairs[k])[i+1].valid = true;
                    }
                }
            }            
        }
    }
    delete cky_parser;

    // allocate memory
    // allocate space for bestC
    bestC = new State3*[seq1_len + 1];
    for (int j1 = 0; j1 < seq1_len; j1++) {
        short low_bound = hmmalign.low_bounds[j1];
        short range = hmmalign.up_bounds[j1] - low_bound + 1;
        bestC[j1] = new State3[range];
        bestC[j1] = bestC[j1] - low_bound;

        for (int j2 = low_bound; j2 <=  hmmalign.up_bounds[j1]; j2++)
            bestC[j1][j2].init();
    }
    bestC[seq1_len] = new State3[1];
    bestC[seq1_len] = bestC[seq1_len] - seq2_len;
    bestC[seq1_len][seq2_len].init();

    // allocate space for bestH/P/Multi/M/M2/
    bestH = new State3***[seq1_len];
    bestP = new State3***[seq1_len];
    bestMulti = new State3***[seq1_len];
    bestM = new State3***[seq1_len];
    bestM2 = new State***[seq1_len];
    
    for (int i1=1; i1<seq1_len; i1++) {
        bestH[i1] = new State3**[seq1_len - i1 + 1];
        bestH[i1] = bestH[i1] - i1;

        bestP[i1] = new State3**[seq1_len - i1 + 1];
        bestP[i1] = bestP[i1] - i1;

        bestMulti[i1] = new State3**[seq1_len - i1 + 1];
        bestMulti[i1] = bestMulti[i1] - i1;

        bestM[i1] = new State3**[seq1_len - i1 + 1];
        bestM[i1] = bestM[i1] - i1;

        bestM2[i1] = new State**[seq1_len - i1 + 1];
        bestM2[i1] = bestM2[i1] - i1;

        int nuci1 = sequences[0].nucs[i1];
        for (int j1 = i1+1; j1<seq1_len; j1++) {
            short i1_lowbound = hmmalign.low_bounds[i1];
            short i1_upbound = hmmalign.up_bounds[i1];
            short i1_range = i1_upbound - i1_lowbound + 1;
            
            int nucj1 = sequences[0].nucs[j1];

            // if (verbose) cout << i1 << " " << j1 << " " << i1_lowbound << " " <<  i1_upbound << " " << i1_range << endl;

            if (_allowed_pairs[nuci1][nucj1]) {
                if (seq1_H_pairs[j1].find(i1) != seq1_H_pairs[j1].end()) {
                    bestH[i1][j1] = new State3*[i1_range];
                    bestH[i1][j1] = bestH[i1][j1] - i1_lowbound; // shift pointer
                }
                if (seq1_P_pairs[j1].find(i1) != seq1_P_pairs[j1].end()) {
                    bestP[i1][j1] = new State3*[i1_range];
                    bestP[i1][j1] = bestP[i1][j1] - i1_lowbound;
                }
                if (seq1_Multi_pairs[j1].find(i1) != seq1_Multi_pairs[j1].end()) {
                    bestMulti[i1][j1] = new State3*[i1_range];
                    bestMulti[i1][j1] = bestMulti[i1][j1] - i1_lowbound;
                }
            }
            if (seq1_M_pairs[j1].find(i1) != seq1_M_pairs[j1].end()) {
                bestM[i1][j1] = new State3*[i1_range];
                bestM[i1][j1] = bestM[i1][j1] - i1_lowbound;
            }
            if (seq1_M2_pairs[j1].find(i1) != seq1_M2_pairs[j1].end()) {
                bestM2[i1][j1] = new State*[i1_range];
                bestM2[i1][j1] = bestM2[i1][j1] - i1_lowbound;
            }
                
            for (int i2 = i1_lowbound; i2 <= i1_upbound; i2++) {
                short j1_lowbound = hmmalign.low_bounds[j1];
                short j1_upbound = hmmalign.up_bounds[j1];
                short j1_range = j1_upbound - j1_lowbound + 1;

                // cout << i1 << " " << j1  << " " << i2 << " " << j1_lowbound << " " <<  j1_upbound << " " << j1_range << endl;

                if (_allowed_pairs[nuci1][nucj1]) {
                    if (seq1_H_pairs[j1].find(i1) != seq1_H_pairs[j1].end()) {
                        bestH[i1][j1][i2] = new State3[j1_range];
                        bestH[i1][j1][i2] = bestH[i1][j1][i2] - j1_lowbound;
                        // initialize
                        for (int j2 = j1_lowbound; j2 <= j1_upbound; j2++) {
                            if (seq2_H_pairs[j2].find(i2) != seq2_H_pairs[j2].end()) {
                                bestH[i1][j1][i2][j2].init();
                            }
                        }
                    }
                    if (seq1_P_pairs[j1].find(i1) != seq1_P_pairs[j1].end()) {
                        bestP[i1][j1][i2] = new State3[j1_range];
                        bestP[i1][j1][i2] = bestP[i1][j1][i2] - j1_lowbound;
                        // initialize
                        for (int j2 = j1_lowbound; j2 <= j1_upbound; j2++) {
                            if (seq2_P_pairs[j2].find(i2) != seq2_P_pairs[j2].end()) {
                                bestP[i1][j1][i2][j2].init();
                            }
                        }
                    }
                    if (seq1_Multi_pairs[j1].find(i1) != seq1_Multi_pairs[j1].end()) {
                        bestMulti[i1][j1][i2] = new State3[j1_range];
                        bestMulti[i1][j1][i2] = bestMulti[i1][j1][i2] - j1_lowbound;
                        // initialize
                        for (int j2 = j1_lowbound; j2 <= j1_upbound; j2++) {
                            if (seq2_Multi_pairs[j2].find(i2) != seq2_Multi_pairs[j2].end()) {
                                bestMulti[i1][j1][i2][j2].init();
                            }
                        }
                    }
                }
                if (seq1_M_pairs[j1].find(i1) != seq1_M_pairs[j1].end()) {
                    bestM[i1][j1][i2] = new State3[j1_range];
                    bestM[i1][j1][i2] = bestM[i1][j1][i2] - j1_lowbound;
                    // initialize
                    for (int j2 = j1_lowbound; j2 <= j1_upbound; j2++) {
                        if (seq2_M_pairs[j2].find(i2) != seq2_M_pairs[j2].end()) {
                            bestM[i1][j1][i2][j2].init();
                        }
                    }
                }
                if (seq1_M2_pairs[j1].find(i1) != seq1_M2_pairs[j1].end()) {
                    bestM2[i1][j1][i2] = new State[j1_range];
                    bestM2[i1][j1][i2] = bestM2[i1][j1][i2] - j1_lowbound;
                    // initialize
                    for (int j2 = j1_lowbound; j2 <= j1_upbound; j2++) {
                        if (seq2_M2_pairs[j2].find(i2) != seq2_M2_pairs[j2].end()) {
                            bestM2[i1][j1][i2][j2].init();
                        }
                    }
                }
            }
        }
    } 
}


#ifdef multilign
void SankoffParser::parse(const vector<string> &seqs, bool limited, const set<pair<int, int>> &allowed_pairs, vector<pair<int, int>> &out_pairs, int num_pairs){
#else
void SankoffParser::parse(const vector<string> &seqs){
#endif
    // test google dense hash map 
    // dense_hash_map<pair<int, int>, int, pair_hash, eqstr> test;
    // pair<int, int> test_key;
    // test.set_empty_key(test_key);
    // test[make_pair(1, 1)] = 2;
    // cout << "test google dense hash map : " << test[make_pair(1, 1)] << endl;

    struct timeval parse_starttime, parse_endtime;
    double parse_elapsed_time;
    // gettimeofday(&parse_starttime, NULL);
    
    // allocate space and pre-computation
    seq1_len = seqs[0].size() + 1;
    seq2_len = seqs[1].size() + 1;
    cout << "sequence length in SankoffParser parse: " << seq1_len << " " << seq2_len << endl;
    prepare(seqs);

    processMem_t mem = GetProcessMemory();
    cout << "VmPeak: " << mem.VmPeak / 1024.0 / 1024.0 << endl;


    // gettimeofday(&parse_endtime, NULL);
    // double parse_elapsed_time;
    // parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    // printf("seqs %d %d only pre-computation time: %f seconds.\n", seq1_len, seq2_len, parse_elapsed_time);

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
                for (int nucj = seq_len-1; nucj>=0; --nucj) {
                    seq->next_pair[nuci][nucj] = next;
                    // cout << i << " " << nuci << " " << nucj << "  " << seq->next_pair[nuci][nucj] << endl;
                    if (_allowed_pairs[nuci][seq->nucs[nucj]]) next = nucj;
                }
            }
        }
    }
    cout << "test next pair" << endl;

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
    // bestC[0][0].set(0, 0, 0, 0, 0, MANNER_NONE, MANNER_NONE, 0.0, 0.0, MANNER_ALN, MANNER_ALN);
    // bestC[1][1].set(1, 0, 0, -v_score_external_unpaired(0, 0), 0, MANNER_NONE, MANNER_C_eq_C_plus_U_INS1, 1.0, weight, MANNER_ALN, MANNER_INS1);
    // bestC[1][0].set(0, 0, 0, 0, -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_INS2, 1.0, weight, MANNER_ALN, MANNER_INS2);
    // bestC[2][1].set(1, 0, 0, -v_score_external_unpaired(0, 0), -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_ALN, 0.0, 0.0, MANNER_ALN, MANNER_ALN);

#else
    bestC[0][0].alnobj.set(0, 0, 0, 0, 0, MANNER_NONE, MANNER_NONE, xlog(1.0), weight*xlog(1.0), MANNER_ALN, MANNER_ALN);

    if (0 >= hmmalign.low_bounds[1] && 0 <= hmmalign.up_bounds[1]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_INS1, 1, 0, true);
        bestC[1][0].ins1obj.set(1, 0, 0, -v_score_external_unpaired(0, 0), 0, MANNER_NONE, MANNER_C_eq_C_plus_U_INS1, alignscore, weight*alignscore, MANNER_ALN, MANNER_INS1);
    }
    if (1 >= hmmalign.low_bounds[0] && 1 <= hmmalign.up_bounds[0]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_INS2, 0, 1, true);
        bestC[0][1].ins2obj.set(0, 0, 0, 0, -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_INS2, alignscore, weight*alignscore, MANNER_ALN, MANNER_INS2);
    }
    if (1 >= hmmalign.low_bounds[1] && 1 <= hmmalign.up_bounds[1]) {
        float alignscore = hmmalign.get_trans_emit_prob0(MANNER_ALN, MANNER_ALN, 1, 1, true);
        bestC[1][1].alnobj.set(1, 0, 0, -v_score_external_unpaired(0, 0), -v_score_external_unpaired(0, 0), MANNER_NONE, MANNER_C_eq_C_plus_U_ALN, alignscore, weight*alignscore, MANNER_ALN, MANNER_ALN);
    }
#endif

    // from left to right
    mem = GetProcessMemory();
    for(int s = 1; s < seq1_len + seq2_len - 1; ++s) {
        // if (s > 2400) verbose = true;

        if (s%10 == 1)
            cout << "s: " << s << " VmPeak: " << mem.VmPeak  / 1024.0 / 1024.0 << endl;

        if (hmmalign.max_j1[s] == -1) continue;

        // boustronphen
        // int center_j1 = int(s * seq1_len / (float)sum_len);
        // int j1 = center_j1 + (1-2*(j%2)) * int((j+1)/2);

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

                // single seq folding subopt
                if (!seq1_H_pairs[j1next][j1].valid) continue;

                j2next = j2;
                while (j2next != -1) {
                    newscore2 = hairpinScore(j2, j2next, seq2);
                    j2next = newscore2.first;
                    if (j2next == -1) break;

                    // single seq folding subopt
                    if (!seq2_H_pairs[j2next][j2].valid) continue;

                    // hmm constraint, right bracket   
                    if (j2next > hmmalign.up_bounds[j1next] || j2next < hmmalign.low_bounds[j1next]) continue;

                    if (verbose) cout << "hairpin candidates: " << j1 << " " << j1next << " " << j2 << " " << j2next << endl;

                    float alignscore = get_hmm_score(j1, j1next, j2, j2next, 2);
                    if (alignscore > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1][j1next][j2][j2next].alnobj,
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_ALN, 
                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                    float ins1score = get_hmm_score(j1, j1next, j2+1, j2next-1, 0);
                    if (ins1score > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1][j1next][j2][j2next].ins1obj,
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_INS1, 
                                        ins1score, MANNER_INS1, MANNER_INS1, weight, verbose);

                    float ins2score = get_hmm_score(j1+1, j1next-1, j2, j2next, 1);
                    if (ins2score > LOG_OF_ZERO)
                        update_if_better(j1, j1next, j2, j2next,
                                        bestH[j1][j1next][j2][j2next].ins2obj,
                                        newscore1.second, newscore2.second, MANNER_NONE, MANNER_H_INS2, 
                                        ins2score, MANNER_INS2, MANNER_INS2, weight, verbose);

                    // 2. generate p(i, j)
                    // only check hmm constraints for right side when convert to P
                    // if (j2 >= hmmalign.low_bounds[j1] && j2 <= hmmalign.up_bounds[j1])
                    if (seq1_P_pairs[j1next][j1].valid && seq2_P_pairs[j2next][j2].valid) {
                        if (verbose) cout << "H to P: " << j1 << " " << j1next << " " << j2 << " " << j2next << endl; 

                        if (alignscore > LOG_OF_ZERO)
                            update_if_better(j1, j1next, j2, j2next,
                                            bestP[j1][j1next][j2][j2next].alnobj,
                                            newscore1.second, newscore2.second, MANNER_H_ALN, MANNER_HAIRPIN_ALN, 
                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                        if (ins1score > LOG_OF_ZERO)
                            update_if_better(j1, j1next, j2, j2next,
                                            bestP[j1][j1next][j2][j2next].ins1obj,
                                            newscore1.second, newscore2.second, MANNER_H_INS1, MANNER_HAIRPIN_INS1, 
                                            ins1score, MANNER_INS1, MANNER_INS1, weight, verbose);

                        if (ins2score > LOG_OF_ZERO)
                            update_if_better(j1, j1next, j2, j2next,
                                            bestP[j1][j1next][j2][j2next].ins2obj,
                                            newscore1.second, newscore2.second, MANNER_H_INS2, MANNER_HAIRPIN_INS2, 
                                            ins2score, MANNER_INS2, MANNER_INS2, weight, verbose);
                    } // H->P(i, j)
                }
            }
        }
        
        // beam of Multi
        // for every state in Multi[j]
        //   1. extend (i, j) to (i, jnext)
        //   2. generate P (i, j)
        if (verbose) cout << "beam of Multi" << endl;
        {
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                int nucj1 = seq1->nucs[j1];
                int i1 = seq1->next_pair[nucj1][0];
                while (i1 < j1 && i1 != -1) { // only consider positions where (i1, j1) can form pairs
                    // single seq folding subopt
                    if (!seq1_Multi_pairs[j1][i1].valid || !seq1_P_pairs[j1][i1].valid) {
                        i1 = seq1->next_pair[nucj1][i1];
                        continue;
                    }

                    int nucj2 = seq2->nucs[j2];
                    int i2 = seq2->next_pair[nucj2][max(0, hmmalign.low_bounds[i1]-1)];
                    int max_i2 = min(j2-1, hmmalign.up_bounds[i1]);
                    while (i2 <= max_i2 && i2 != -1) { // only consider positions where (i2, j2) can form pairs
                        // single seq folding subopt
                        if (!seq2_Multi_pairs[j2][i2].valid || !seq2_P_pairs[j2][i2].valid) {
                            i2 = seq2->next_pair[nucj2][i2];
                            continue;
                        }

                        State3& state3 = bestMulti[i1][j1][i2][j2];

                        // 2. generate P (i, j)
                        // single seq folding subopt
                        {
                            if (verbose) cout << "Multi to P: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;                            

                            int score1 = multiloop2Pscore(i1, j1, seq1);
                            int score2 = multiloop2Pscore(i2, j2, seq2);

                            if (state3.alnobj.manner != MANNER_NONE)
                                update_if_better(i1, j1, i2, j2, bestP[i1][j1][i2][j2].alnobj,
                                                state3.alnobj.seq1foldscore + score1, 
                                                state3.alnobj.seq2foldscore + score2,
                                                state3.alnobj.manner, MANNER_P_eq_MULTI_ALN,
                                                state3.alnobj.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                            if (state3.ins1obj.manner != MANNER_NONE)
                                update_if_better(i1, j1, i2, j2, bestP[i1][j1][i2][j2].ins1obj,
                                                state3.ins1obj.seq1foldscore + score1, 
                                                state3.ins1obj.seq2foldscore + score2,
                                                state3.ins1obj.manner, MANNER_P_eq_MULTI_INS1,
                                                state3.ins1obj.alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);
                            if (state3.ins2obj.manner != MANNER_NONE)
                                update_if_better(i1, j1, i2, j2, bestP[i1][j1][i2][j2].ins2obj,
                                                state3.ins2obj.seq1foldscore + score1, 
                                                state3.ins2obj.seq2foldscore + score2,
                                                state3.ins2obj.manner, MANNER_P_eq_MULTI_INS2,
                                                state3.ins2obj.alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                        } // 2. generate P (i, j)

                        /*
                        State state;
                        for (int m=0; m<3; m++){
                            switch (m)
                            {
                                case 0:
                                    state = state3.ins1obj;
                                    break;
                                case 1:
                                    state = state3.ins2obj;
                                    break;
                                case 2:
                                    state = state3.alnobj;
                                    break;
                                
                                default:
                                    break;
                            }
                            if (state.manner == MANNER_NONE) continue;
                        
                            // 1. extend (i, j) to (i, jnext)
                            {
                                tuple<int, int, char, int> result1 = multiloopUnpairedScore(i1, j1, seq1, &state.trace1);
                                tuple<int, int, char, int> result2 = multiloopUnpairedScore(i2, j2, seq2, &state.trace2);

                                int j1next = get<0>(result1);
                                int j2next = get<0>(result2);

                                int new_seq1_l2 = get<3>(result1);
                                int new_seq2_l2 = get<3>(result2);

                                if (verbose) cout << "beamMulti jnext: " << m << " "  << i1  <<  " " << j1 << " " << i2 << " " << j2 << " "  << j1next  <<  " " << j2next << endl; 
                                
                                float alignscore;
                                if (m == 2) {
                                    if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) { 
                                        // single seq folding subopt
                                        if (seq1_Multi_pairs[j1next][i1] && seq2_Multi_pairs[j2next][i2]) {
                                            
                                            alignscore = get_hmm_score_right(j1, j1next, j2, j2next, 2, 2);
                                            alignscore = xlog_mul(state.alignscore, alignscore);
                                            
                                            update_if_better(i1, j1next, i2, j2next, bestMulti[i1][j1next][i2][j2next].alnobj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                        } 
                                    }
                                    if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) {
                                        // single seq folding subopt
                                        if (seq1_Multi_pairs[j1next][i1]) {

                                            alignscore = get_hmm_score_right(j1, j1next, j2, j2, 2, 2);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1next, i2, j2,  bestMulti[i1][j1next][i2][j2].alnobj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                                        }
                                    }
                                    if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) {
                                        // single seq folding subopt
                                        if (seq2_Multi_pairs[j2next][i2]) {
                                            
                                            alignscore = get_hmm_score_right(j1, j1, j2, j2next, 2, 2);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1, i2, j2next, bestMulti[i1][j1][i2][j2next].alnobj,
                                                            state.seq1foldscore,
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_ALN,
                                                            state.trace1.paddings.l1, state.trace1.paddings.l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose); 
                                        }
                                    }
                                }
                                        
                                else if (m == 0){
                                    if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
                                        // single seq folding subopt
                                        if (seq1_Multi_pairs[j1next][i1] && seq2_Multi_pairs[j2next][i2]) {
                                            
                                            alignscore = get_hmm_score_right(j1, j1next, j2, j2next-1, 0, 0);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1next, i2, j2next, bestMulti[i1][j1next][i2][j2next].ins1obj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                                        }
                                    }
                                    if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) {
                                        // single seq folding subopt
                                        if (seq1_Multi_pairs[j1next][i1]) {

                                            alignscore = get_hmm_score_right(j1, j1next, j2, j2-1, 0, 0);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1next, i2, j2, bestMulti[i1][j1next][i2][j2].ins1obj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                            alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                                        }
                                    }
                                    if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) {
                                        // single seq folding subopt
                                        if (seq2_Multi_pairs[j2next][i2]) {

                                            alignscore = get_hmm_score_right(j1, j1, j2, j2next-1, 0, 0);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1, i2, j2next, bestMulti[i1][j1][i2][j2next].ins1obj,
                                                            state.seq1foldscore,
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS1,
                                                            state.trace1.paddings.l1, state.trace1.paddings.l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_INS1, MANNER_INS1, weight, verbose); 
                                        }
                                    }
                                }       
                                else if (m == 1) {
                                    if (j1next != -1 && j2next != -1 && j2next >= hmmalign.low_bounds[j1next] && j2next <= hmmalign.up_bounds[j1next]) {
                                        // single seq folding subopt
                                        if (seq1_Multi_pairs[j1next][i1] && seq2_Multi_pairs[j2next][i2]) {

                                            alignscore = get_hmm_score_right(j1, j1next-1, j2, j2next, 1, 1);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1next, i2, j2next, bestMulti[i1][j1next][i2][j2next].ins2obj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore + get<1>(result2),
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, new_seq2_l2,
                                                            alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                                        }
                                    }
                                    if (j1next != -1 && j2 >= hmmalign.low_bounds[j1next] && j2 <= hmmalign.up_bounds[j1next]) {
                                        // single seq folding subopt
                                        if (seq1_Multi_pairs[j1next][i1]) {

                                            alignscore = get_hmm_score_right(j1, j1next-1, j2, j2, 1, 1);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1next, i2, j2,  bestMulti[i1][j1next][i2][j2].ins2obj,
                                                            state.seq1foldscore + get<1>(result1),
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_MULTI_eq_MULTI_plus_U_INS2,
                                                            state.trace1.paddings.l1, new_seq1_l2,
                                                            state.trace2.paddings.l1, state.trace2.paddings.l2,
                                                            alignscore, MANNER_INS2, MANNER_INS2, weight, verbose); 
                                        }
                                    }
                                    if (j2next != -1 && j2next >= hmmalign.low_bounds[j1] && j2next <= hmmalign.up_bounds[j1]) {
                                        // single seq folding subopt
                                        if (seq2_Multi_pairs[j2next][i2]) {

                                            alignscore = get_hmm_score_right(j1, j1-1, j2, j2next, 1, 1);
                                            alignscore = xlog_mul(state.alignscore, alignscore);

                                            update_if_better(i1, j1, i2, j2next, bestMulti[i1][j1][i2][j2next].ins2obj,
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
                            // only check hmm constrain when convert to P? check in dynalign
                            // if (j2 >= hmmalign.low_bounds[j1] && j2 <= hmmalign.up_bounds[j1]) 
                            // single seq folding subopt
                            if (seq1_P_pairs[j1][i1] && seq2_P_pairs[j2][i2]){
                                if (verbose) cout << "Multi to P: " << m <<  " " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;                            

                                int score1 = multiloop2Pscore(i1, j1, seq1);
                                int score2 = multiloop2Pscore(i2, j2, seq2);

                                if (m == 2)
                                    update_if_better(i1, j1, i2, j2, bestP[i1][j1][i2][j2].alnobj,
                                                    state.seq1foldscore + score1, 
                                                    state.seq2foldscore + score2,
                                                    state.manner, MANNER_P_eq_MULTI_ALN,
                                                    state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                                else if (m == 0)
                                    update_if_better(i1, j1, i2, j2, bestP[i1][j1][i2][j2].ins1obj,
                                                    state.seq1foldscore + score1, 
                                                    state.seq2foldscore + score2,
                                                    state.manner, MANNER_P_eq_MULTI_INS1,
                                                    state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                                else if (m == 1)
                                    update_if_better(i1, j1, i2, j2, bestP[i1][j1][i2][j2].ins2obj,
                                                    state.seq1foldscore + score1, 
                                                    state.seq2foldscore + score2,
                                                    state.manner, MANNER_P_eq_MULTI_INS2,
                                                    state.alignscore, state.startHMMstate, state.endHMMstate, weight, verbose);
                            } // 2. generate P (i, j)
                        } // loop state ALN/INS1/INS2
                        */
                        i2 = seq2->next_pair[nucj2][i2];
                    } // while loop i2
                    i1 = seq1->next_pair[nucj1][i1];
                } // while loop i1
            } // for loop j1
        } // beam of Multi

        
        // beam of P
        // for every state in P[j]
        //   1. generate new helix/bulge
        //   2. M = P
        //   3. M2 = M + P
        //   4. C = C + P
        if (verbose) cout << "beam of P" << endl;
        {
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                int nucj1 = seq1->nucs[j1];
                int i1 = seq1->next_pair[nucj1][0];
                while (i1 != -1 && i1 < j1) { // only consider positions where (i1, j1) can form pairs
                    // single seq folding subopt
                    if (!seq1_P_pairs[j1][i1].valid) {
                        i1 = seq1->next_pair[nucj1][i1];
                        continue;
                    }
                
                    int nucj2 = seq2->nucs[j2];
                    int i2 = seq2->next_pair[nucj2][max(1, hmmalign.low_bounds[i1]) - 1];
                    int max_i2 = min(j2-1, hmmalign.up_bounds[i1]);
                    while (i2 != -1 && i2 <= max_i2) { // only consider positions where (i2, j2) can form pairs
                        // single seq folding subopt
                        if (!seq2_P_pairs[j2][i2].valid) {
                            i2 = seq2->next_pair[nucj2][i2];
                            continue;
                        }

                        State3& state3 = bestP[i1][j1][i2][j2];

                        // 1. generate new helix / single_branch
                        // new state is of shape p..i..j..q
                        // Note: p >= 1
                        {
                            int nuci1 = seq1->nucs[i1];
                            int nuci1_1 = seq1->nucs[i1 - 1];
                            int nucj1 = seq1->nucs[j1];
                            int nucj1p1 = seq1->nucs[j1 + 1];

                            int nuci2 = seq2->nucs[i2];
                            int nuci2_1 = seq2->nucs[i2 - 1];
                            int nucj2 = seq2->nucs[j2];
                            int nucj2p1 = seq2->nucs[j2 + 1];

                            // p1<i1, q1>j1; p2<i2, q2>j2
                            // ALN state 
                            State& state = state3.alnobj;
                            if (state.manner != MANNER_NONE) {
                                int min_p1 = max(i1 - MAX_LOOP_LEN - 1, 1); // SINGLE_MAX_LEN
                                int m = 2;
                                for (int p1 = i1-1; p1 >= min_p1; --p1) {
                                    int nucp1 = seq1->nucs[p1];
                                    int nucp1p1 = seq1->nucs[p1 + 1];
                                    int q1 = seq1->next_pair[nucp1][j1];

                                    // int i1_p1 = i1 - p1;
                                    // while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                                    while (q1 != -1 && (q1 - j1 - 1 <= MAX_LOOP_LEN)) {
                                        // single seq folding subopt
                                        if (!seq1_P_pairs[q1][p1].valid) {
                                            q1 = seq1->next_pair[nucp1][q1];
                                            continue;
                                        }

                                        int nucq1 = seq1->nucs[q1];
                                        int nucq1_1 = seq1->nucs[q1 - 1];
                                        // int p2p1 = P2PScore(p1,q1,i1,j1,nucp1,nucp1p1,nucq1_1,nucq1,nuci1_1,nuci1,nucj1,nucj1p1); 
                                        // assert (p2p1 == seq1_internal[i1][j1][p1][q1]);
                                        int p2p1 = seq1_internal[i1][j1][p1][q1];

                                        int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                                        int min_p2 = max(hmmalign.low_bounds[p1], max(1, i2 - MAX_LOOP_LEN - 1));
                                        for (int p2 = max_p2; p2 >= min_p2; --p2) { // for-loop seq2 p2
                                            int nucp2 = seq2->nucs[p2];
                                            int nucp2p1 = seq2->nucs[p2 + 1];
                                            int q2 = seq2->next_pair[nucp2][j2];

                                            // speed up
                                            if (q2 == -1 || q2 > hmmalign.up_bounds[q1] || (q2 - j2 - 1 > MAX_LOOP_LEN)) continue;
                                            if (q2 <= hmmalign.low_bounds[q1]-1)
                                                q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1]-1]; 

                                            int i2_p2 = i2 - p2;
                                            // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                            while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= MAX_LOOP_LEN)) {
                                                // single seq folding subopt
                                                if (!seq2_P_pairs[q2][p2].valid) {
                                                    q2 = seq2->next_pair[nucp2][q2];
                                                    continue;
                                                }

                                                int nucq2 = seq2->nucs[q2];
                                                int nucq2_1 = seq2->nucs[q2 - 1];
                                                int q2_q1 = q2-hmmalign.low_bounds[q1];

                                                // int p2p2 = P2PScore(p2,q2,i2,j2,nucp2,nucp2p1,nucq2_1,nucq2,nuci2_1,nuci2,nucj2,nucj2p1);
                                                // assert (p2p2 == seq2_internal[i2][j2][p2][q2]);
                                                int p2p2 = seq2_internal[i2][j2][p2][q2];

                                                if (verbose) cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;

                                                pair<float, HMMManner> pre_alignscore, post_alignscore;
                                                float pre_align_trans, post_align_trans, alignscore;
                                                // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj), si == sj == ALN/INS1/INS2
                                                {
                                                    // si == sj == ALN
                                                    // align cost = (p, i) + (i, j) + (j, q)
                                                    pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                                                    post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                                                    
                                                    alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                    alignscore = xlog_mul(alignscore, post_align_trans);
                                    
                                                    if (alignscore > LOG_OF_ZERO)
                                                        update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].alnobj, // bestPALN[q1+q2][get_key1(q1, p1, p2)],  // helix, one branch using one state MANNER_SINGLE
                                                                        state.seq1foldscore + p2p1,
                                                                        state.seq2foldscore + p2p2,
                                                                        state.manner, MANNER_SINGLE_ALN, 
                                                                        static_cast<char>(i1 - p1), q1 - j1,
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                                                    // si == sj == INS1
                                                    // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                                    pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                                                    post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0);
                                                    
                                                    alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                    alignscore = xlog_mul(alignscore, post_align_trans); 

                                                    if (alignscore > LOG_OF_ZERO)
                                                        update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins1obj, // bestPINS1[q1+q2][get_key1(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                        state.seq1foldscore + p2p1,
                                                                        state.seq2foldscore + p2p2,
                                                                        state.manner, MANNER_SINGLE_INS1, 
                                                                        static_cast<char>(i1 - p1), q1 - j1,
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);

                                                    // si == sj == INS2
                                                    // align cost = (p, i) + (i, j) + (j, q);
                                                    pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, m);
                                                    post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, m, 1);
                                                    
                                                    alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                    alignscore = xlog_mul(alignscore, post_align_trans);

                                                    if (alignscore > LOG_OF_ZERO)
                                                        update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins2obj, // bestPINS2[q1+q2][get_key1(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                        state.seq1foldscore + p2p1,
                                                                        state.seq2foldscore + p2p2,
                                                                        state.manner, MANNER_SINGLE_INS2, 
                                                                        static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                                        static_cast<char>(i2 - p2), q2 - j2,
                                                                        alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                                } // ALN -> ALN/INS1/INS2

                                                q2 = seq2->next_pair[nucp2][q2];
                                            } // while loop q2
                                        } // for loop p2
                                        q1 = seq1->next_pair[nucp1][q1];
                                    } // while loop q1   
                                } // for loop p1

                                // 2. M = P
                                // accessible pairs in multiloop must be aligned
                                // single seq folding subopt
                                if (seq1_M_pairs[j1][i1].valid && seq2_M_pairs[j2][i2].valid) {
                                    if (verbose) cout << "P to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                                    int newscore1 = branch_score(i1, j1, seq1);
                                    int newscore2 = branch_score(i2, j2, seq2);
                                    update_if_better(i1, j1, i2, j2, bestM[i1][j1][i2][j2].alnobj, // bestM[key1][newi2][newj2], // beamM[j2][make_pair(i1, i2)], 
                                                    state.seq1foldscore + newscore1,
                                                    state.seq2foldscore + newscore2, 
                                                    state.manner, MANNER_M_eq_P,
                                                    state.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                } // 2. M = P

                                int k1 = i1 - 1;
                                int k2 = i2 - 1;
                                if (k2 >= hmmalign.low_bounds[k1] && k2 <= hmmalign.up_bounds[k1]) {
                                    // 3. M2 = M + P 
                                    // accessible pairs in multiloop must be aligned
                                    if (k1 > 1 && k2 > 1) {
                                        int newscore1 = branch_score(i1, j1, seq1) + state.seq1foldscore;
                                        int newscore2 = branch_score(i2, j2, seq2) + state.seq2foldscore;

                                        for (int newi1 = 1; newi1 < k1; newi1++) { // TODO: k1 - 3, no sharp turn
                                            // single seq folding subopt
                                            if (!seq1_M_pairs[k1][newi1].valid || !seq1_M2_pairs[j1][newi1].valid) continue;

                                            for (int newi2 = hmmalign.low_bounds[newi1]; newi2 <= min(k2-1, hmmalign.up_bounds[newi1]); newi2++) {
                                                
                                                // single seq folding subopt
                                                if (!seq2_M_pairs[k2][newi2].valid || !seq2_M2_pairs[j2][newi2].valid) continue;

                                                if (verbose) cout << "M2=M+P: " << i1 << " " << j1 << " "  << i2 << " " << j2 <<  " " << newi1 << " " << j1 << " "  << newi2 << " " << j2  << endl;

                                                State3& mstate3 = bestM[newi1][k1][newi2][k2];
                                                {
                                                    State& mstate = mstate3.alnobj;
                                                    if (mstate.endHMMstate != HMMMANNER_NONE) {
                                                        float alignscore = xlog_mul(mstate.alignscore, hmmalign.trans_probs[2][2]); 
                                                        alignscore = xlog_mul(alignscore, state.alignscore);
                                                        
                                                        if (alignscore > LOG_OF_ZERO)
                                                            update_if_better(newi1, j1, newi2, j2, bestM2[newi1][j1][newi2][j2],
                                                                            mstate.seq1foldscore + newscore1,
                                                                            mstate.seq2foldscore + newscore2,
                                                                            mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                                            alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                                                    }
                                                }
                                                {
                                                    State& mstate = mstate3.ins1obj;
                                                    if (mstate.endHMMstate != HMMMANNER_NONE) {
                                                        float alignscore = xlog_mul(mstate.alignscore, hmmalign.trans_probs[0][2]); 
                                                        alignscore = xlog_mul(alignscore, state.alignscore);
                                                        
                                                        if (alignscore > LOG_OF_ZERO)
                                                            update_if_better(newi1, j1, newi2, j2, bestM2[newi1][j1][newi2][j2],
                                                                            mstate.seq1foldscore + newscore1,
                                                                            mstate.seq2foldscore + newscore2,
                                                                            mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                                            alignscore, MANNER_INS1, MANNER_ALN, weight, verbose);
                                                    }
                                                }
                                                {
                                                    State& mstate = mstate3.ins2obj;
                                                    if (mstate.endHMMstate != HMMMANNER_NONE) {
                                                        float alignscore = xlog_mul(mstate.alignscore, hmmalign.trans_probs[1][2]); 
                                                        alignscore = xlog_mul(alignscore, state.alignscore);
                                                        
                                                        if (alignscore > LOG_OF_ZERO)
                                                            update_if_better(newi1, j1, newi2, j2, bestM2[newi1][j1][newi2][j2],
                                                                            mstate.seq1foldscore + newscore1,
                                                                            mstate.seq2foldscore + newscore2,
                                                                            mstate.manner, MANNER_M2_eq_M_plus_P, k1, k2,
                                                                            alignscore, MANNER_INS2, MANNER_ALN, weight, verbose);
                                                    }
                                                }
                                            }
                                        }
                                    } // 3. M2 = M + P 
                                    
                                    // 4. C = C + P
                                    // external pairs must be aligned
                                    if (k1 >= 0 && k2 >= 0) {
                                        int newscore1 = external_paired_score(k1, j1, seq1);
                                        int newscore2 = external_paired_score(k2, j2, seq2);
                                        {
                                            State& prefix_C = bestC[k1][k2].alnobj;
                                            if (prefix_C.endHMMstate != HMMMANNER_NONE) {
                                                
                                                if (verbose) cout << "C+P: "<< i1 << " " << j1 << " "  << i2 << " " << j2  <<  " " << state.manner << " " << state.premanner << " " << newscore1 << " " << newscore2 << endl;

                                                float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[2][2]); 
                                                alignscore = xlog_mul(alignscore, state.alignscore); 

                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(0, j1, 0, j2, bestC[j1][j2].alnobj,
                                                                    prefix_C.seq1foldscore + state.seq1foldscore + newscore1,
                                                                    prefix_C.seq2foldscore + state.seq2foldscore + newscore2, 
                                                                    prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                                    k1, k2,
                                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);  
                                            }
                                        }
                                        {
                                            State& prefix_C = bestC[k1][k2].ins1obj;
                                            if (prefix_C.endHMMstate != HMMMANNER_NONE) {
                                                if (verbose) cout << "C+P: "<< i1 << " " << j1 << " "  << i2 << " " << j2  <<  " " << state.manner << " " << state.premanner << " " << newscore1 << " " << newscore2 << endl;

                                                float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[0][2]); 
                                                alignscore = xlog_mul(alignscore, state.alignscore); 

                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(0, j1, 0, j2, bestC[j1][j2].alnobj,
                                                                    prefix_C.seq1foldscore + state.seq1foldscore + newscore1,
                                                                    prefix_C.seq2foldscore + state.seq2foldscore + newscore2, 
                                                                    prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                                    k1, k2,
                                                                    alignscore, MANNER_INS1, MANNER_ALN, weight, verbose);  
                                            }
                                        }
                                        {
                                            State& prefix_C = bestC[k1][k2].ins2obj;
                                            if (prefix_C.endHMMstate != HMMMANNER_NONE) {
                                                if (verbose) cout << "C+P: "<< i1 << " " << j1 << " "  << i2 << " " << j2  <<  " " << state.manner << " " << state.premanner << " " << newscore1 << " " << newscore2 << endl;

                                                float alignscore = xlog_mul(prefix_C.alignscore, hmmalign.trans_probs[1][2]); 
                                                alignscore = xlog_mul(alignscore, state.alignscore); 

                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(0, j1, 0, j2, bestC[j1][j2].alnobj,
                                                                    prefix_C.seq1foldscore + state.seq1foldscore + newscore1,
                                                                    prefix_C.seq2foldscore + state.seq2foldscore + newscore2, 
                                                                    prefix_C.manner, MANNER_C_eq_C_plus_P,
                                                                    k1, k2,
                                                                    alignscore, MANNER_INS2, MANNER_ALN, weight, verbose);  
                                            }
                                        }
                                                 
                                    } // 4. C = C + P
                                } // 3. M2 = M + P and 4. C = C + P
                            } // ALN state 

                            State& in2state = state3.ins2obj;
                            if (in2state.manner != MANNER_NONE) {
                                int p1 = i1;
                                int q1 = j1;

                                int nucp1 = seq1->nucs[p1];
                                int nucp1p1 = seq1->nucs[p1 + 1];
                                int nucq1 = seq1->nucs[q1];
                                int nucq1_1 = seq1->nucs[q1 - 1];
                                int p2p1 = 0; 

                                int m = 1;
                                int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                                int min_p2 = max(hmmalign.low_bounds[p1], max(1, i2 - MAX_LOOP_LEN - 1));
                                for (int p2 = max_p2; p2 >= min_p2; --p2) { // for-loop seq2 p2 
                                    int nucp2 = seq2->nucs[p2];
                                    int nucp2p1 = seq2->nucs[p2 + 1];
                                    int q2 = seq2->next_pair[nucp2][j2];

                                    // speed up
                                    // if (q2 <= max(j2, hmmalign.low_bounds[q1]-1))
                                    //     q2 = seq2->next_pair[nucp2][max(j2, hmmalign.low_bounds[q1]-1)]; 
                                    if (q2 == -1 || q2 > hmmalign.up_bounds[q1] || (q2 - j2 - 1 > MAX_LOOP_LEN)) continue;
                                    if (q2 <= hmmalign.low_bounds[q1]-1)
                                        q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1]-1]; 

                                    // int i2_p2 = i2 - p2;
                                    // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                    while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= MAX_LOOP_LEN)) {
                                        // single seq folding subopt
                                        if (!seq2_P_pairs[q2][p2].valid){
                                            q2 = seq2->next_pair[nucp2][q2];
                                            continue;
                                        }

                                        int nucq2 = seq2->nucs[q2];
                                        int nucq2_1 = seq2->nucs[q2 - 1];
                                        // int p2p2 = P2PScore(p2,q2,i2,j2,nucp2,nucp2p1,nucq2_1,nucq2,nuci2_1,nuci2,nucj2,nucj2p1);
                                        // assert (p2p2 == seq2_internal[i2][j2][p2][q2]);
                                        int p2p2 = seq2_internal[i2][j2][p2][q2];

                                        if (verbose) cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;

                                        pair<float, HMMManner> pre_alignscore, post_alignscore;
                                        float pre_align_trans, post_align_trans, alignscore;
                                        // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj), si == sj == ALN/INS1/INS2
                                        {
                                            // si == sj == ALN
                                            // align cost = (p, i) + (i, j) + (j, q)
                                            pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                                            
                                            alignscore = xlog_mul(pre_align_trans, in2state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans);
                            
                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].alnobj, // bestPALN[q1+q2][get_key1(q1, p1, p2)],  // helix, one branch using one state MANNER_SINGLE
                                                                in2state.seq1foldscore + p2p1,
                                                                in2state.seq2foldscore + p2p2,
                                                                in2state.manner, MANNER_SINGLE_ALN, 
                                                                static_cast<char>(i1 - p1), q1 - j1,
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                                            // si == sj == INS1
                                            // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                            pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0);
                                            
                                            alignscore = xlog_mul(pre_align_trans, in2state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans); 

                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins1obj, // bestPINS1[q1+q2][get_key1(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                in2state.seq1foldscore + p2p1,
                                                                in2state.seq2foldscore + p2p2,
                                                                in2state.manner, MANNER_SINGLE_INS1, 
                                                                static_cast<char>(i1 - p1), q1 - j1,
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);

                                            // si == sj == INS2
                                            // align cost = (p, i) + (i, j) + (j, q);
                                            pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, m);
                                            post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, m, 1);
                                            
                                            alignscore = xlog_mul(pre_align_trans, in2state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans);

                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins2obj, // bestPINS2[q1+q2][get_key1(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                in2state.seq1foldscore + p2p1,
                                                                in2state.seq2foldscore + p2p2,
                                                                in2state.manner, MANNER_SINGLE_INS2, 
                                                                static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                        } 
                                        q2 = seq2->next_pair[nucp2][q2];
                                    } // while loop q2
                                } // for loop p2
                            } // INS2 state 

                            // INS1 state 
                            State& ins1state = state3.ins1obj;
                            if (ins1state.manner != MANNER_NONE) {
                                // fixed p2 q2
                                int p2 = i2;
                                int q2 = j2;

                                int nucp2 = seq2->nucs[p2];
                                int nucp2p1 = seq2->nucs[p2 + 1];
                                int nucq2 = seq2->nucs[q2];
                                int nucq2_1 = seq2->nucs[q2 - 1];
                                int p2p2 = 0;

                                int min_p1 = max(i1 - MAX_LOOP_LEN - 1, 1); // SINGLE_MAX_LEN
                                int m = 0;
                                for (int p1 = i1-1; p1 >= min_p1; --p1) {
                                    if (p2 < hmmalign.low_bounds[p1] || p2 > hmmalign.up_bounds[p1]) continue; // TODO

                                    int nucp1 = seq1->nucs[p1];
                                    int nucp1p1 = seq1->nucs[p1 + 1];
                                    int q1 = seq1->next_pair[nucp1][j1];

                                    int i1_p1 = i1 - p1;
                                    // while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                                    while (q1 != -1 && (q1 - j1 - 1 <= MAX_LOOP_LEN)) {
                                        if (q2 < hmmalign.low_bounds[q1] || q2 > hmmalign.up_bounds[q1]) { // TODO
                                            q1 = seq1->next_pair[nucp1][q1];
                                            continue;
                                        }

                                        // single seq folding subopt
                                        if (!seq1_P_pairs[q1][p1].valid) {
                                            q1 = seq1->next_pair[nucp1][q1];
                                            continue;
                                        }

                                        int nucq1 = seq1->nucs[q1];
                                        int nucq1_1 = seq1->nucs[q1 - 1];
                                        // int p2p1 = P2PScore(p1,q1,i1,j1,nucp1,nucp1p1,nucq1_1,nucq1,nuci1_1,nuci1,nucj1,nucj1p1); 
                                        // assert (p2p1 == seq1_internal[i1][j1][p1][q1]);
                                        int p2p1 = seq1_internal[i1][j1][p1][q1];

                                        if (verbose) cout << "P2P: " << m << " " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << p1 << " " << q1 << " " << p2 << " " << q2 << " " <<  p2p1 << " " << p2p2 << endl;

                                        pair<float, HMMManner> pre_alignscore, post_alignscore;
                                        float pre_align_trans, post_align_trans, alignscore;
                                        // (i1,j1;i2,j2;si,sj)->(p1,q1;p2,q2;si,sj), si == sj == ALN/INS1/INS2
                                        {
                                            // si == sj == ALN
                                            // align cost = (p, i) + (i, j) + (j, q)
                                            pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, m);
                                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2, m, 2);
                                            alignscore = xlog_mul(pre_align_trans, ins1state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans);
                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].alnobj, // bestPALN[q1+q2][get_key1(q1, p1, p2)],  // helix, one branch using one state MANNER_SINGLE
                                                                ins1state.seq1foldscore + p2p1,
                                                                ins1state.seq2foldscore + p2p2,
                                                                ins1state.manner, MANNER_SINGLE_ALN, 
                                                                static_cast<char>(i1 - p1), q1 - j1,
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                                            // si == sj == INS1
                                            // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                            pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, m);
                                            post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, m, 0); 
                                            alignscore = xlog_mul(pre_align_trans, ins1state.alignscore);
                                            alignscore = xlog_mul(alignscore, post_align_trans); 
                                            if (alignscore > LOG_OF_ZERO)
                                                update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins1obj, // bestPINS1[q1+q2][get_key1(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
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
                                                update_if_better(p1, q1, p2, q2, bestP[p1][q1][p2][q2].ins2obj, // bestPINS2[q1+q2][get_key1(q1, p1, p2)], // helix, one branch using one state MANNER_SINGLE
                                                                ins1state.seq1foldscore + p2p1,
                                                                ins1state.seq2foldscore + p2p2,
                                                                ins1state.manner, MANNER_SINGLE_INS2, 
                                                                static_cast<char>(i1 - p1), q1 - j1, // char?? int?? 
                                                                static_cast<char>(i2 - p2), q2 - j2,
                                                                alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                        } // ALN -> ALN/INS1/INS2
    
                                        q1 = seq1->next_pair[nucp1][q1];
                                    } // while loop q1   
                                } // for loop p1
                            } // INS2 state 
                        } // 1. generate new helix / single_branch

                        i2 = seq2->next_pair[nucj2][i2];
                    } // while loop i2 
                    
                    i1 = seq1->next_pair[nucj1][i1];
                } // while loop i1 
            } // for loop j1
        } 
        
        // beam of M2
        // for every state in M2[j]
        //   1. multi-loop  (by extending M2 on the left)
        //   2. M = M2
        if (verbose) cout << "beam of M2" << endl;
        {
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (int i1 = 1; i1 < j1; i1++) {
                    // single seq folding subopt
                    if (!seq1_M2_pairs[j1][i1].valid) continue;

                    for (int i2 = hmmalign.low_bounds[i1]; i2 <= min(j2-1, hmmalign.up_bounds[i1]); i2++) {
                        // single seq folding subopt
                        if (!seq2_M2_pairs[j2][i2].valid) continue;

                        if (verbose) cout << i1 << " " << j1 << " " << i2 << " " << j2 << " " << hmmalign.low_bounds[i1] << " " << hmmalign.up_bounds[i1] << " " << hmmalign.low_bounds[j1] << " " << hmmalign.up_bounds[j1] << endl ;

                        State& state = bestM2[i1][j1][i2][j2];
                        if (state.manner == MANNER_NONE) continue;

                        // 1. multi-loop
                        {
                            for (int p1 = i1-1; p1 >= max(i1 - SINGLE_MAX_LEN - 1, 1); --p1) {
                                int nucp1 = seq1->nucs[p1];
                                int q1 = seq1->next_pair[nucp1][j1];
                                
                                int i1_p1 = i1 - p1;
                                while (q1 != -1 && (i1_p1 + (q1 - j1) - 2 <= SINGLE_MAX_LEN)) {
                                // while (q1 != -1 && (q1 - j1 - 1 <= SINGLE_MAX_LEN)) {
                                    // single seq folding subopt
                                    if (!seq1_Multi_pairs[q1][p1].valid) {
                                        q1 = seq1->next_pair[nucp1][q1];
                                        continue;
                                    }

                                    // the current shape is p..i M2 j ..q
                                    int max_p2 = min(i2-1, hmmalign.up_bounds[p1]);
                                    int min_p2 = max(max(1, i2 - SINGLE_MAX_LEN - 1), hmmalign.low_bounds[p1]);
                                    for (int p2 = max_p2; p2 >= min_p2; --p2) {
                                        int nucp2 = seq2->nucs[p2];
                                        int q2 = seq2->next_pair[nucp2][j2];

                                        // speed up
                                        int i2_p2 = i2 - p2;
                                        if (q2 == -1 || q2 > hmmalign.up_bounds[q1] || (i2_p2 + (q2 - j2) - 2 > SINGLE_MAX_LEN)) continue;
                                        if (q2 < hmmalign.low_bounds[q1])
                                            q2 = seq2->next_pair[nucp2][hmmalign.low_bounds[q1]-1]; 
                                        
                                        while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (i2_p2 + (q2 - j2) - 2 <= SINGLE_MAX_LEN)) {
                                        // while (q2 <= hmmalign.up_bounds[q1] && q2 != -1 && (q2 - j2 - 1 <= SINGLE_MAX_LEN)) {
                                            // single seq folding subopt
                                            if (seq2_Multi_pairs[q2][p2].valid) {
                                                int newscore1 = multi_unpaired_score2(i1, j1, p1, q1, seq1);
                                                int newscore2 = multi_unpaired_score2(i2, j2, p2, q2, seq2);

                                                float pre_align_trans, post_align_trans, alignscore;

                                                if (verbose) cout << "M2 to Multi: " << i1 << " " << j1 << " " << i2 << " " << j2  << " " << p1 << " " << q1 << " " << p2 << " " << q2 << endl;

                                                // update bestMultiALN
                                                pre_align_trans = get_hmm_score_left(p1, i1, p2, i2, 2, 2);
                                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2, 2, 2);

                                                alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);

                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestMulti[p1][q1][p2][q2].alnobj,
                                                                    state.seq1foldscore + newscore1,
                                                                    state.seq2foldscore + newscore2, 
                                                                    state.manner, MANNER_MULTI_ALN,
                                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);

                                                // update bestMultiINS1
                                                // align cost = (p, i) + (i, j) + (j, q); p2+1, q2-1
                                                pre_align_trans = get_hmm_score_left(p1, i1, p2+1, i2, 0, 2);
                                                post_align_trans = get_hmm_score_right(j1, q1, j2, q2-1, 2, 0); 
                                                
                                                alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);

                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestMulti[p1][q1][p2][q2].ins1obj,
                                                                    state.seq1foldscore + newscore1,
                                                                    state.seq2foldscore + newscore2, 
                                                                    state.manner, MANNER_MULTI_INS1,
                                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_INS1, MANNER_INS1, weight, verbose);

                                                // update bestMultiINS2
                                                // align cost = (p, i) + (i, j) + (j, q); p1+1, q1-1
                                                pre_align_trans = get_hmm_score_left(p1+1, i1, p2, i2, 1, 2);
                                                post_align_trans = get_hmm_score_right(j1, q1-1, j2, q2, 2, 1); 
                                                
                                                alignscore = xlog_mul(pre_align_trans, state.alignscore);
                                                alignscore = xlog_mul(alignscore, post_align_trans);

                                                if (alignscore > LOG_OF_ZERO)
                                                    update_if_better(p1, q1, p2, q2, bestMulti[p1][q1][p2][q2].ins2obj,
                                                                    state.seq1foldscore + newscore1,
                                                                    state.seq2foldscore + newscore2, 
                                                                    state.manner, MANNER_MULTI_INS2,
                                                                    static_cast<char>(i1 - p1), q1 - j1,
                                                                    static_cast<char>(i2 - p2), q2 - j2,
                                                                    alignscore, MANNER_INS2, MANNER_INS2, weight, verbose);
                                            }
                                            q2 = seq2->next_pair[nucp2][q2];
                                        } // while loop enumerate q2
                                    } // for loop enumerate p2
                                    q1 = seq1->next_pair[nucp1][q1];
                                } // q1
                            } // p1
                        }// 1. multi-loop

                        // 2. M = M2
                        // single seq folding subopt
                        if (seq1_M_pairs[j1][i1].valid && seq2_M_pairs[j2][i2].valid){
                            if (verbose) cout << "M2 to M: " << i1 << " " << j1 << " " << i2 << " " << j2 << endl;
                            update_if_better(i1, j1, i2, j2, bestM[i1][j1][i2][j2].alnobj, // bestM[key1][newi2][newj2], //beamM[j2][make_pair(i1, i2)], 
                                            state.seq1foldscore,
                                            state.seq2foldscore,
                                            state.manner, MANNER_M_eq_M2,
                                            state.alignscore, MANNER_ALN, MANNER_ALN, weight, verbose);
                        } // 2. M = M2

                    } // for loop i2
                } // for loop i1
            } // for loop j1
        }

        // beam of M
        // for every state in M[j]
        //   1. M = M + unpaired 
        if (verbose) cout << "beam of M" << endl;
        {
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                for (int i1 = 1; i1 < j1; i1++) {
                    // single seq folding subopt
                    if (!seq1_M_pairs[j1][i1].valid) continue;

                    for (int i2 = hmmalign.low_bounds[i1]; i2 <= min(j2-1, hmmalign.up_bounds[i1]); i2++) {
                        // single seq folding subopt
                        if (!seq2_M_pairs[j2][i2].valid) continue;

                        State3& state3 = bestM[i1][j1][i2][j2];
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

                            if (verbose) cout << "M+U: " << i1 << " " << j1 << " " << i2 << " " << j2 << " " << state.i1 << " " << state.i2 << " " << state.j1 << endl;
                            // M+U: 14 35 18 42

                            float trans_emit_prob;
                            if (j1 < seq1->seq_len - 1 && j2 < seq2->seq_len - 1){
                                // hmm constraint
                                if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]){
                                    // single seq folding subopt
                                    if (seq1_M_pairs[j1+1][i1].valid && seq2_M_pairs[j2+1][i2].valid){
                                        trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_ALN, j1+1, j2+1, true);
                                        float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1+1, i2, j2+1, bestM[i1][j1+1][i2][j2+1].alnobj, // bestM[key1p1][newi2][newj2p1], // [s+2][j2+1][make_pair(i1, i2)],
                                                            state.seq1foldscore,
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_M_eq_M_plus_U_ALN,
                                                            alignscore, pre_manner, MANNER_ALN, weight, verbose);
                                    }
                                }
                            }
                            if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) {
                                // hmm constraint
                                if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) {
                                    // single seq folding subopt
                                    if (seq1_M_pairs[j1+1][i1].valid){
                                        trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_INS1, j1+1, j2, true);
                                        float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1+1, i2, j2, bestM[i1][j1+1][i2][j2].ins1obj, // bestM[key1p1][newi2][newj2p1], // [s+1][j2][make_pair(i1, i2)],
                                                            state.seq1foldscore,
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_M_eq_M_plus_U_INS1,
                                                            alignscore, pre_manner, MANNER_INS1, weight, verbose);
                                    }
                                }
                            }
                            if (j1 <= seq1->seq_len - 1 && j2 < seq2->seq_len - 1) {
                                // hmm constraint
                                if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) {
                                    // single seq folding subopt
                                    if (seq2_M_pairs[j2+1][i2].valid){
                                        trans_emit_prob = hmmalign.get_trans_emit_prob0(pre_manner, MANNER_INS2, j1, j2+1, true);
                                        float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                                        if (alignscore > LOG_OF_ZERO)
                                            update_if_better(i1, j1, i2, j2+1, bestM[i1][j1][i2][j2+1].ins2obj, // bestM[key1][newi2][newj2p1], //s+1][j2+1][make_pair(i1, i2)],
                                                            state.seq1foldscore,
                                                            state.seq2foldscore,
                                                            state.manner, MANNER_M_eq_M_plus_U_INS2,
                                                            alignscore, pre_manner, MANNER_INS2, weight, verbose);
                                    }
                                }
                            }
                        }
                    } // for loop i2
                } // for loop i1
            } // for loop j1
        }

        // beam of C
        // C = C + U
        if (verbose) cout << "beam of C" << endl;
        {
            for (int j1=hmmalign.min_j1[s]; j1<=hmmalign.max_j1[s]; j1++) {
                int j2 = s - j1;
                if (j2 < hmmalign.low_bounds[j1] || j2 > hmmalign.up_bounds[j1]) continue;

                int newscore1 = external_unpaired_score(j1 + 1, seq1);
                int newscore2 = external_unpaired_score(j2 + 1, seq2);

                State3& state3 = bestC[j1][j2];
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

                    if (verbose) cout << "C+U: " << m << " " << j1 << " "  << j2  << " " <<  j1 << " "  << j2  << " " << state.endHMMstate << endl;

                    float trans_emit_prob; //  alignscore;
                    if (j1 == seq1->seq_len - 1 && j2 == seq2->seq_len - 1) {
                        // cout << m << " " << j1 << " "  << j2  << " " << state.score << endl;

                        trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_ALN, j1+1, j2+1, true);
                        float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                        if (alignscore > LOG_OF_ZERO)
                            update_if_better(0, j1+1, 0, j2+1, bestC[j1+1][j2+1].alnobj, 
                                            state.seq1foldscore,
                                            state.seq2foldscore,
                                            state.manner, MANNER_C_eq_C_plus_U_ALN,
                                            alignscore, state.endHMMstate, MANNER_ALN, weight, verbose);
                        continue;
                    }

                    if ((j1 < seq1->seq_len - 1 && j2 < seq2->seq_len - 1)){ // ALN
                        // hmm constraint
                        if (j2+1 >= hmmalign.low_bounds[j1+1] && j2+1 <= hmmalign.up_bounds[j1+1]) {

                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_ALN, j1+1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1+1, 0, j2+1, bestC[j1+1][j2+1].alnobj, 
                                                state.seq1foldscore + newscore1,
                                                state.seq2foldscore + newscore2,
                                                state.manner, MANNER_C_eq_C_plus_U_ALN,
                                                alignscore, state.endHMMstate, MANNER_ALN, weight, verbose);
                        }
                    }
                    if (j1 < seq1->seq_len - 1 && j2 <= seq2->seq_len - 1) { // INS1
                        // hmm constraint
                        if (j2 >= hmmalign.low_bounds[j1+1] && j2 <= hmmalign.up_bounds[j1+1]) {

                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_INS1, j1+1, j2, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1+1, 0, j2, bestC[j1+1][j2].ins1obj,
                                                state.seq1foldscore + newscore1,
                                                state.seq2foldscore,
                                                state.manner, MANNER_C_eq_C_plus_U_INS1,
                                                alignscore, state.endHMMstate, MANNER_INS1, weight, verbose);
                        }
                    }
                    if (j2 < seq2->seq_len - 1 && j1 <= seq1->seq_len - 1) { // INS2
                        // hmm constraint
                        if (j2+1 >= hmmalign.low_bounds[j1] && j2+1 <= hmmalign.up_bounds[j1]) {

                            trans_emit_prob = hmmalign.get_trans_emit_prob0(state.endHMMstate, MANNER_INS2, j1, j2+1, true);
                            float alignscore = xlog_mul(state.alignscore, trans_emit_prob);
                            if (alignscore > LOG_OF_ZERO)
                                update_if_better(0, j1, 0, j2+1, bestC[j1][j2+1].ins2obj,
                                                state.seq1foldscore,
                                                state.seq2foldscore + newscore2,
                                                state.manner, MANNER_C_eq_C_plus_U_INS2,
                                                alignscore, state.endHMMstate, MANNER_INS2, weight, verbose);
                        }
                    }

                }           
            }     
        }

    } // for loop s
    cout << "end of loop s" << endl;

    gettimeofday(&parse_endtime, NULL);
    parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    printf("seqs %d %d only parse time: %f seconds.\n", seq1_len, seq2_len, parse_elapsed_time);  

    // // processMem_t mem;
    // // mem = GetProcessMemory();
    // // cout << "VmPeak: " << mem.VmPeak << endl;

    State &state = bestC[seq1_len][seq2_len].alnobj;
    float bestscore = state.score;
    cout << "inside: " << state.score << " " << state.seq1foldscore << " " << state.seq2foldscore << " " << state.alignscore << endl;
    // // cout << "inside: " << state.manner << " " << state.startHMMstate << " " << state.endHMMstate << endl;

    // backtrace
    tuple<string, string, string, string> ret = get_parentheses(*seq1, *seq2, hmmalign);
    cout << get<0>(ret) <<endl;
    cout << get<1>(ret) <<endl;
    cout << get<2>(ret) <<endl;
    cout << get<3>(ret) <<endl;

/*
#ifdef multilign
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
*/
    /*************************/
    // clear allocated space
    for (int i1=1; i1<seq1_len; i1++) {
        int nuci1 = sequences[0].nucs[i1];
        for (int j1 = i1+1; j1<seq1_len; j1++) {
            short i1_lowbound = hmmalign.low_bounds[i1];
            short i1_upbound = hmmalign.up_bounds[i1];
            short i1_range = i1_upbound - i1_lowbound + 1;
            
            int nucj1 = sequences[0].nucs[j1];

            for (int i2 = i1_lowbound; i2 <= i1_upbound; i2++) {
                short j1_lowbound = hmmalign.low_bounds[j1];
                
                if (_allowed_pairs[nuci1][nucj1]) {
                    if (seq1_H_pairs[j1][i1].valid) {
                        bestH[i1][j1][i2] = bestH[i1][j1][i2] + j1_lowbound;
                        delete[] bestH[i1][j1][i2];   
                    }
                    if (seq1_P_pairs[j1][i1].valid) {
                        bestP[i1][j1][i2] = bestP[i1][j1][i2] + j1_lowbound;
                        delete[] bestP[i1][j1][i2];
                    }
                    if (seq1_Multi_pairs[j1][i1].valid) {
                        bestMulti[i1][j1][i2] = bestMulti[i1][j1][i2] + j1_lowbound;
                        delete[] bestMulti[i1][j1][i2];
                    }
                }
                if (seq1_M_pairs[j1][i1].valid) {
                    bestM[i1][j1][i2] = bestM[i1][j1][i2] + j1_lowbound;
                    delete[] bestM[i1][j1][i2];
                }
                if (seq1_M2_pairs[j1][i1].valid) {
                    bestM2[i1][j1][i2] = bestM2[i1][j1][i2] + j1_lowbound;
                    delete[] bestM2[i1][j1][i2];
                }
            }
            if (_allowed_pairs[nuci1][nucj1]) {
                if (seq1_H_pairs[j1][i1].valid) {
                    bestH[i1][j1] = bestH[i1][j1] + i1_lowbound; // shift pointer
                    delete[] bestH[i1][j1];
                }
                if (seq1_P_pairs[j1][i1].valid) {
                    bestP[i1][j1] = bestP[i1][j1] + i1_lowbound;
                    delete[] bestP[i1][j1];
                }
                if (seq1_Multi_pairs[j1][i1].valid) {
                    bestMulti[i1][j1] = bestMulti[i1][j1] + i1_lowbound;
                    delete[] bestMulti[i1][j1];
                }
            }
            if (seq1_M_pairs[j1][i1].valid) {
                bestM[i1][j1] = bestM[i1][j1] + i1_lowbound;
                delete[] bestM[i1][j1];
            }
            if (seq1_M2_pairs[j1][i1].valid) {
                bestM2[i1][j1] = bestM2[i1][j1] + i1_lowbound;
                delete[] bestM2[i1][j1];
            }
        }
        
        bestH[i1] = bestH[i1] + i1;
        delete[] bestH[i1];
        bestP[i1] = bestP[i1] + i1;
        delete[] bestP[i1];
        bestMulti[i1] = bestMulti[i1] + i1;
        delete[] bestMulti[i1];
        bestM[i1] = bestM[i1] + i1;
        delete[] bestM[i1];
        bestM2[i1] = bestM2[i1] + i1;
        delete[] bestM2[i1];
    } 
    delete[] bestH;
    delete[] bestP;
    delete[] bestMulti;
    delete[] bestM;
    delete[] bestM2;
    // cout << "done." << endl;

    for (int j1 = 0; j1 < seq1_len; j1++) {
        short low_bound = hmmalign.low_bounds[j1];
        bestC[j1] = bestC[j1] + low_bound;
        delete[] bestC[j1];
    }
    bestC[seq1_len] = bestC[seq1_len] + seq2_len;
    delete[] bestC;

    // internal loop score
    for (int i1=2; i1<seq1_len-3; i1++) {
        short mini = max(1, i1-MAX_LOOP_LEN-1);
        int nuci = sequences[0].nucs[i1];
        for (int j1 = i1+4; j1<seq1_len-1; j1++) {
            int nucj = sequences[0].nucs[j1];
            if (_allowed_pairs[nuci][nucj]) {
                short minj = j1+1;
                short maxj = min(minj+MAX_LOOP_LEN+1, seq1_len);
                for (int k=mini; k<i1; k++) {
                    seq1_internal[i1][j1][k] = seq1_internal[i1][j1][k] + minj; 
                    delete[] seq1_internal[i1][j1][k];
                }
                seq1_internal[i1][j1] = seq1_internal[i1][j1] + mini;
                delete[] seq1_internal[i1][j1];
            }
        }
        seq1_internal[i1] = seq1_internal[i1] + (i1+4);
        delete[] seq1_internal[i1];
    }
    delete[] seq1_internal;
    // cout << "done." << endl;

    for (int i2=2; i2<seq2_len-3; i2++) {
        int nuci = sequences[1].nucs[i2];
        short mini = max(1, i2-MAX_LOOP_LEN-1);
        for (int j2 = i2+4; j2<seq2_len-1; j2++) {
            int nucj = sequences[1].nucs[j2];
            if (_allowed_pairs[nuci][nucj]) {
                short minj = j2+1;
                short maxj = min(minj+MAX_LOOP_LEN+1, seq2_len);
                for (int k=mini; k<i2; k++) {
                    seq2_internal[i2][j2][k] = seq2_internal[i2][j2][k] + minj;
                    delete[] seq2_internal[i2][j2][k];
                }
                seq2_internal[i2][j2] = seq2_internal[i2][j2] + mini;
                delete[] seq2_internal[i2][j2];    
            }   
        }
        seq2_internal[i2] = seq2_internal[i2] + (i2+4); 
        delete[] seq2_internal[i2];
    }
    delete[] seq2_internal;
    // cout << "done." << endl;

    // cout << "clear allocated space in alignment" << endl;
    hmmalign.clear(true);
}

SankoffParser::SankoffParser(float aln_weight, int beam_size, int LFbeam, int LAbeam, float energy_diff, bool is_verbose)
    :weight(aln_weight),
     beam(beam_size),
     lfbeam(LFbeam),
     alnbeam(LAbeam),
     max_energy_diff(energy_diff),
     verbose(is_verbose){

    cout << "beam : " << beam << " lfbeam: " << lfbeam << " alnbeam: " << alnbeam << " max_energy_diff: " << max_energy_diff << endl;  
    
    initialize();
}
