/*
 *HMMAlign.cpp*
 The main code for HMMAlign

 author: Sizhen Li
 created by: 08/2021
*/

#include <string>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <string.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <set>

#include "HMMAlign.h"

using namespace std;

void BeamAlign::traceback(vector<char> &aln1, vector<char> &aln2, HMMManner endmanner){
    int i_seq1 = seq1_len - 1;
    int i_seq2 = seq2_len - 1;
    int step = i_seq1 + i_seq2;
    HMMManner cur_manner = endmanner;

    while (true) {
        if (i_seq1 == 0 && i_seq2 == 0) break;
        // cout << "traceback, cur_manner: " << i_seq1 << " " << i_seq2 << " " << cur_manner << endl;
        switch (cur_manner) {
            case 3: // ALIGN_ALN:
                aln1.push_back(GET_NUC_HMM(seq1[start1+i_seq1]));
                aln2.push_back(GET_NUC_HMM(seq2[start2+i_seq2]));
                cur_manner = bestALN[step][i_seq2].pre;
                step = step - 2;
                i_seq1 --;
                i_seq2 --;
                break;
            case 1: // ALIGN_INS1:
                aln1.push_back(GET_NUC_HMM(seq1[start1+i_seq1]));
                aln2.push_back('-');
                cur_manner = bestINS1[step][i_seq2].pre;
                step = step - 1;
                i_seq1 --;
                break;
            case 2: //ALIGN_INS2:
                aln1.push_back('-');
                aln2.push_back(GET_NUC_HMM(seq2[start2+i_seq2]));
                cur_manner = bestINS2[step][i_seq2].pre;
                step = step - 1;
                i_seq2 --;
                break;
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d: manner %d\n", i_seq1, i_seq2, cur_manner); fflush(stdout);
                assert(false);
        }
    }
}

void BeamAlign::traceback2(int i1, int j1, int i2, int j2, vector<char> &aln1, vector<char> &aln2, HMMManner endmanner){
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1; 
    
    int i_seq1 = seq1len - 1;
    int i_seq2 = seq2len - 1;
    int step = i_seq1 + i_seq2;
    HMMManner cur_manner = endmanner;

    while (true) {
        if (i_seq1 == 0 && i_seq2 == 0) break;
        // cout << "traceback, cur_manner: " << i_seq1 << " " << i_seq2 << " " << cur_manner << endl;
        switch (cur_manner) {
            case 3: // ALIGN_ALN:
                aln1.push_back(GET_NUC_HMM(seq1[start1+i_seq1]));
                aln2.push_back(GET_NUC_HMM(seq2[start2+i_seq2]));
                cur_manner = local_scores[step][i_seq2].alnobj.pre;
                step = step - 2;
                i_seq1 --;
                i_seq2 --;
                break;
            case 1: // ALIGN_INS1:
                aln1.push_back(GET_NUC_HMM(seq1[start1+i_seq1]));
                aln2.push_back('-');
                cur_manner = local_scores[step][i_seq2].ins1obj.pre;
                step = step - 1;
                i_seq1 --;
                break;
            case 2: //ALIGN_INS2:
                aln1.push_back('-');
                aln2.push_back(GET_NUC_HMM(seq2[start2+i_seq2]));
                cur_manner = local_scores[step][i_seq2].ins2obj.pre;
                step = step - 1;
                i_seq2 --;
                break;
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d: manner %d\n", i_seq1, i_seq2, cur_manner); fflush(stdout);
                assert(false);
        }
    }
}

float BeamAlign::get_aln_similarity(char* &seq1_aln_line, char* &seq2_aln_line, char gap_symbol){
    int n_match_pos = 0;
	for(int cnt = 0; cnt < (int)strlen(seq1_aln_line); cnt++)
	{
		if(seq1_aln_line[cnt] != gap_symbol && seq1_aln_line[cnt] != 'N' && seq1_aln_line[cnt] == seq2_aln_line[cnt])
		{
			n_match_pos++;
		}
	}

	int n_total_pos = 0;
	for(int cnt = 0; cnt < (int)strlen(seq1_aln_line); cnt++)
	{
		if(seq1_aln_line[cnt] == gap_symbol && seq2_aln_line[cnt] == gap_symbol)
		{
		}
		else
		{
			n_total_pos++;
		}
	}

	return( (float)n_match_pos / n_total_pos );
}

void BeamAlign::load_init_params(){
    // Copy transition matrix.
	trans_probs = (float**)malloc(sizeof(float*) * (N_STATES + 2));
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		trans_probs[cnt1] = (float*)malloc(sizeof(float) * (N_STATES + 2));
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			trans_probs[cnt1][cnt2] = xlog(0.0f);
		} // cnt2 loop 
	} // cnt1 loop

	// Copy emission probabilities.
	emission_probs = (float**)malloc(sizeof(float*) * (N_OUTPUTS + 2));
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		emission_probs[cnt1] = (float*)malloc(sizeof(float) * (N_STATES + 2));
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			emission_probs[cnt1][cnt2] = xlog(0.0f);
		} // cnt2 loop
	} // cnt1 loop

	fam_hmm_pars = (float*)malloc(sizeof(float) * (N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES + 2));
	fam_thresholds = (float*)malloc(sizeof(float) * (N_BINZ + 2));
}
int BeamAlign::get_bin_index(float similarity, int n_bins){
    if(similarity == 1.0)
	{
		return(n_bins - 1);
	}
	else
	{
		return((int)(n_bins * similarity));
	}
}
void BeamAlign::set_parameters_by_sim(float similarity)
{
    load_init_params(); // allocal space

    // load parameters from file
    char *data_dir = getcwd(NULL, 0);
    char* phmm_pars_file = (char*)malloc(sizeof(char*) * (strlen(data_dir) + strlen("src/data_tables/fam_hmm_pars.dat") + 2));
    sprintf(phmm_pars_file, "%s/%s", data_dir, "src/data_tables/fam_hmm_pars.dat");

    // Read the parameters file. This parameters file contains 10 different parameter sets and thresholds.
	FILE* fam_par_fp = fopen(phmm_pars_file, "r");
    if(fam_par_fp == NULL)
	{
		double* p=0;
		*p = 0;
		printf("Cannot find phmm parameters file, exiting @ %s(%d).\n", __FILE__, __LINE__);
		exit(0);
	}

    // Load all parameters from file.
	for(int cnt = 0; cnt < N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES; cnt++)
	{   
		int tmp = fscanf(fam_par_fp, "%f", &fam_hmm_pars[cnt]);
		// printf("%d: %.3f\n", cnt, fam_hmm_pars[cnt]);
	}
	// printf("\n\n\nReading thresholds!!!\n");
	// Read thresholds.
	for(int cnt = 0; cnt < N_BINZ; cnt++)
	{
		int tmp = fscanf(fam_par_fp, "%f", &fam_thresholds[cnt]);
		// printf("%d: %f\n", cnt, fam_thresholds[cnt]);
	}
	fclose(fam_par_fp);
    free(phmm_pars_file);
    // cout << "test" << endl;

	//int fam_par_set_index = (int)(similarity * (double)N_BINZ);
	int fam_par_set_index = get_bin_index(similarity, N_BINZ);
	// Load emission probabilities.
	// Each parameter set is (N_STATES + N_OUTPUTS) * N_STATES doubles long.
	int start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * fam_par_set_index;
	float* par_ptr = fam_hmm_pars + start_linear_index;
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			//emission_probs[cnt1][cnt2] = *(par_ptr + cnt1 * N_STATES + cnt2);
			emission_probs[cnt1][cnt2] = xlog(par_ptr[cnt1 * N_STATES + cnt2]);
            // cout << cnt1 << " " << cnt2 << " " << emission_probs[cnt1][cnt2] << endl;
		}
	}
	start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * fam_par_set_index + N_STATES * N_OUTPUTS;
	par_ptr = fam_hmm_pars + start_linear_index;
	// Load trans probabilities.
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			//trans_probs[cnt1][cnt2] = *(par_ptr + cnt1 * N_STATES + cnt2);
			trans_probs[cnt1][cnt2] = xlog(par_ptr[cnt1 * N_STATES + cnt2]);
            // cout << cnt1 << " " << cnt2 << " " << trans_probs[cnt1][cnt2] << endl;
		}
    }

    // free(fam_hmm_pars);
    // free(fam_thresholds);
}

float BeamAlign::get_trans_emit_prob0(int prev_state, int current_state, int i, int k, bool new_pars){
    prev_state--;
    current_state--;

    // cout << prev_state << " " << current_state << " " << i << " " << " " << k << endl;

    if (prev_state < 0 || current_state < 0) return xlog(0.0); // N.B.

    float trans_prob; 
    if (new_pars)
        trans_prob = trans_probs[prev_state][current_state];
    else 
        trans_prob = xlog(ML_trans_probs[prev_state][current_state]);
        
    // get_emit_prob
    int i_sym = seq1[i];
    int k_sym = seq2[k];

    // N character
    if (i_sym == 4) i_sym = 0; // rand() % 4;
    if (k_sym == 4) k_sym = 0; // rand() % 4;

    // Fix symbols to gaps in case of of insertions.
    // Gap is coded into value 4 in the emission table.
    if(current_state == 0 || k == 0)  k_sym = 4; 
    // else k_sym = seq2[k];

    if(current_state == 1 || i == 0) i_sym = 4;
    // else i_sym = seq1[i];

    // Compute the symbol index into emission table using the coded nucleotide values:
    // A->0, C->1, G->2, U->3, T->3, .->4
    // This defines a counting system in base of 5. (25 values.)
    // There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
    int sym_index = i_sym * 5 + k_sym;

    // The indices correspond to the end symbol
	if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}

    float emit_prob;
    if (new_pars)
        emit_prob = emission_probs[sym_index][current_state];
    else
        emit_prob = xlog(ML_emit_probs[sym_index][current_state]);


    // cout << prev_state << " " << current_state << " " << i << " " << " " << k << " " << emit_prob << " " << trans_prob << endl;

    return(xlog_mul(emit_prob, trans_prob));
}

// float BeamAlign::get_trans_emit_prob1(int current_state, int next_state, int i, int k, bool new_pars){
//     current_state--;
//     next_state--;

//     // cout << prev_state << " " << current_state << " " << i << " " << " " << k << endl;

//     if (current_state < 0 || next_state < 0) return xlog(0.0); // N.B.

//     float trans_prob; 
//     if (new_pars)
//         trans_prob = trans_probs[current_state][next_state];
//     else 
//         trans_prob = xlog(ML_trans_probs[current_state][next_state]);
        
//     // get_emit_prob
//     int i_sym = seq1[i];
//     int k_sym = seq2[k];

//     // N character
//     if (i_sym == 4) i_sym = 0; // rand() % 4;
//     if (k_sym == 4) k_sym = 0; // rand() % 4;

//     // Fix symbols to gaps in case of of insertions.
//     // Gap is coded into value 4 in the emission table.
//     if(current_state == 0 || k == 0)  k_sym = 4; 
//     // else k_sym = seq2[k];

//     if(current_state == 1 || i == 0) i_sym = 4;
//     // else i_sym = seq1[i];

//     // Compute the symbol index into emission table using the coded nucleotide values:
//     // A->0, C->1, G->2, U->3, T->3, .->4
//     // This defines a counting system in base of 5. (25 values.)
//     // There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
//     int sym_index = i_sym * 5 + k_sym;

//     // The indices correspond to the end symbol?
// 	// if(i == seq1_len && k == seq2_len)
// 	// {
// 	// 	sym_index = 26;
// 	// }

//     float emit_prob;
//     if (new_pars)
//         emit_prob = emission_probs[sym_index][current_state];
//     else
//         emit_prob = xlog(ML_emit_probs[sym_index][current_state]);


//     // cout << prev_state << " " << current_state << " " << i << " " << " " << k << " " << emit_prob << " " << trans_prob << endl;

//     return(xlog_mul(emit_prob, trans_prob));
// }

float BeamAlign::get_trans_emit_prob(int prev_state, int current_state, int i, int k, bool new_pars){
    prev_state--;
    current_state--;

    // cout << prev_state << " " << current_state << " " << i << " " << " " << k << endl;

    if (prev_state < 0 || current_state < 0) return xlog(0.0); // N.B.

    float trans_prob; 
    if (new_pars)
        trans_prob = trans_probs[prev_state][current_state];
    else 
        trans_prob = xlog(ML_trans_probs[prev_state][current_state]);

    // get_emit_prob
    int i_sym = seq1[i];
    int k_sym = seq2[k];

    // N character
    if (i_sym == 4) i_sym = 0; // rand() % 4;
    if (k_sym == 4) k_sym = 0; // rand() % 4;

    // Fix symbols to gaps in case of of insertions.
    // Gap is coded into value 4 in the emission table.
    if(current_state == 0 || k == 0)  k_sym = 4; 
    // else k_sym = seq2[k];

    if(current_state == 1 || i == 0) i_sym = 4;
    // else i_sym = seq1[i];

    // Compute the symbol index into emission table using the coded nucleotide values:
    // A->0, C->1, G->2, U->3, T->3, .->4
    // This defines a counting system in base of 5. (25 values.)
    // There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
    int sym_index = i_sym * 5 + k_sym;

    // The indices correspond to the end symbol?
	if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}

    float emit_prob;
    if (new_pars)
        emit_prob = emission_probs[sym_index][current_state];
    else
        emit_prob = xlog(ML_emit_probs[sym_index][current_state]);


    // cout << prev_state << " " << current_state << " " << i << " " << " " << k << " " << emit_prob << " " << trans_prob << endl;

    return(xlog_mul(emit_prob, trans_prob));
}

float BeamAlign::get_trans_emit_prob(int prev_state, int current_state, int i_1, int k_1, int i, int k, HMMManner s1, HMMManner s2, bool new_pars){
    prev_state--;
    current_state--;

    // cout << "seq1_len, seq2_len: " << seq1_len << " " << seq2_len << endl;
    // cout << prev_state << " " << current_state << " " << i << " " << k << " " << i+start1 << " " << k+start2 <<endl;

    if (prev_state < 0 || current_state < 0) return xlog(0.0); // N.B.

    float trans_prob; 
    if (new_pars)
        trans_prob = trans_probs[prev_state][current_state];
    else 
        trans_prob = xlog(ML_trans_probs[prev_state][current_state]);

    if (i_1 == 0 && k_1 == 0){ 
        if (s1 != HMMMANNER_NONE && current_state != s1 - 1) 
            return xlog(0.0); // if specify the first state, current_state must be s1
        // if (s1 != HMMMANNER_NONE)
        trans_prob = xlog(1.0); // do not include transition prob from prev_state
    }

    if (i+start1 == seq1_len && k+start2 == seq2_len){
        return xlog(1.0); // do not include transition prob from prev_state
    } 
    if (i+start1 == seq1_len || k+start2 == seq2_len) { // TODO: debug
        return xlog(0.0);
    }
    
    int i_sym = seq1[i+start1];
    int k_sym = seq2[k+start2];

    // N character
    if (i_sym == 4) i_sym = 0; // rand() % 4;
    if (k_sym == 4) k_sym = 0; // rand() % 4;

    // Fix symbols to gaps in case of of insertions.
    // Gap is coded into value 4 in the emission table.
    if(current_state == 0)  k_sym = 4; 
    // else k_sym = seq2[k+start2];

    if(current_state == 1) i_sym = 4;
    // else i_sym = seq1[i+start1];

    // Compute the symbol index into emission table using the coded nucleotide values:
    // A->0, C->1, G->2, U->3, T->3, .->4
    // This defines a counting system in base of 5. (25 values.)
    // There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
    int sym_index = i_sym * 5 + k_sym;

    if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}

    float emit_prob;
    if (new_pars)
        emit_prob = emission_probs[sym_index][current_state];
    else
        emit_prob = xlog(ML_emit_probs[sym_index][current_state]);

    return(xlog_mul(emit_prob, trans_prob));
}

float BeamAlign::get_trans_emit_prob_left(int current_state, int next_state, int i, int k, HMMManner s1, HMMManner s2, bool new_pars){
    current_state--;
    next_state--;
    
    // cout << "seq1_len, seq2_len: " << seq1_len << " " << seq2_len << endl;
    // cout << current_state << " " << next_state << " " << i << " " << " " << k << endl;

    if (current_state < 0 || next_state < 0) return xlog(0.0); // N.B.

    float trans_prob; 
    if (new_pars)
        trans_prob = trans_probs[current_state][next_state];
    else 
        trans_prob = xlog(ML_trans_probs[current_state][next_state]);

    if (i == 0 && k == 0){ 
        if (s1 != HMMMANNER_NONE && next_state != s1 - 1) 
            return xlog(0.0); // if specify the first state, current_state must be s1
        // if (s1 != HMMMANNER_NONE)
        return xlog(1.0); // do not include transition prob from prev_state
    }

    if (i+start1 == seq1_len && k+start2 == seq2_len){
        return xlog(1.0); // do not include transition prob from prev_state
    }
    if (i+start1 == seq1_len || k+start2 == seq2_len) { // TODO: debug
        return xlog(0.0);
    }
    
    int i_sym = seq1[i+start1];
    int k_sym = seq2[k+start2];

    // N character
    if (i_sym == 4) i_sym = 0; // rand() % 4;
    if (k_sym == 4) k_sym = 0; // rand() % 4;

    // Fix symbols to gaps in case of of insertions.
    // Gap is coded into value 4 in the emission table.
    if(current_state == 0)  k_sym = 4; 
    // else k_sym = seq2[k+start2];

    if(current_state == 1) i_sym = 4;
    // else i_sym = seq1[i+start1];

    // Compute the symbol index into emission table using the coded nucleotide values:
    // A->0, C->1, G->2, U->3, T->3, .->4
    // This defines a counting system in base of 5. (25 values.)
    // There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
    int sym_index = i_sym * 5 + k_sym;

    if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}

    float emit_prob;
    if (new_pars)
        emit_prob = emission_probs[sym_index][current_state];
    else
        emit_prob = xlog(ML_emit_probs[sym_index][current_state]);

    return(xlog_mul(emit_prob, trans_prob));
}

float BeamAlign::get_trans_emit_prob_right(int prev_state, int current_state, int i_1, int k_1, int i, int k, HMMManner s1, HMMManner s2, bool new_pars){
    prev_state--;
    current_state--;
    
    // cout << "seq1_len, seq2_len: " << seq1_len << " " << seq2_len << endl;
    // cout << prev_state << " " << current_state << " " << i << " " << " " << k << endl;

    if (prev_state < 0 || current_state < 0) return xlog(0.0); // N.B.

    float trans_prob; 
    if (new_pars)
        trans_prob = trans_probs[prev_state][current_state];
    else 
        trans_prob = xlog(ML_trans_probs[prev_state][current_state]);

    if (i_1 == 0 && k_1 == 0){ 
        if (s1 != HMMMANNER_NONE && current_state != s1 - 1) 
            return xlog(0.0); // if specify the first state, current_state must be s1
        
        // trans_prob = xlog(1.0); // do not include transition prob from prev_state
        return xlog(1.0);
    }

    if (i+start1 == seq1_len && k+start2 == seq2_len){
        return xlog(1.0); // do not include transition prob from prev_state
    }
    if (i+start1 == seq1_len || k+start2 == seq2_len) { // TODO: debug
        return xlog(0.0);
    }
    
    int i_sym = seq1[i+start1];
    int k_sym = seq2[k+start2];

    // N character
    if (i_sym == 4) i_sym = 0; // rand() % 4;
    if (k_sym == 4) k_sym = 0; // rand() % 4;

    // Fix symbols to gaps in case of of insertions.
    // Gap is coded into value 4 in the emission table.
    if(current_state == 0)  k_sym = 4; 
    // else k_sym = seq2[k+start2];

    if(current_state == 1) i_sym = 4;
    // else i_sym = seq1[i+start1];

    // Compute the symbol index into emission table using the coded nucleotide values:
    // A->0, C->1, G->2, U->3, T->3, .->4
    // This defines a counting system in base of 5. (25 values.)
    // There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
    int sym_index = i_sym * 5 + k_sym;

    if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}

    float emit_prob;
    if (new_pars)
        emit_prob = emission_probs[sym_index][current_state];
    else
        emit_prob = xlog(ML_emit_probs[sym_index][current_state]);

    // if (i == seq1_len && k == seq2_len)
    // cout << prev_state << " " << current_state << " " << start1+i << " " << start2+k << " " << emit_prob << " " << trans_prob << " " << xlog_mul(emit_prob, trans_prob) << endl;

    return(xlog_mul(emit_prob, trans_prob));
}

bool BeamAlign::update_if_better(AlignState &state, float newscore, HMMManner pre_manner, int step, int i, int k, HMMManner start_manner) {
    if (state.ml < newscore){
        // cout << "better and update: " << pre_manner << " " << step << " " << i << " " << k  << " " <<  state.ml << " " << newscore << endl;
        // state.set(newscore, pre_manner, manner, step, i, k, start_manner);
        state.set(newscore, pre_manner, step, i, k, start_manner);

        // cout << state.step << " " << state.i << " " << state.k << " " << state.ml << endl;
        return true;
    }
    return false;
};

bool BeamAlign::update_if_better_backward(AlignState &state, float newscore) {
    if (state.beta < newscore){
        // cout << "better and update: " <<  state.beta << " " << newscore << endl;
        state.set_backward(newscore);

        return true;
    }
    return false;
};

void BeamAlign::update(AlignState &state, float newscore, HMMManner pre_manner, int step, int i, int k) {
    // cout << "update: " << pre_manner << " " << i << " " << k  << " " <<  state.ml << " " << newscore << endl;
    state.set(newscore, pre_manner, step, i, k);
};

void BeamAlign::beam_prune(unordered_map<int, AlignState> &beamstep){
    scores.clear();
    for (auto &item : beamstep) {
        int ik = item.first;
        AlignState &cand = item.second;
        scores.push_back(make_pair(cand.ml, ik));
    }
    if (scores.size() <= beam) return;
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    // return threshold;
}

void BeamAlign::prepare(int j1, int j2) {
    // bestALN.clear();
    // bestALN.resize(j1 + j2 + 2);
    // bestINS1.clear();
    // bestINS1.resize(j1 + j2 + 2);
    // bestINS2.clear();
    // bestINS2.resize(j1 + j2 + 2);

    int sum_len = j1+j2+2;
    bestALN = new unordered_map<int, AlignState>[sum_len];
    bestINS1 = new unordered_map<int, AlignState>[sum_len];
    bestINS2 = new unordered_map<int, AlignState>[sum_len];
}

void BeamAlign::cal_align_prob(float threshold){
    float aln_prob, ins1_prob, ins2_prob;
   
    for (int s = 0; s < seq1_len + seq2_len; s++){
        for(auto &item : bestALN[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            aln_prob = xlog_div(xlog_mul(state.ml, state.beta), forward_score);
            if (aln_prob > float(-9.91152)) 
                aln_env[i][k].prob = xlog_sum(aln_env[i][k].prob, aln_prob);
        }

        for(auto &item : bestINS1[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            ins1_prob = xlog_div(xlog_mul(state.ml, state.beta), forward_score);
            if (ins1_prob > float(-9.91152)) 
                aln_env[i][k].prob = xlog_sum(aln_env[i][k].prob, ins1_prob);
        }

        for(auto &item : bestINS2[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            ins2_prob = xlog_div(xlog_mul(state.ml, state.beta), forward_score);
            if (ins2_prob > float(-9.91152)) 
                aln_env[i][k].prob = xlog_sum(aln_env[i][k].prob, ins2_prob);
        }
    }
    // cout << "after calculating probs." << endl;

    vector<int> cands;
    up_bounds.resize(seq1_len+1, -1);
    low_bounds.resize(seq1_len+1, seq2_len);

    max_j1.resize(seq1_len + seq2_len, -1);
    min_j1.resize(seq1_len + seq2_len, seq1_len);

    max_range = 0;
    int sum_range = 0;
    for (int i = 0; i < seq1_len; i++){
        int upbound = -1;
        int lowbound = seq2_len;
        cands.clear();
        for(auto &item : aln_env[i]){
            int k = item.first;
            float prob = item.second.prob;
            // cout << i << " " << k << " " << prob  << " " << threshold<< endl;
            if (prob < threshold) {
                cands.push_back(k);
                continue;
            }
            
            lowbound = min(k, lowbound);
            upbound = max(k, upbound);

            // runtime debug
            // lowbound = i;
            // upbound = i;
        }

        if (lowbound > upbound) {
            low_bounds[i] = 0;
            up_bounds[i] = seq2_len - 1;
        } else {
            low_bounds[i] = lowbound;
            up_bounds[i] = upbound;
        }

        // cout << i << " " << low_bounds[i] << " " << up_bounds[i] << endl;
        for (int k = low_bounds[i]; k <= up_bounds[i]; k++) {
            max_j1[i+k] = max(max_j1[i+k], i);
            min_j1[i+k] = min(min_j1[i+k], i);
        }

        max_range = max(max_range, up_bounds[i] - low_bounds[i] + 1);
        sum_range += up_bounds[i] - low_bounds[i] + 1;

        for (auto &k : cands) {
            aln_env[i].erase(k);
        }
    }

    // for (int s = 1; s < seq1_len + seq2_len; s++){
    //     cout << s << " " << min_j1[s] << " " << max_j1[s] << endl;
    // }

    cout << "max_range: " << max_range << endl;
    cout << "average range: " << sum_range / float(seq1_len) << endl;

    // boundary case
    low_bounds[seq1_len] = seq2_len;
    up_bounds[seq1_len] = seq2_len;
}

void BeamAlign::forward() {
    // seq1_len = seq1.size();
    // seq2_len = seq2.size();

    prepare(seq1_len, seq2_len);

    // initial state
    bestALN[0][0].ml = xlog(1.0);
    bestALN[0][0].step = 0;
    bestALN[0][0].i = 0;
    bestALN[0][0].k = 0;

    float trans_emit_prob, newscore;
    for(int s = 0; s < seq1_len + seq2_len; ++s){
        
        unordered_map<int, AlignState>& beamALN = bestALN[s];
        unordered_map<int, AlignState>& beamINS1 = bestINS1[s];
        unordered_map<int, AlignState>& beamINS2 = bestINS2[s];

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};
        for (int pre_m=0; pre_m < 3; pre_m++){
            // unordered_map<keys, State, key_hash, key_equal>& beamH = *beams[i_beams];
            unordered_map<int, AlignState>& beamstep = *beams[pre_m];

            // beam prune
            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);

            HMMManner manner;
            switch (pre_m) {
                case 2:
                    manner = MANNER_ALN;
                    break;
                case 1:
                    manner = MANNER_INS2;
                    break;
                case 0:
                    manner = MANNER_INS1;
                    break;
                default:
                    manner = HMMMANNER_NONE;
                    break;
            } // end swtich pre_m
            
            for (auto &item : beamstep) {
                AlignState state = item.second;
                int i = state.i;
                int k = state.k;
                float ml = state.ml;
                if (i + k != s) continue;

                // cout << "s, i, k: " << s << " " << i << " " << k << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = MANNER_ALN; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step = s + 2;
                        // if ((next_i > seq1_len) || (next_k > seq2_len)) continue;

                        if ((next_i < seq1_len && next_k < seq2_len) || (next_i == seq1_len && next_k == seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, true);
                            newscore = xlog_sum(bestALN[next_step][next_k].ml, xlog_mul(state.ml, trans_emit_prob));
                            update(bestALN[next_step][next_k], newscore, manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = MANNER_INS1; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step = s + 1;
                        if ((next_i >= seq1_len) || (next_k >= seq2_len)) continue;

                        trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, true);
                        newscore = xlog_sum(bestINS1[next_step][next_k].ml, xlog_mul(state.ml, trans_emit_prob));
                        update(bestINS1[next_step][next_k], newscore, manner, next_step, next_i, next_k);

                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = MANNER_INS2; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step = s + 1;
                        if ((next_i > seq1_len) || (next_k >= seq2_len)) continue;

                        trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, true);
                        newscore = xlog_sum(bestINS2[next_step][next_k].ml, xlog_mul(state.ml, trans_emit_prob));
                        update(bestINS2[next_step][next_k], newscore, manner, next_step, next_i, next_k);

                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }

    }

    forward_score = bestALN[seq1_len+seq2_len][seq2_len].ml;
    printf("forward score: %.10f\n", forward_score);
}

void BeamAlign::backward() {
    bestALN[seq1_len + seq2_len][seq2_len].beta = xlog(1.0);
    // bestALN[seq1_len + seq2_len][seq2_len].manner = 3; //ALIGN_ALN;
    bestALN[seq1_len + seq2_len][seq2_len].step = seq1_len + seq2_len;
    bestALN[seq1_len + seq2_len][seq2_len].i = seq1_len;
    bestALN[seq1_len + seq2_len][seq2_len].k = seq2_len;

    for(int s = seq1_len + seq2_len - 2; s >= 0; --s) {
        float trans_emit_prob;
        unordered_map<int, AlignState> &beamINS1 = bestINS1[s];
        unordered_map<int, AlignState> &beamINS2 = bestINS2[s];
        unordered_map<int, AlignState> &beamALN = bestALN[s];
        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN}; // N.B.

        for (int pre_m=0; pre_m < beams.size(); pre_m++){
            unordered_map<int, AlignState> &beamstep = *beams[pre_m]; // N.B.
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                int i = state.i;
                int k = state.k;
                int step = state.step;
                // int manner = state.manner;
                if (i + k != s) continue;

                HMMManner manner;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        break;
                    case 1:
                        manner = MANNER_INS2;
                        break;
                    case 0:
                        manner = MANNER_INS1;
                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                int next_manner;
                int next_i, next_k, next_step;
                for (int m = 1; m <= 3; m++){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;

                        if (((next_i == seq1_len) && (next_k == seq2_len)) || ((next_i < seq1_len) && (next_k < seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, true);
                            state.beta = xlog_sum(state.beta, xlog_mul(bestALN[next_step][next_k].beta, trans_emit_prob));
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, true);
                            state.beta = xlog_sum(state.beta, xlog_mul(bestINS1[next_step][next_k].beta, trans_emit_prob));
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, true);
                            state.beta = xlog_sum(state.beta, xlog_mul(bestINS2[next_step][next_k].beta, trans_emit_prob));
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    float back_score = bestALN[0][0].beta;
    printf("backward score: %.10f\n", back_score);
    // return back_score;
}

float BeamAlign::viterbi_path(bool newpara) {
    seq1_len = seq1.size();
    seq2_len = seq2.size();

    prepare(seq1_len, seq2_len);

    // initial state
    bestALN[0][0].ml = xlog(1.0);
    bestALN[0][0].pre = MANNER_ALN;
    bestALN[0][0].step = 0;
    bestALN[0][0].i = 0;
    bestALN[0][0].k = 0;

    start1=0;
    start2=0;

    float trans_emit_prob;
    // processMem_t mem;
    for(int s = 0; s < seq1_len + seq2_len; ++s){

        unordered_map<int, AlignState>& beamALN = bestALN[s];
        unordered_map<int, AlignState>& beamINS1 = bestINS1[s];
        unordered_map<int, AlignState>& beamINS2 = bestINS2[s];

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};
        for (int pre_m=0; pre_m < 3; pre_m++){
            // unordered_map<keys, State, key_hash, key_equal>& beamH = *beams[i_beams];
            unordered_map<int, AlignState>& beamstep = *beams[pre_m];

            // beam prune
            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);

            HMMManner manner;
            switch (pre_m) {
                case 2:
                    manner = MANNER_ALN;
                    break;
                case 1:
                    manner = MANNER_INS2;
                    break;
                case 0:
                    manner = MANNER_INS1;
                    break;
                default:
                    manner = HMMMANNER_NONE;
                    break;
            } // end swtich pre_m
            
            for (auto &item : beamstep) {
                AlignState state = item.second;
                int i = state.i;
                int k = state.k;
                float ml = state.ml;
                if (i + k != s) continue; //TODO: 0, 0 debug

                // cout << "s, i, k: " << s << " " << i << " " << k << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = MANNER_ALN; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step = s + 2;
                        // if ((next_i > seq1_len) || (next_k > seq2_len)) continue;

                        if ((next_i < seq1_len && next_k < seq2_len) || (next_i == seq1_len && next_k == seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                            // cout << ml << " " << trans_emit_prob << " " << bestALN[next_step][next_k].ml << endl;
                            update_if_better(bestALN[next_step][next_k], xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = MANNER_INS1; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step = s + 1;
                        if ((next_i >= seq1_len) || (next_k >= seq2_len)) continue;

                        trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                        update_if_better(bestINS1[next_step][next_k], xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k);
                        
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = MANNER_INS2; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step = s + 1;
                        if ((next_i >= seq1_len) || (next_k >= seq2_len)) continue;

                        trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                        update_if_better(bestINS2[next_step][next_k], xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k);
                        
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }

    }
    // cout << "ml: " << bestALN[seq1_len + seq2_len][seq2_len].ml << endl;
    float ml = bestALN[seq1_len + seq2_len][seq2_len].ml;
    cout << "viterbi path score: " << ml << endl;

    if (newpara) {
        // get similarity again
        vector<char> aln1, aln2;
        traceback(aln1, aln2, bestALN[seq1_len + seq2_len][seq2_len].pre);

        int aln_len = aln1.size();
        assert (aln1.size() == aln2.size());

        char* seq1_aln_line_str = (char*)malloc(sizeof(char) * (aln_len+2)); 
        char* seq2_aln_line_str = (char*)malloc(sizeof(char) * (aln_len+2));
        for (int i = aln_len - 1; i >= 0; --i){
            seq1_aln_line_str[aln_len - 1 - i] = aln1[i];
            seq2_aln_line_str[aln_len - 1 - i] = aln2[i];
        }
        seq1_aln_line_str[aln_len] = 0;
        seq2_aln_line_str[aln_len] = 0;
        cout << seq1_aln_line_str << endl;
        cout << seq2_aln_line_str << endl;

        float similarity = get_aln_similarity(seq1_aln_line_str, seq2_aln_line_str);
        cout << "similarity with new paras: " << similarity << endl;

        aln1.clear();
        aln2.clear();
        free(seq1_aln_line_str);
        free(seq2_aln_line_str);

        // backward
        viterbi_backward(newpara);

        return ml;
    }

    // traceback
    vector<char> aln1, aln2;
    traceback(aln1, aln2, bestALN[seq1_len + seq2_len][seq2_len].pre);

    int aln_len = aln1.size();
    assert (aln1.size() == aln2.size());

    char* seq1_aln_line_str = (char*)malloc(sizeof(char) * (aln_len+2)); 
    char* seq2_aln_line_str = (char*)malloc(sizeof(char) * (aln_len+2));
    for (int i = aln_len - 1; i >= 0; --i){
        seq1_aln_line_str[aln_len - 1 - i] = aln1[i];
        seq2_aln_line_str[aln_len - 1 - i] = aln2[i];
    }
    seq1_aln_line_str[aln_len] = 0;
    seq2_aln_line_str[aln_len] = 0;
    cout << seq1_aln_line_str << endl;
    cout << seq2_aln_line_str << endl;

    float similarity = get_aln_similarity(seq1_aln_line_str, seq2_aln_line_str);
    cout << "similarity with initial paras: " << similarity << endl;

    aln1.clear();
    aln2.clear();
    free(seq1_aln_line_str);
    free(seq2_aln_line_str);

    // load new parameters
    set_parameters_by_sim(similarity);

    forward();
    backward();

    float threshold = fam_thresholds[get_bin_index(similarity, N_BINZ)];
    cout << "threshold: " << threshold << endl;

    aln_env.clear();
    aln_env.resize(seq1_len + 1);
    cal_align_prob(threshold); // prune
    // cout << "after calculate align probs." << endl;

    delete[] bestALN;
    delete[] bestINS1;
    delete[] bestINS2;

    return ml;
}

float BeamAlign::evalulate(bool newpara, const std::set<pair<int, int>> &align_pairs) {
    seq1_len = seq1.size();
    seq2_len = seq2.size();

    prepare(seq1_len, seq2_len);

    // initial state
    bestALN[0][0].ml = xlog(1.0);
    bestALN[0][0].step = 0;
    bestALN[0][0].i = 0;
    bestALN[0][0].k = 0;

    start1=0;
    start2=0;

    float trans_emit_prob;
    // processMem_t mem;
    for(int s = 0; s < seq1_len + seq2_len; ++s){

        unordered_map<int, AlignState>& beamALN = bestALN[s];
        unordered_map<int, AlignState>& beamINS1 = bestINS1[s];
        unordered_map<int, AlignState>& beamINS2 = bestINS2[s];

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};
        for (int pre_m=0; pre_m < 3; pre_m++){
            // unordered_map<keys, State, key_hash, key_equal>& beamH = *beams[i_beams];
            unordered_map<int, AlignState>& beamstep = *beams[pre_m];

            // beam prune
            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);

            HMMManner manner;
            switch (pre_m) {
                case 2:
                    manner = MANNER_ALN;
                    break;
                case 1:
                    manner = MANNER_INS2;
                    break;
                case 0:
                    manner = MANNER_INS1;
                    break;
                default:
                    manner = HMMMANNER_NONE;
                    break;
            } // end swtich pre_m
            
            for (auto &item : beamstep) {
                AlignState state = item.second;
                int i = state.i;
                int k = state.k;
                float ml = state.ml;
                if (i + k != s) continue; //TODO: 0, 0 debug

                // cout << "s, i, k: " << s << " " << i << " " << k << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = MANNER_ALN; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step = s + 2;
                        // if ((next_i > seq1_len) || (next_k > seq2_len)) continue;

                        // constrain
                        if (next_i < seq1_len && next_k < seq2_len){
                            if (align_pairs.find(make_pair(next_i, next_k)) == align_pairs.end()) continue;
                        }

                        if ((next_i < seq1_len && next_k < seq2_len) || (next_i == seq1_len && next_k == seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                            // cout << ml << " " << trans_emit_prob << " " << bestALN[next_step][next_k].ml << endl;
                            update_if_better(bestALN[next_step][next_k], xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = MANNER_INS1; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step = s + 1;
                        if ((next_i >= seq1_len) || (next_k >= seq2_len)) continue;

                        // constrain
                        if (align_pairs.find(make_pair(next_i, -1)) == align_pairs.end()) continue;

                        trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                        update_if_better(bestINS1[next_step][next_k], xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k);
                        
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = MANNER_INS2; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step = s + 1;
                        if ((next_i >= seq1_len) || (next_k >= seq2_len)) continue;

                        // constrain
                        if (align_pairs.find(make_pair(-1, next_k)) == align_pairs.end()) continue;

                        trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                        update_if_better(bestINS2[next_step][next_k], xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k);
                        
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }

    }
    // cout << "ml: " << bestALN[seq1_len + seq2_len][seq2_len].ml << endl;
    float ml = bestALN[seq1_len + seq2_len][seq2_len].ml;
    cout << "viterbi path score: " << ml << endl;

    return ml;
}

void BeamAlign::viterbi_backward(bool newpara) {
    bestALN[seq1_len + seq2_len][seq2_len].beta = xlog(1.0);
    bestALN[seq1_len + seq2_len][seq2_len].step = seq1_len + seq2_len;
    bestALN[seq1_len + seq2_len][seq2_len].i = seq1_len;
    bestALN[seq1_len + seq2_len][seq2_len].k = seq2_len;

    for(int s = seq1_len + seq2_len - 2; s >= 0; --s) {
        float trans_emit_prob;
        unordered_map<int, AlignState> &beamINS1 = bestINS1[s];
        unordered_map<int, AlignState> &beamINS2 = bestINS2[s];
        unordered_map<int, AlignState> &beamALN = bestALN[s];
        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN}; // N.B.

        for (int pre_m=0; pre_m < beams.size(); pre_m++){
            unordered_map<int, AlignState> &beamstep = *beams[pre_m]; // N.B.
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                int i = state.i;
                int k = state.k;
                int step = state.step;
                // int manner = state.manner;
                if (i + k != s) continue;

                HMMManner manner;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        break;
                    case 1:
                        manner = MANNER_INS2;
                        break;
                    case 0:
                        manner = MANNER_INS1;
                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                int next_manner;
                int next_i, next_k, next_step;
                for (int m = 1; m <= 3; m++){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;

                        if (next_i < seq1_len && next_k < seq2_len) {
                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                        } 

                        if (((next_i == seq1_len) && (next_k == seq2_len)) || ((next_i < seq1_len) && (next_k < seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                            update_if_better_backward(state, xlog_mul(bestALN[next_step][next_k].beta, trans_emit_prob));
                            // state.beta = xlog_sum(state.beta, xlog_mul(bestALN[next_step][next_k].beta, trans_emit_prob));

                            // debug
                            // cout << pre_m << " " << i << " " << k << " " << state.ml << " " << state.beta << " " << xlog_mul(state.ml, state.beta) << endl;
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;

                        // check whether in valid region
                        // loop the region slightly to make sure arrive at the end position with state s2
                        if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                        if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                            update_if_better_backward(state, xlog_mul(bestINS1[next_step][next_k].beta, trans_emit_prob));
                            // state.beta = xlog_sum(state.beta, xlog_mul(bestINS1[next_step][next_k].beta, trans_emit_prob));
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;

                        // check whether in valid region
                        // loop the region slightly to make sure arrive at the end position with state s2
                        if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                        if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, next_i, next_k, newpara);
                            update_if_better_backward(state, xlog_mul(bestINS2[next_step][next_k].beta, trans_emit_prob));
                            // state.beta = xlog_sum(state.beta, xlog_mul(bestINS2[next_step][next_k].beta, trans_emit_prob));
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    float back_score = bestALN[0][0].beta;
    printf("backward score: %.10f\n", back_score);
    // return back_score;
}

float BeamAlign::viterbi_path_local(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose) {    
    // from zero (start point): include first transition prob ALN->ALN/INS1/INS2
    // partial: don't include first transition prob; set i1-- i2-- as start point
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1;
    
    local_scores.clear();
    local_scores.resize(seq1len + seq2len + 2);

    start1 = i1;
    start2 = i2;
    // end1 = j1; //  + 1;
    // end2 = j2; //  + 1;

    // initial state
    local_scores[0][0].i = 0;
    local_scores[0][0].k = 0;
    local_scores[0][0].alnobj.ml = xlog(1.0);
    
    float trans_emit_prob;
    for(int s = 0; s < seq1len + seq2len; ++s){
        unordered_map<int, AlignState3>& beamstep = local_scores[s];
        for (auto &item : beamstep) {
            AlignState3 state3 = item.second;
            int i = state3.i;
            int k = state3.k;
            if (i + k != s) continue; // lisiz TODO: DEBUG

            for (int pre_m=0; pre_m < 3; pre_m++){
                HMMManner manner, start_manner;
                float ml;
                AlignState state;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        state = state3.alnobj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_ALN;
                        else start_manner = state.start_manner;

                        break;
                    case 1:
                        manner = MANNER_INS2;
                        state = state3.ins2obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS2;
                        else start_manner = state.start_manner;

                        break;
                    case 0:
                        manner = MANNER_INS1;
                        state = state3.ins1obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS1;
                        else start_manner = state.start_manner;

                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                // check 
                if (manner == HMMMANNER_NONE) continue; // wrong manner
                ml = state.ml;
                if (ml <= xlog(0)) continue; // invalid state

                // if (s1 == HMMMANNER_NONE) 
                // cout << s << " " << i << " " << k << " " << ml << " " << seq1len << " " << seq2len << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                bool update;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                        case 3: // ALIGN_ALN:
                            next_manner = MANNER_ALN; // ALIGN_ALN;
                            next_i = i + 1;
                            next_k = k + 1;
                            next_step = s + 2;

                            // if (next_i < seq1_len && next_k < seq2_len) {
                            //     // check whether in valid region
                            //     // loop the region slightly to make sure arrive at the end position with state s2
                            //     if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            //     if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            // } 
                            
                            if ((next_i < seq1len && next_k < seq2len) || (next_i == seq1len && next_k == seq2len)) { // seq2_len may out of reach of seq1_len
                                trans_emit_prob = get_trans_emit_prob(manner, next_manner, i, k, next_i, next_k, s1, HMMMANNER_NONE, true);
                                update = update_if_better(local_scores[next_step][next_k].alnobj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                                if (update) {
                                    local_scores[next_step][next_k].i = next_i;
                                    local_scores[next_step][next_k].k = next_k;
                                }
                            }

                            // cout << m << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] << " " << up_bounds[start1 + next_i] << endl;
                            break;

                        case 1: // ALIGN_INS1:
                            next_manner = MANNER_INS1; // ALIGN_INS1;
                            next_i = i + 1;
                            next_k = k;
                            next_step = s + 1;
                            
                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            // if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            // if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, i, k, next_i, next_k, s1, HMMMANNER_NONE, true);
                            update = update_if_better(local_scores[next_step][next_k].ins1obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            break;

                        case 2: // ALIGN_INS2:
                            next_manner = MANNER_INS2; // ALIGN_INS2;
                            next_i = i;
                            next_k = k + 1;
                            next_step = s + 1;

                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            // if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            // if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            
                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, i, k, next_i, next_k, s1, HMMMANNER_NONE, true);
                            update = update_if_better(local_scores[next_step][next_k].ins2obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            break;
                        
                        default:
                            break;
                    }
                }
            }
        }
    }
    // cout << "end." << endl;

    switch (s2)
    {
        case 1:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].ins1obj.ml;
            break;

        case 2:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].ins2obj.ml;
            break;

        case 3:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].alnobj.ml;
            break;
        
        default:
            return local_scores[seq1len + seq2len][seq2len].alnobj.ml; // TODO
            break;
    }

    
}

void BeamAlign::viterbi_path_local22(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose) {  
    int org_i1 = i1;
    int org_i2 = i2;

    // from zero (start point): include first transition prob ALN->ALN/INS1/INS2
    // partial: don't include first transition prob; set i1-- i2-- as start point
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1;

    start1 = i1;
    start2 = i2;
    // end1 = j1; //  + 1;
    // end2 = j2; //  + 1;

    // initial state
    local_scores.clear();
    local_scores.resize(seq1len + seq2len + 2);

    local_scores[0][0].i = 0;
    local_scores[0][0].k = 0;
    local_scores[0][0].alnobj.ml = xlog(1.0);
    
    float trans_emit_prob;
    for(int s = 0; s < seq1len + seq2len; ++s){
        // cout << "s: " << s << " " << i1 << " " << i2 <<  endl;
        unordered_map<int, AlignState3>& beamstep = local_scores[s];
        for (auto &item : beamstep) {
            AlignState3 state3 = item.second;
            int i = state3.i;
            int k = state3.k;
            if (i + k != s) continue; // lisiz TODO: DEBUG

            for (int pre_m=0; pre_m < 3; pre_m++){
                HMMManner manner, start_manner;
                float ml;
                AlignState state;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        state = state3.alnobj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_ALN;
                        else start_manner = state.start_manner;

                        break;
                    case 1:
                        manner = MANNER_INS2;
                        state = state3.ins2obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS2;
                        else start_manner = state.start_manner;

                        break;
                    case 0:
                        manner = MANNER_INS1;
                        state = state3.ins1obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS1;
                        else start_manner = state.start_manner;

                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                // check 
                if (manner == HMMMANNER_NONE) continue; // wrong manner
                ml = state.ml;
                if (ml <= xlog(0)) continue; // invalid state

                // if (s1 == HMMMANNER_NONE) 
                // cout << s << " " << start1 << " " << start2 << " " << i << " " << k << " " << start1 + i << " " << start2 + k << " " << seq1len << " " << seq2len << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                bool update;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                        case 3: // ALIGN_ALN:
                            next_manner = MANNER_ALN; // ALIGN_ALN;
                            next_i = i + 1;
                            next_k = k + 1;
                            next_step = s + 2;

                            if (next_i < seq1len && next_k < seq2len) {
                                // check whether in valid region
                                // loop the region slightly to make sure arrive at the end position with state s2
                                if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                                if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            } 
                            
                            if ((next_i < seq1len && next_k < seq2len) || (next_i == seq1len && next_k == seq2len)) {
                                trans_emit_prob = get_trans_emit_prob(manner, next_manner, i, k, next_i, next_k, s1, HMMMANNER_NONE, true);
                                update = update_if_better(local_scores[next_step][next_k].alnobj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                                if (update) {
                                    local_scores[next_step][next_k].i = next_i;
                                    local_scores[next_step][next_k].k = next_k;
                                }
                            }

                            // cout << m << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] << " " << up_bounds[start1 + next_i] << endl;
                            break;

                        case 1: // ALIGN_INS1:
                            next_manner = MANNER_INS1; // ALIGN_INS1;
                            next_i = i + 1;
                            next_k = k;
                            next_step = s + 1;
                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, i, k, next_i, next_k, s1, HMMMANNER_NONE, true);
                            update = update_if_better(local_scores[next_step][next_k].ins1obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            break;

                        case 2: // ALIGN_INS2:
                            next_manner = MANNER_INS2; // ALIGN_INS2;
                            next_i = i;
                            next_k = k + 1;
                            next_step = s + 1;

                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            
                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob(manner, next_manner, i, k, next_i, next_k, s1, HMMMANNER_NONE, true);
                            update = update_if_better(local_scores[next_step][next_k].ins2obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            break;
                        
                        default:
                            break;
                    }
                }

                // save alignment score
                if (start2+k >= low_bounds[start1+i] && start2+k <= up_bounds[start1+i]) {
                    if (pre_m + 1 == s1) {
                        // cout << s1 << " " << org_i1 << " " << org_i2 << " " << i << " " << k << endl;
                        // all_local_scores[s1-1][start1+1][start2+1-low_bounds[start1+1]][i][k] = ml;
                        all_local_scores[s1-1][org_i1][org_i2][start1+i][start2+k] = ml;
                    }
                }
                
            }
        }
    }
}

float BeamAlign::viterbi_path_local_left(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose) {
    // from zero (start point): include first transition prob ALN->ALN/INS1/INS2
    // partial: don't include first transition prob; set i1-- i2-- as start point
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1;
    
    local_scores.clear();
    local_scores.resize(seq1len + seq2len + 2);

    start1 = i1;
    start2 = i2;
    // end1 = j1; //  + 1;
    // end2 = j2; //  + 1;

    // initial state
    local_scores[0][0].i = 0;
    local_scores[0][0].k = 0;
    local_scores[0][0].alnobj.ml = xlog(1.0);
    
    float trans_emit_prob;
    for(int s = 0; s < seq1len + seq2len; ++s){
        unordered_map<int, AlignState3>& beamstep = local_scores[s];
        for (auto &item : beamstep) {
            AlignState3 state3 = item.second;
            int i = state3.i;
            int k = state3.k;
            if (i + k != s) continue; // lisiz TODO: DEBUG

            for (int pre_m=0; pre_m < 3; pre_m++){
                HMMManner manner, start_manner;
                float ml;
                AlignState state;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        state = state3.alnobj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_ALN;
                        else start_manner = state.start_manner;

                        break;
                    case 1:
                        manner = MANNER_INS2;
                        state = state3.ins2obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS2;
                        else start_manner = state.start_manner;

                        break;
                    case 0:
                        manner = MANNER_INS1;
                        state = state3.ins1obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS1;
                        else start_manner = state.start_manner;

                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                // check 
                if (manner == HMMMANNER_NONE) continue; // wrong manner
                ml = state.ml;
                if (ml <= xlog(0)) continue; // invalid state

                // cout << s << " " << i << " " << k << " " << start1 + i << " " << start2 + k << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                bool update;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                        case 3: // ALIGN_ALN:
                            next_manner = MANNER_ALN; // ALIGN_ALN;
                            next_i = i + 1;
                            next_k = k + 1;
                            next_step = s + 2;

                            // if (next_i < seq1_len && next_k < seq2_len) {
                            //     // check whether in valid region
                            //     // loop the region slightly to make sure arrive at the end position with state s2
                            //     if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            //     if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            // } 
                            
                            if ((next_i < seq1len && next_k < seq2len) || (next_i == seq1len && next_k == seq2len)) { // seq2_len may out of reach of seq1_len
                                trans_emit_prob = get_trans_emit_prob_left(manner, next_manner, i, k, s1, s2, true);
                                // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                                update = update_if_better(local_scores[next_step][next_k].alnobj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                                if (update) { // lisiz TODO: better solution
                                    local_scores[next_step][next_k].i = next_i;
                                    local_scores[next_step][next_k].k = next_k;
                                }
                            }

                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 1: // ALIGN_INS1:
                            next_manner = MANNER_INS1; // ALIGN_INS1;
                            next_i = i + 1;
                            next_k = k;
                            next_step = s + 1;
                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            // cout << m << " "   << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                        
                            // if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            // if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_left(manner, next_manner, i, k, s1, s2, true);
                            // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                            update = update_if_better(local_scores[next_step][next_k].ins1obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 2: // ALIGN_INS2:
                            next_manner = MANNER_INS2; // ALIGN_INS2;
                            next_i = i;
                            next_k = k + 1;
                            next_step = s + 1;

                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            // if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            // if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            
                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_left(manner, next_manner, i, k, s1, s2, true);
                            update = update_if_better(local_scores[next_step][next_k].ins2obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;
                        
                        default:
                            break;
                    }
                }
            }
        }
    }

    // traceback for debug
    // vector<char> aln1, aln2;
    // cout << local_scores[seq1_len + seq2_len][seq2_len].alnobj.pre << endl;
    // if (local_scores[seq1_len + seq2_len][seq2_len].alnobj.pre != HMMMANNER_NONE){
    //     traceback2(aln1, aln2, local_scores[seq1_len + seq2_len][seq2_len].alnobj.pre);
    //     reverse(aln1.begin(), aln1.end());
    //     reverse(aln2.begin(), aln2.end());

    //     string seq1aln1(aln1.begin(), aln1.end());
    //     string seq2aln1(aln2.begin(), aln2.end());

    //     cout << seq1aln1 << " " << seq2aln1 << endl;

    //     aln1.clear();
    //     aln2.clear();
    // }
    
    switch (s2)
    {
        case 1:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].ins1obj.ml;
            break;

        case 2:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].ins2obj.ml;
            break;

        case 3:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].alnobj.ml;
            break;
        
        default:
            return local_scores[seq1len + seq2len][seq2len].alnobj.ml; // TODO
            break;
    }

    
}


void BeamAlign::viterbi_path_local_left22(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose) {
    int org_i1 = i1;
    int org_i2 = i2;

    // from zero (start point): include first transition prob ALN->ALN/INS1/INS2
    // partial: don't include first transition prob; set i1-- i2-- as start point
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1;
    
    // prepare(seq1_len, seq2_len);
    local_scores.clear();
    local_scores.resize(seq1len + seq2len + 2);

    start1 = i1;
    start2 = i2;
    // end1 = j1; //  + 1;
    // end2 = j2; //  + 1;

    // initial state
    local_scores[0][0].i = 0;
    local_scores[0][0].k = 0;
    local_scores[0][0].alnobj.ml = xlog(1.0);
    
    float trans_emit_prob;
    for(int s = 0; s < seq1len + seq2len; ++s){
        unordered_map<int, AlignState3>& beamstep = local_scores[s];
        for (auto &item : beamstep) {
            AlignState3 state3 = item.second;
            int i = state3.i;
            int k = state3.k;
            if (i + k != s) continue; // lisiz TODO: DEBUG

            for (int pre_m=0; pre_m < 3; pre_m++){
                HMMManner manner, start_manner;
                float ml;
                AlignState state;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        state = state3.alnobj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_ALN;
                        else start_manner = state.start_manner;

                        break;
                    case 1:
                        manner = MANNER_INS2;
                        state = state3.ins2obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS2;
                        else start_manner = state.start_manner;

                        break;
                    case 0:
                        manner = MANNER_INS1;
                        state = state3.ins1obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS1;
                        else start_manner = state.start_manner;

                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                // check 
                if (manner == HMMMANNER_NONE) continue; // wrong manner
                ml = state.ml;
                if (ml <= xlog(0)) continue; // invalid state
                if ((i >= seq1len) || (k >= seq2len)) continue; // boundary case

                // cout << s << " " << i << " " << k << " " << start1 + i << " " << start2 + k << " " << seq1_len << " " << seq2_len << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                bool update;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                        case 3: // ALIGN_ALN:
                            next_manner = MANNER_ALN; // ALIGN_ALN;
                            next_i = i + 1;
                            next_k = k + 1;
                            next_step = s + 2;

                            if (next_i < seq1len && next_k < seq2len) {
                                // check whether in valid region
                                // loop the region slightly to make sure arrive at the end position with state s2
                                if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                                if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            } 
                            
                            if ((next_i < seq1len && next_k < seq2len) || (next_i == seq1len && next_k == seq2len)) { // seq2_len may out of reach of seq1_len
                                trans_emit_prob = get_trans_emit_prob_left(manner, next_manner, i, k, s1, s2, true);
                                // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                                update = update_if_better(local_scores[next_step][next_k].alnobj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                                if (update) { // lisiz TODO: better solution
                                    local_scores[next_step][next_k].i = next_i;
                                    local_scores[next_step][next_k].k = next_k;
                                }
                            }

                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 1: // ALIGN_INS1:
                            next_manner = MANNER_INS1; // ALIGN_INS1;
                            next_i = i + 1;
                            next_k = k;
                            next_step = s + 1;
                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            // cout << m << " "   << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                        
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_left(manner, next_manner, i, k, s1, s2, true);
                            // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                            update = update_if_better(local_scores[next_step][next_k].ins1obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 2: // ALIGN_INS2:
                            next_manner = MANNER_INS2; // ALIGN_INS2;
                            next_i = i;
                            next_k = k + 1;
                            next_step = s + 1;

                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            
                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_left(manner, next_manner, i, k, s1, s2, true);
                            update = update_if_better(local_scores[next_step][next_k].ins2obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;
                        
                        default:
                            break;
                    }
                }
                
                // save alignment score
                // cout << s1-1 << " " << pre_m << " " << start1+1 << " " << start2+1-low_bounds[start1+1] << endl;
                if (start2+k >= low_bounds[start1+i] && start2+k <= up_bounds[start1+i]) {
                    // left_local_scores[s1-1][pre_m][start1+1][start2+1-low_bounds[start1+1]][i][k] = ml;
                    // cout << s1 << " " << s2 << " " << org_i1 << " " << org_i2 << " " << i << " " << k << endl;
                    left_local_scores[s1-1][pre_m][org_i1][org_i2][start1+i][start2+k] = ml;
                }
            }
        }
    }

}

float BeamAlign::viterbi_path_local_right(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool limited, bool verbose) {
    // cout << i1 << " " << j1 << " " << i2 << " " << j2 << " " << s1 << " " << s2 <<  endl;

    // from zero (start point): include first transition prob ALN->ALN/INS1/INS2
    // partial: don't include first transition prob; set i1-- i2-- as start point
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1;
    
    local_scores.clear();
    local_scores.resize(seq1len + seq2len + 2);

    start1 = i1;
    start2 = i2;
    // end1 = j1; //  + 1;
    // end2 = j2; //  + 1;

    // initial state
    local_scores[0][0].i = 0;
    local_scores[0][0].k = 0;
    local_scores[0][0].alnobj.ml = xlog(1.0);
    
    float trans_emit_prob;
    for(int s = 0; s < seq1len + seq2len; ++s){
        unordered_map<int, AlignState3>& beamstep = local_scores[s];
        for (auto &item : beamstep) {
            AlignState3 state3 = item.second;
            int i = state3.i;
            int k = state3.k;
            if (i + k != s) continue; // lisiz TODO: DEBUG

            for (int pre_m=0; pre_m < 3; pre_m++){
                HMMManner manner, start_manner;
                float ml;
                AlignState state;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        state = state3.alnobj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_ALN;
                        else start_manner = state.start_manner;

                        break;
                    case 1:
                        manner = MANNER_INS2;
                        state = state3.ins2obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS2;
                        else start_manner = state.start_manner;

                        break;
                    case 0:
                        manner = MANNER_INS1;
                        state = state3.ins1obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS1;
                        else start_manner = state.start_manner;

                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                // check 
                if (manner == HMMMANNER_NONE) continue; // wrong manner
                ml = state.ml;
                if (ml <= xlog(0)) continue; // invalid state

                // cout << s << " " << i << " " << k << " " << start1 + i << " " << start2 + k << endl;

                HMMManner next_manner;
                int next_i, next_k, next_step;
                bool update;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                        case 3: // ALIGN_ALN:
                            next_manner = MANNER_ALN; // ALIGN_ALN;
                            next_i = i + 1;
                            next_k = k + 1;
                            next_step = s + 2;

                            // check whether in valid region
                            if (limited && next_i < seq1len && next_k < seq2len) {
                                if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                                if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            }
                            
                            if ((next_i < seq1len && next_k < seq2len) || (next_i == seq1len && next_k == seq2len)) {
                                trans_emit_prob = get_trans_emit_prob_right(manner, next_manner, i, k, next_i, next_k, s1, s2, true);
                                // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                                update = update_if_better(local_scores[next_step][next_k].alnobj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                                if (update) { // lisiz TODO: better solution
                                    local_scores[next_step][next_k].i = next_i;
                                    local_scores[next_step][next_k].k = next_k;
                                }
                            }

                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 1: // ALIGN_INS1:
                            next_manner = MANNER_INS1; // ALIGN_INS1;
                            next_i = i + 1;
                            next_k = k;
                            next_step = s + 1;
                            
                            // check whether in valid region
                            if (limited) {
                                if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                                if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            }

                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_right(manner, next_manner, i, k, next_i, next_k, s1, s2, true);
                            // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                            update = update_if_better(local_scores[next_step][next_k].ins1obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 2: // ALIGN_INS2:
                            next_manner = MANNER_INS2; // ALIGN_INS2;
                            next_i = i;
                            next_k = k + 1;
                            next_step = s + 1;

                            // check whether in valid region
                            if (limited) {
                                if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                                if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            }      
                            
                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_right(manner, next_manner, i, k, next_i, next_k, s1, s2, true);
                            update = update_if_better(local_scores[next_step][next_k].ins2obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;
                        
                        default:
                            break;
                    }
                }
            }
        }
    }

    // traceback for debug
    // vector<char> aln1, aln2;
    // cout << local_scores[seq1_len + seq2_len][seq2_len].alnobj.pre << endl;
    // if (local_scores[seq1_len + seq2_len][seq2_len].alnobj.pre != HMMMANNER_NONE){
    //     traceback2(aln1, aln2, local_scores[seq1_len + seq2_len][seq2_len].alnobj.pre);
    //     reverse(aln1.begin(), aln1.end());
    //     reverse(aln2.begin(), aln2.end());
    //     string seq1aln1(aln1.begin(), aln1.end());
    //     string seq2aln1(aln2.begin(), aln2.end());
    //     cout << seq1aln1 << " " << seq2aln1 << endl;
    //     aln1.clear();
    //     aln2.clear();
    // }
    

    switch (s2)
    {
        case 1:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].ins1obj.ml;
            break;

        case 2:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].ins2obj.ml;
            break;

        case 3:
            return local_scores[seq1len + seq2len - 2][seq2len - 1].alnobj.ml;
            break;
        
        default:
            return local_scores[seq1len + seq2len][seq2len].alnobj.ml;
            break;
    }  
}

void BeamAlign::viterbi_path_local_right22(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose) {
    int org_i1 = i1;
    int org_i2 = i2;

    // from zero (start point): include first transition prob ALN->ALN/INS1/INS2
    // partial: don't include first transition prob; set i1-- i2-- as start point
    if (i1>0) i1--; 
    if (i2>0) i2--;

    // prepare
    int seq1len = j1 - i1 + 1;
    int seq2len = j2 - i2 + 1;
    
    local_scores.clear();
    local_scores.resize(seq1len + seq2len + 2);

    start1 = i1;
    start2 = i2;
    // end1 = j1; //  + 1;
    // end2 = j2; //  + 1;

    // initial state
    local_scores[0][0].i = 0;
    local_scores[0][0].k = 0;
    local_scores[0][0].alnobj.ml = xlog(1.0);
    
    float trans_emit_prob;
    for(int s = 0; s < seq1len + seq2len; ++s){
        unordered_map<int, AlignState3>& beamstep = local_scores[s];
        for (auto &item : beamstep) {
            AlignState3 state3 = item.second;
            int i = state3.i;
            int k = state3.k;
            if (i + k != s) continue; // lisiz TODO: DEBUG

            for (int pre_m=0; pre_m < 3; pre_m++){
                HMMManner manner, start_manner;
                float ml;
                AlignState state;
                switch (pre_m) {
                    case 2:
                        manner = MANNER_ALN;
                        state = state3.alnobj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_ALN;
                        else start_manner = state.start_manner;

                        break;
                    case 1:
                        manner = MANNER_INS2;
                        state = state3.ins2obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS2;
                        else start_manner = state.start_manner;

                        break;
                    case 0:
                        manner = MANNER_INS1;
                        state = state3.ins1obj;
                        // start manner
                        if ((s + i) == 0) start_manner = MANNER_INS1;
                        else start_manner = state.start_manner;

                        break;
                    default:
                        manner = HMMMANNER_NONE;
                        break;
                } // end swtich pre_m

                // check 
                if (manner == HMMMANNER_NONE) continue; // wrong manner
                ml = state.ml;
                if (ml <= xlog(0)) continue; // invalid state
                if ((i >= seq1len) || (k >= seq2len)) continue; // boundary case

                HMMManner next_manner;
                int next_i, next_k, next_step;
                bool update;
                for (int m = 3; m >= 1; m--){
                    switch (m)
                    {
                        case 3: // ALIGN_ALN:
                            next_manner = MANNER_ALN; // ALIGN_ALN;
                            next_i = i + 1;
                            next_k = k + 1;
                            next_step = s + 2;

                            if (next_i < seq1len && next_k < seq2len) {
                                // check whether in valid region
                                // loop the region slightly to make sure arrive at the end position with state s2
                                if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                                if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            } 
                            
                            if ((next_i < seq1len && next_k < seq2len) || (next_i == seq1len && next_k == seq2len)) { // seq2_len may out of reach of seq1_len
                                trans_emit_prob = get_trans_emit_prob_right(manner, next_manner, i, k, next_i, next_k, s1, s2, true);
                                // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                                update = update_if_better(local_scores[next_step][next_k].alnobj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                                if (update) { // lisiz TODO: better solution
                                    local_scores[next_step][next_k].i = next_i;
                                    local_scores[next_step][next_k].k = next_k;
                                }
                            }

                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 1: // ALIGN_INS1:
                            next_manner = MANNER_INS1; // ALIGN_INS1;
                            next_i = i + 1;
                            next_k = k;
                            next_step = s + 1;
                            
                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            // cout << m << " "   << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 

                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_right(manner, next_manner, i, k, next_i, next_k, s1, s2, true);
                            // cout << trans_emit_prob << " " << local_scores[next_step][next_k].alnobj.ml << " " << xlog_mul(ml, trans_emit_prob) << endl;
                            update = update_if_better(local_scores[next_step][next_k].ins1obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;

                        case 2: // ALIGN_INS2:
                            next_manner = MANNER_INS2; // ALIGN_INS2;
                            next_i = i;
                            next_k = k + 1;
                            next_step = s + 1;

                            // check whether in valid region
                            // loop the region slightly to make sure arrive at the end position with state s2
                            if ((start2 + next_k) < low_bounds[start1 + next_i]) continue;
                            if ((start2 + next_k) > up_bounds[start1 + next_i]) continue; 
                            
                            // boundary case
                            if ((next_i >= seq1len) || (next_k >= seq2len)) continue;

                            trans_emit_prob = get_trans_emit_prob_right(manner, next_manner, i, k, next_i, next_k, s1, s2, true);
                            update = update_if_better(local_scores[next_step][next_k].ins2obj, xlog_mul(ml, trans_emit_prob), manner, next_step, next_i, next_k, start_manner);
                            if (update) {
                                local_scores[next_step][next_k].i = next_i;
                                local_scores[next_step][next_k].k = next_k;
                            }
                            
                            // cout << m << " " << update << " "  << next_i << " " << next_k << " " << start1 + next_i << " " <<  start2 + next_k << " " << low_bounds[start1 + next_i] - 1 << " " << up_bounds[start1 + next_i] + 1 << endl;
                            break;
                        
                        default:
                            break;
                    }
                }

                // save alignment score
                // cout << s1-1 << " " << pre_m << " " << start1+1 << " " << start2+1-low_bounds[start1+1] << endl;
                if (start2+k >= low_bounds[start1+i] && start2+k <= up_bounds[start1+i]) {
                    // right_local_scores[s1-1][pre_m][start1+1][start2+1-low_bounds[start1+1]][i][k] = ml;
                    // cout << s1-1 << " " << pre_m << " " << org_i1 << " " << org_i2 << " " << start1+i << " " << start2+k << endl;
                    right_local_scores[s1-1][pre_m][org_i1][org_i2][start1+i][start2+k] = ml;
                }
            }
        }
    }

}


void BeamAlign::viterbi_path_all_locals()
{
    cout << "viterbi_path_all_locals start" << endl;
    // int seq1len = seq1_len;
    // int seq2len = seq2_len;

    all_local_scores = new float****[3];
    left_local_scores = new float*****[3];
    right_local_scores = new float*****[3];
    for (int s1 = 0; s1 < 3; s1++){
        all_local_scores[s1] = new float***[seq1_len];
        for (int i1 = 1; i1 < seq1_len; i1++) {
            int i2_low = max(1, low_bounds[i1]);
            int i2_up = up_bounds[i1];
            int i2_range = i2_up - i2_low + 1;
            // cout << "all_local_scores i2_range: " << i2_low << " " << i2_up << " " << i2_range << endl;

            all_local_scores[s1][i1] = new float**[i2_range];
            all_local_scores[s1][i1] = all_local_scores[s1][i1] - i2_low;

            int j1_range = min(35, seq1_len - i1 + 1);
            // cout << "all_local_scores j1_range: " << i1+1 << " " << i1+j1_range << " " << j1_range << endl;
            for (int i2 = i2_low; i2 <= i2_up; i2++) {
                all_local_scores[s1][i1][i2] = new float*[j1_range];
                all_local_scores[s1][i1][i2] = all_local_scores[s1][i1][i2] - (i1-1);

                int j2_range = min(35, seq2_len - i2 + 1); 
                for (int j1 = i1-1; j1 < i1+j1_range-1; j1++) {
                    int j2_low = max(i2-1, low_bounds[j1]);
                    int j2_up = min(i2+j2_range-1, up_bounds[j1]);
                    if (j2_up < j2_low) continue;

                    int new_j2_range = j2_up - j2_low + 1;
                    // cout << "all_local_scores new_j2_range: " << j2_low << " " << j2_up << " " << new_j2_range << endl;

                    all_local_scores[s1][i1][i2][j1] = new float[new_j2_range];
                    all_local_scores[s1][i1][i2][j1] =  all_local_scores[s1][i1][i2][j1] - j2_low;

                    for (int j2 = j2_low; j2 <= j2_up; j2++) {
                        // cout << "all_local_scores: " <<  s1 << " " << i1 << " " << i2  << " " << j1 << " " << j2 << " " << i2_range << " " << j1_range << " " << new_j2_range << endl;

                        all_local_scores[s1][i1][i2][j1][j2] = xlog(0);
                    }
                }
            }
        }
        

        left_local_scores[s1] = new float****[3];
        right_local_scores[s1] = new float****[3];
        for (int s2 = 0; s2 < 3; s2++){
            left_local_scores[s1][s2] = new float***[seq1_len];
            right_local_scores[s1][s2] = new float***[seq1_len];

            for (int i1 = 1; i1 < seq1_len; i1++) {
                int i2_low = max(1, low_bounds[i1]);
                int i2_up = up_bounds[i1];
                int i2_range = i2_up - i2_low + 1;
                // cout << "i2_range: " << i2_low << " " << i2_up << " " << i2_range << endl;

                left_local_scores[s1][s2][i1] = new float**[i2_range];
                right_local_scores[s1][s2][i1] = new float**[i2_range];

                left_local_scores[s1][s2][i1] = left_local_scores[s1][s2][i1] - i2_low;
                right_local_scores[s1][s2][i1] = right_local_scores[s1][s2][i1] - i2_low;
                
                // cout << i1 << " " << low_bounds[i1] <<  " " << up_bounds[i1] << " " << i2_range << endl;

                int j1_range = min(35, seq1_len - i1 + 1);
                // cout << "j1_range: " << i1+1 << " " << i1+j1_range << " " << j1_range << endl;
                for (int i2 = i2_low; i2 <= i2_up; i2++) {
                    left_local_scores[s1][s2][i1][i2] = new float*[j1_range];
                    right_local_scores[s1][s2][i1][i2] = new float*[j1_range];

                    left_local_scores[s1][s2][i1][i2] = left_local_scores[s1][s2][i1][i2] - (i1-1);
                    right_local_scores[s1][s2][i1][i2] = right_local_scores[s1][s2][i1][i2] - (i1-1);
                    
                    int j2_range = min(35, seq2_len - i2 + 1); // TODO, save more, change data structure
                    for (int j1 = i1-1; j1 < i1+j1_range-1; j1++) {
                        int j2_low = max(i2-1, low_bounds[j1]);
                        int j2_up = min(i2+j2_range-1, up_bounds[j1]);
                        if (j2_up < j2_low) continue;
                        int new_j2_range = j2_up - j2_low + 1;
                        
                        left_local_scores[s1][s2][i1][i2][j1] = new float[new_j2_range];
                        right_local_scores[s1][s2][i1][i2][j1] = new float[new_j2_range];

                        left_local_scores[s1][s2][i1][i2][j1] = left_local_scores[s1][s2][i1][i2][j1] - j2_low;
                        right_local_scores[s1][s2][i1][i2][j1] =  right_local_scores[s1][s2][i1][i2][j1] - j2_low;
                        
                        for (int j2 = j2_low; j2 <= j2_up; j2++) {
                            // cout << s1 << " " << s2 << " " << i1 << " " << i2  << " " << j1 << " " << j2 << " " << i2_range << " " << j1_range << " " << new_j2_range << endl;

                            left_local_scores[s1][s2][i1][i2][j1][j2] = xlog(0);
                            right_local_scores[s1][s2][i1][i2][j1][j2] = xlog(0);
                        }
                    }
                }
            }
        }
    }
    
    // pre-compute
    for (int s1 = 0; s1 < 3; s1++){
        HMMManner s1manner;
        switch (s1)
        {
            case 0:
                s1manner = MANNER_INS1;
                break;
            case 1:
                s1manner = MANNER_INS2;
                break;
            case 2:
                s1manner = MANNER_ALN;
                break;
            
            default:
                break;
        }

        for (int i1 = 1; i1 < seq1_len; i1++) {
            for (int i2 = max(1, low_bounds[i1]); i2 <= up_bounds[i1]; i2++){
                // cout << s1 << " " << i1 << " " << i2 << " " << min(seq1_len-1, i1+31) << " " << min(seq2_len-1, i2+31) <<  endl;
                // cout << "test00" << endl;
                viterbi_path_local22(i1, min(seq1_len-1, i1+31), i2, min(seq2_len-1, i2+31), s1manner, HMMMANNER_NONE);
                // cout << "test0" << endl;
                viterbi_path_local_left22(i1, min(seq1_len-1, i1+31), i2, min(seq2_len-1, i2+31), s1manner, HMMMANNER_NONE);
                // cout << "test1" << endl;
                viterbi_path_local_right22(i1, min(seq1_len-1, i1+31), i2, min(seq2_len-1, i2+31), s1manner, HMMMANNER_NONE);
            }
        }
    }

}

void BeamAlign::set(int beam_size, vector<int> &seq1_nuc_types, vector<int> &seq2_nuc_types)
{
    cout << "aln beam size: " << beam_size << endl;
    beam = beam_size;
    seq1 = seq1_nuc_types;
    seq2 = seq2_nuc_types;
}

BeamAlign::BeamAlign(int beam_size,
                     bool is_eval,
                     bool is_verbose)
    : beam(beam_size), 
      eval(is_eval), 
      verbose(is_verbose){}

void BeamAlign::clear(bool local_path){
    for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		free(trans_probs[cnt1]);
	} // cnt1 loop
    free(trans_probs);

	// Copy emission probabilities.
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		free(emission_probs[cnt1]);
	} // cnt1 loop
    free(emission_probs);

    free(fam_hmm_pars);
    free(fam_thresholds);

    if (local_path) {
        for (int s1 = 0; s1 < 3; s1++){
            for (int i1 = 1; i1 < seq1_len; i1++) {
                int i2_low = max(1, low_bounds[i1]);
                int i2_up = up_bounds[i1];
                int i2_range = i2_up - i2_low + 1;

                int j1_range = min(35, seq1_len - i1 + 1);
                for (int i2 = i2_low; i2 <= i2_up; i2++) {
                    int j2_range = min(35, seq2_len - i2 + 1); 
                    for (int j1 = i1-1; j1 < i1+j1_range-1; j1++) {
                        int j2_low = max(i2-1, low_bounds[j1]);
                        int j2_up = min(i2+j2_range-1, up_bounds[j1]);
                        if (j2_up < j2_low) continue;
                        all_local_scores[s1][i1][i2][j1] =  all_local_scores[s1][i1][i2][j1] + j2_low;
                        delete[] all_local_scores[s1][i1][i2][j1];
                    }
                    all_local_scores[s1][i1][i2] = all_local_scores[s1][i1][i2] + (i1-1);
                    delete[] all_local_scores[s1][i1][i2];
                }
                all_local_scores[s1][i1] = all_local_scores[s1][i1] + i2_low;
                delete[] all_local_scores[s1][i1];
            
            }
            delete[] all_local_scores[s1];
        
            for (int s2 = 0; s2 < 3; s2++){
                for (int i1 = 1; i1 < seq1_len; i1++) {
                    int i2_low = max(1, low_bounds[i1]);
                    int i2_up = up_bounds[i1];
                    int i2_range = i2_up - i2_low + 1;

                    int j1_range = min(35, seq1_len - i1 + 1);
                    for (int i2 = i2_low; i2 <= i2_up; i2++) {
                        int j2_range = min(35, seq2_len - i2 + 1); // TODO, save more, change data structure
                        for (int j1 = i1-1; j1 < i1+j1_range-1; j1++) {
                            int j2_low = max(i2-1, low_bounds[j1]);
                            int j2_up = min(i2+j2_range-1, up_bounds[j1]);
                            if (j2_up < j2_low) continue;
                            left_local_scores[s1][s2][i1][i2][j1] = left_local_scores[s1][s2][i1][i2][j1] + j2_low;
                            right_local_scores[s1][s2][i1][i2][j1] =  right_local_scores[s1][s2][i1][i2][j1] + j2_low;
                        
                            delete[] left_local_scores[s1][s2][i1][i2][j1];
                            delete[] right_local_scores[s1][s2][i1][i2][j1];
                        }
                        left_local_scores[s1][s2][i1][i2] = left_local_scores[s1][s2][i1][i2] + (i1-1);
                        right_local_scores[s1][s2][i1][i2] = right_local_scores[s1][s2][i1][i2] + (i1-1);

                        delete[] left_local_scores[s1][s2][i1][i2];
                        delete[] right_local_scores[s1][s2][i1][i2];
                    }
                    left_local_scores[s1][s2][i1] = left_local_scores[s1][s2][i1] + i2_low;
                    right_local_scores[s1][s2][i1] = right_local_scores[s1][s2][i1] + i2_low;
                        
                    delete[] left_local_scores[s1][s2][i1];
                    delete[] right_local_scores[s1][s2][i1];
                }
                delete[] left_local_scores[s1][s2];
                delete[] right_local_scores[s1][s2];
            }
            delete[] left_local_scores[s1];
            delete[] right_local_scores[s1];
        }
        
        delete[] all_local_scores;
        delete[] left_local_scores;
        delete[] right_local_scores;
    }
}