/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h> 

#include "LinearPartition.h"
#include "bpp.cpp"

#define SPECIAL_HP

using namespace std;

// unsigned long quickselect_partition(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper) {
//     pf_type pivot = scores[upper].first;
//     while (lower < upper) {
//         while (scores[lower].first < pivot) ++lower;
//         while (scores[upper].first > pivot) --upper;
//         if (scores[lower].first == scores[upper].first) ++lower;
//         else if (lower < upper) swap(scores[lower], scores[upper]);
//     }
//     return upper;
// }

// // in-place quick-select
// pf_type quickselect(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
//     if ( lower == upper ) return scores[lower].first;
//     unsigned long split = quickselect_partition(scores, lower, upper);
//     unsigned long length = split - lower + 1;
//     if (length == k) return scores[split].first;
//     else if (k  < length) return quickselect(scores, lower, split-1, k);
//     else return quickselect(scores, split+1, upper, k - length);
// }


pf_type BeamCKYPFParser::beam_prune(std::unordered_map<int, PFState> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        PFState &cand = item.second;
        int k = i - 1;
        pf_type newalpha = (k >= 0 ? bestC[k].alpha : pf_type(0.0)) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN_PF;
    pf_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}


void BeamCKYPFParser::prepare(unsigned len) {
    seq_length = len;

    nucs = new int[seq_length];
    bestC = new PFState[seq_length];
    bestH = new unordered_map<int, PFState>[seq_length];
    bestP = new unordered_map<int, PFState>[seq_length];
    bestM = new unordered_map<int, PFState>[seq_length];
    bestM2 = new unordered_map<int, PFState>[seq_length];
    bestMulti = new unordered_map<int, PFState>[seq_length];
    
    scores.reserve(seq_length);
}

void BeamCKYPFParser::postprocess() {

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
}

void BeamCKYPFParser::parse(const string& seq) {
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    vector<int> next_pair[NOTON];
    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            // next_pair
            next_pair[nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = seq_length-1; j >=0; --j) {
                next_pair[nuci][j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }
        }
    }

#ifdef SPECIAL_HP
#ifdef lv
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
#endif

#ifdef lv
        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;
#else
        if(seq_length > 0) Fast_LogPlusEquals(bestC[0].alpha, score_external_unpaired(0, 0));
        if(seq_length > 1) Fast_LogPlusEquals(bestC[1].alpha, score_external_unpaired(0, 1));
#endif

    value_type newscore;
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, PFState>& beamstepH = bestH[j];
        unordered_map<int, PFState>& beamstepMulti = bestMulti[j];
        unordered_map<int, PFState>& beamstepP = bestP[j];
        unordered_map<int, PFState>& beamstepM2 = bestM2[j];
        unordered_map<int, PFState>& beamstepM = bestM[j];
        PFState& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                // for nucj put H(j, j_next) into H[j_next]
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;
#ifdef lv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
#else
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore);
#endif
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
                for (auto &item : beamstepH) {
                    int i = item.first;
                    PFState &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)=
#ifdef lv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-i-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[i];
                        else if (jnext-i-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[i];
                        else if (jnext-i-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[i];
#endif
                        newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
#else
                        newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore);
#endif
                    }

                    // 2. generate p(i, j)
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            for(auto& item : beamstepMulti) {
                int i = item.first;
                PFState& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
#ifdef lv
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
#else
                        newscore = score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lv
                    newscore = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore/kT);
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore);
#endif
                }
            }
        }

        // beam of P
        {   
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            for(auto& item : beamstepP) {
                int i = item.first;
                PFState& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
#ifndef lv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lv
                                newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1);

                                // SHAPE for Vienna only
                                if (use_shape)
                                {
                                    newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                }


                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore/kT);
#else
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);

#endif
                            } else {
                                // single branch
#ifdef lv
                                newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore/kT);
#else
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j,
                                                                       nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);
#endif
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
#ifdef lv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore/kT);
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lv
                    newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = state.alpha + newscore/kT;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = state.alpha + newscore;
#endif
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        PFState& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        PFState& prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lv
                        newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore/kT);      
#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore);
#endif
                    } else {
#ifdef lv
                        newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore/kT);       
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore);
#endif
                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            for(auto& item : beamstepM2) {
                int i = item.first;
                PFState& state = item.second;

                // 1. multi-loop
                for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                    int nucp = nucs[p];
                    int q = next_pair[nucp][j];
                    if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lv
                    Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha);      

#else
                    newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                    Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha + newscore);      
#endif
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
            }
        }

        // beam of M
        {
            if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                PFState& state = item.second;
                if (j < seq_length-1) {
#ifdef lv
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha + newscore); 
#endif
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lv
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
                    
#else
                newscore = score_external_unpaired(j+1, j+1);
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha + newscore); 
#endif
            }
        }
    }  // end of for-loo j

    PFState& viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    // unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;

#ifdef lv
    fprintf(stderr,"Free Energy of Ensemble: %.2f kcal/mol\n", -kT * viterbi.alpha / 100.0);
#else
    fprintf(stderr,"Log Partition Coefficient: %.5f\n", viterbi.alpha);
#endif

    if(is_verbose) fprintf(stderr,"Partition Function Calculation Time: %.2f seconds.\n", parse_elapsed_time);

    fflush(stdout);

    // lhuang
    if(pf_only && !forest_file.empty()) dump_forest(seq, true); // inside-only forest

    if(!pf_only){
        outside(next_pair);
    	if (!forest_file.empty())
    	  dump_forest(seq, false); // inside-outside forest
        
        // cal_PairProb(viterbi);

        if (mea_) PairProb_MEA(seq);

        if (threshknot_) ThreshKnot(seq);
    }
    // postprocess();
    return;
}

void BeamCKYPFParser::print_states(FILE *fptr, unordered_map<int, PFState>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        PFState & state = item.second;
        if (inside_only) fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
            fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

void BeamCKYPFParser::dump_forest(const string & seq, bool inside_only) {  
    printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
    FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
    fprintf(fptr, "%s\n", seq.c_str());
    int n = seq.length(), j;
    for (j = 0; j < n; j++) {
        if (inside_only) fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
        else fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
    }
    double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
    for (j = 0; j < n; j++) 
        print_states(fptr, bestP[j], j, "P", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM[j], j, "M", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

BeamCKYPFParser::BeamCKYPFParser(int beam_size,
                                bool nosharpturn,
                                bool verbose,
                                string bppfile,
                                string bppfileindex,
                                bool pfonly,
                                float bppcutoff,
                                string forestfile,
                                bool mea,
                                float MEA_gamma,
                                string MEA_file_index,
                                bool MEA_bpseq,
                                bool ThreshKnot,
                                float ThreshKnot_threshold,
                                string ThreshKnot_file_index,
                                string shape_file_path)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      bpp_file(bppfile),
      bpp_file_index(bppfileindex),
      pf_only(pfonly),
      bpp_cutoff(bppcutoff),
      forest_file(forestfile), 
      mea_(mea),
      gamma(MEA_gamma),
      mea_file_index(MEA_file_index),
      bpseq(MEA_bpseq),
      threshknot_(ThreshKnot),
      threshknot_threshold(ThreshKnot_threshold),
      threshknot_file_index(ThreshKnot_file_index){
#ifdef lv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif

    if (shape_file_path != "" ){
        use_shape = true;
        int position;
        string data;

        double temp_after_mb_shape;

        ifstream in(shape_file_path);

        if (!in.good()){
            cout<<"Reading SHAPE file error!"<<endl;
            assert(false);
        }

        // actually, we can combine the SHAPE_data and the energy_stack together
        while (!(in >> position >> data).fail()) {
            // cout<<"position data "<< int(position)<<endl<<data<<endl;
            // assert(int(position) == SHAPE_data.size() + 1);
            // cout<<"data "<<data<<endl;
            if (isdigit(int(data[0])) == 0){
                SHAPE_data.push_back(double((-1.000000)));
            }

            else {
                SHAPE_data.push_back(stod(data));
            }
            

        }

        for (int i = 0; i<SHAPE_data.size(); i++){
            temp_after_mb_shape = SHAPE_data[i] < 0 ? 0. : (m * log(SHAPE_data[i] + 1) + b);

            pseudo_energy_stack.push_back((int)roundf(temp_after_mb_shape * 100.));

            assert(pseudo_energy_stack.size() == i + 1 );

            // cout<<"pseudo energy "<<i<<' '<<SHAPE_data[i]<<' '<<temp_after_mb_shape<<' '<<pseudo_energy_stack[i]<<' '<<pseudo_energy_stack.size()<<endl;

        }
    }



}

// int main(int argc, char** argv){

//     struct timeval total_starttime, total_endtime;
//     gettimeofday(&total_starttime, NULL);

//     int beamsize = 100;
//     bool sharpturn = false;
//     bool is_verbose = false;
//     string bpp_file;
//     string bpp_prefix;
//     bool pf_only = false;
//     float bpp_cutoff = 0.0;
//     string forest_file;

//     float MEA_gamma = 3.0;
//     bool mea = false;
//     bool MEA_bpseq = false;
//     string MEA_prefix;
//     float ThreshKnot_threshold = 0.3;
//     bool ThreshKnot = false;
//     string ThresKnot_prefix;


//     // SHAPE
//     string shape_file_path = "";



//     if (argc > 1) {
//         beamsize = atoi(argv[1]);
//         sharpturn = atoi(argv[2]) == 1;
//         is_verbose = atoi(argv[3]) == 1;
//         bpp_file = argv[4];
//         bpp_prefix = argv[5];
//         pf_only = atoi(argv[6]) == 1;
//         bpp_cutoff = atof(argv[7]);
//     	forest_file = argv[8];
//         mea = atoi(argv[9]) == 1;
//         MEA_gamma = atof(argv[10]);
//         ThreshKnot = atoi(argv[11]) == 1;
//         ThreshKnot_threshold = atof(argv[12]);
//         ThresKnot_prefix = argv[13];
//         MEA_prefix = argv[14];
//         MEA_bpseq = atoi(argv[15]) == 1;
//         shape_file_path = argv[16];
//     }


//     if (is_verbose) printf("beam size: %d\n", beamsize);

//     // variables for decoding
//     int num=0, total_len = 0;
//     unsigned long long total_states = 0;
//     double total_score = .0;
//     double total_time = .0;

//     int seq_index = 0;
//     string bpp_file_index = "";
//     string ThreshKnot_file_index = "";
//     string MEA_file_index = "";

//     for (string seq; getline(cin, seq);) {
//         if (seq.length() == 0)
//             continue;

//         if (seq[0] == ';' || seq[0] == '>') {
//             printf("%s\n", seq.c_str());
//             if (!bpp_file.empty()) {
//                 FILE *fptr = fopen(bpp_file.c_str(), "a"); 
//                 if (fptr == NULL) { 
//                     printf("Could not open file!\n"); 
//                     return 0; 
//                 }
//                 fprintf(fptr, "%s\n", seq.c_str());
//                 fclose(fptr); 
//             }
//             continue;
//         }

//         if (!isalpha(seq[0])){
//             printf("Unrecognized sequence: %s\n", seq.c_str());
//             continue;
//         }

//         seq_index ++;
//         if (!bpp_prefix.empty()) bpp_file_index = bpp_prefix + to_string(seq_index);

//         if (!ThresKnot_prefix.empty()) ThreshKnot_file_index = ThresKnot_prefix + to_string(seq_index);

//         if (!MEA_prefix.empty()) MEA_file_index = MEA_prefix + to_string(seq_index);
        
//         // convert to uppercase
//         transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

//         // convert T to U
//         replace(seq.begin(), seq.end(), 'T', 'U');

//         // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
//         BeamCKYPFParser parser(beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);

//         // BeamCKYPFParser::DecoderResult result = parser.parse(seq);
//         parser.parse(seq);
//     }

//     gettimeofday(&total_endtime, NULL);
//     double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

//     if(is_verbose) fprintf(stderr,"Total Time: %.2f seconds.\n", total_elapsed_time);

//     return 0;
// }
