/*
 *HMMAlignInterface.cpp*
 The main code for HMMAlignInterface: Linear-Time simultaneous alignment of two sequences 

 author: Sizhen Li
 created by: 03/2022
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <stack>
#include <sys/time.h>
#include <algorithm>

#include "HMMAlignInterface.h"
#include "Utils/utility.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){
    vector<string> seqs;
    int beam_size;
    bool is_eval=false;
    bool is_verbose=false;

    if (argc > 1) {
        beam_size = atoi(argv[1]);
        is_eval = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
    }
    cout << "is_eval: " << is_eval << endl;

    int seq1len, seq2len;

    if (!is_eval) {
        for (string seq; getline(cin, seq);) {
            if (seq.length() == 0)
                continue;

            if (seq[0] == ';' || seq[0] == '>') {
                printf("%s\n", seq.c_str());
                continue;
            }

            // convert to uppercase
            transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

            // convert T to U
            replace(seq.begin(), seq.end(), 'T', 'U');

            printf("%s\n", seq.c_str());

            seqs.push_back(seq);

            if (seqs.size() == 1) {
                seq1len = seq.size();
            }  
            else if (seqs.size() == 2) {
                seq2len = seq.size();

                struct timeval parse_starttime, parse_endtime;
                gettimeofday(&parse_starttime, NULL);
                
                BeamAlign parser(beam_size, is_eval, is_verbose);
                parser.viterbi_path(false);

                cout << "recompute forward-backward viterbi with new parameters" << endl;
                float aln_viterbi = parser.viterbi_path(true);
                    
                gettimeofday(&parse_endtime, NULL);
                double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
                printf("seqs %d %d time: %f seconds.\n", seq1len, seq2len, parse_elapsed_time);
            }
        }
    } else { // eval mode
        int lineIndex = 0;
        vector<string> alnseqs;
        for (string seq; getline(cin, seq);) {
            if (seq.length() == 0)
                continue;

            if (seq[0] == ';' || seq[0] == '>') {
                printf("%s\n", seq.c_str());
            }

            else if (lineIndex % 3 == 1) { // seq
                // convert to uppercase
                transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
                // convert T to U
                replace(seq.begin(), seq.end(), 'T', 'U');

                // printf("%s\n", seq.c_str());

                seqs.push_back(seq);
            }

            else { // aligned seq
                // convert to uppercase
                transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
                // convert T to U
                replace(seq.begin(), seq.end(), 'T', 'U');

                printf("%s\n", seq.c_str());

                alnseqs.push_back(seq);
            } 
            lineIndex ++;
            // cout << seqs.size() << " " << alnseqs.size() << endl;
            if (seqs.size()==2 && alnseqs.size()==2) break; // stop read from file
        } // end of read from file
        
        // starttime
        struct timeval parse_starttime, parse_endtime;
        gettimeofday(&parse_starttime, NULL);

        // string to in vector
        vector<int> seq1_nucs, seq2_nucs;
        seq1len = seqs[0].size() + 1;
        seq1_nucs.resize(seq1len);
        for (int j=0; j<seq1len; j++){
            if (j == 0) seq1_nucs[j] = 4;
            else seq1_nucs[j] = GET_ACGU_NUM(seqs[0][j-1]);
        }
        seq2len = seqs[1].size() + 1;
        seq2_nucs.resize(seq2len);
        for (int j=0; j<seq2len; j++){
            if (j == 0) seq2_nucs[j] = 4;
            else seq2_nucs[j] = GET_ACGU_NUM(seqs[1][j-1]);
        }

        // constraint align pair
        assert(alnseqs[0].size() == alnseqs[1].size());
        int aln_len = alnseqs[0].size();
        cout << "aln length: " << aln_len << endl;
        int seq1_pos = 0;
        int seq2_pos = 0;
        set<pair<int, int>> align_pairs;
        for (int j=0; j<aln_len; j++) {
            char aln1_nuc = alnseqs[0][j];
            char aln2_nuc = alnseqs[1][j];
            // cout << j << " " << aln1_nuc << " " << aln2_nuc << endl;

            if (aln1_nuc == '-' && aln2_nuc == '-') break; // TODO: add error info

            if (aln1_nuc != '-' && aln2_nuc != '-') {
                seq1_pos ++;
                seq2_pos ++;
                align_pairs.insert(make_pair(seq1_pos, seq2_pos));
                // cout << "aln: " <<  aln1_nuc << " " << aln2_nuc << " " << seq1_pos << " " << seq2_pos <<  endl;
            } 
            else if (aln1_nuc != '-' && aln2_nuc == '-') {
                seq1_pos ++;
                align_pairs.insert(make_pair(seq1_pos, -1));
                // cout << "ins1: " <<  aln1_nuc << " " << aln2_nuc << " " << seq1_pos << " " << seq2_pos <<  endl;
            }
            else { // aln1_nuc == "-" && aln2_nuc != "-"
                seq2_pos ++;
                align_pairs.insert(make_pair(-1, seq2_pos));
                // cout << "ins2: " <<  aln1_nuc << " " << aln2_nuc << " " << seq1_pos << " " << seq2_pos <<  endl;
            }
        }


        BeamAlign parser(beam_size, is_eval, is_verbose);
        parser.set(beam_size, seq1_nucs, seq2_nucs);
        parser.viterbi_path(false);

        // realign use new parameters
        parser.viterbi_path(true);

        cout << "eval..." << endl;
        float aln_viterbi = parser.evalulate(true, align_pairs);
            
        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        printf("seqs %d %d time: %f seconds.\n", seq1len, seq2len, parse_elapsed_time);
    }
    
    return 0;
}