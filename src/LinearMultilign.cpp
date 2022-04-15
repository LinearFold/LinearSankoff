/*
 *LinearMultilign.cpp*
 The main code for LinearMultilign: Linear-Time simultaneous folding and alignment of multiple sequences 

 author: Sizhen Li
 created by: 09/2021
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <stack>
#include <sys/time.h>

#include "LinearMultilign.h"

// BeamMultilignParser::BeamMultilignParser(int beam_size, int beam_size2)
//     :beam(beam_size),
//      alnbeam(beam_size2){}

int main(int argc, char** argv){
    int beam_size, beam_size2;
    float aln_weight;
    bool is_verbose;

    if (argc > 1) {
        beam_size = atoi(argv[1]);
        beam_size2 = atoi(argv[2]);
        aln_weight = atof(argv[3]);
        is_verbose = atoi(argv[4]) == 1;
    }

    // BeamMultilignParser parser(beam_size, beam_size2);
    cout << "beam size: " << beam_size << " " << beam_size2 << endl;

    // load sequences
    vector<string> seqs; 
    int avg_seq_len = 0;
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
        avg_seq_len += seq.size();
    }
    int num_seqs = seqs.size();
    avg_seq_len /= num_seqs;
    printf("number of sequences: %d\n", num_seqs);
    cout << "average sequence length: " << avg_seq_len << endl;

    // progressively call dynalign
    BeamSankoffParser* parser = new BeamSankoffParser(beam_size, beam_size2, aln_weight, is_verbose);
    vector<string> twoseqs; 
    int i_seq1 = 0;
    int seq1len = seqs[i_seq1].size();

    // limited allowed base pairs
    set<pair<int, int>> allowed_pairs;  
    vector<pair<int, int>> out_pairs; 
    for (int i_iter=0; i_iter<2; i_iter++) {
        for (int i_seq2=1; i_seq2<num_seqs; i_seq2++) {
            // collect two seqs
            cout << "========== " << i_iter << " ========== " << i_seq1 << " ========== " << i_seq2 << " ==========" << endl;
            // cout << i_iter << " " << i_seq1 << " " << i_seq2 << endl;
            int seq2len = seqs[i_seq2].size();
            // cout << "seq length: " << seq1len << " " << seq2len << endl;
            twoseqs.push_back(seqs[i_seq1]);
            twoseqs.push_back(seqs[i_seq2]); 

            // call dynalign
            if (i_iter==0 && i_seq2 == 1) // no constraint
                parser->parse(twoseqs, false, allowed_pairs, out_pairs, avg_seq_len); 
            else // constraint folding
                parser->parse(twoseqs, true, allowed_pairs, out_pairs, avg_seq_len); 

            // reset allowed_pairs
            allowed_pairs.clear();
            for (const auto& item : out_pairs) {
                allowed_pairs.insert(item);
                // cout << "allowed_pairs: " << get<0>(item) << " " << get<1>(item) << endl;
            }
            out_pairs.clear();
            twoseqs.clear();
        }
    }
    cout << "Done!" << endl;
    delete parser;
}
