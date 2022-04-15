/*
 *LinearSankoff_Interface.cpp*
 The main code for LinearSankoff_Interface: Linear-Time simultaneous folding and alignment of two sequences 

 author: Sizhen Li
 created by: 03/2022
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <sys/time.h>

#include "LinearSankoffInterface.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
    vector<string> seqs;
    int beam_size, beam_size2;
    float aln_weight;
    float energy_diff;
    bool is_verbose=false;

    if (argc > 1) {
        beam_size = atoi(argv[1]);
        beam_size2 = atoi(argv[2]);
        aln_weight = atof(argv[3]);
        energy_diff = atof(argv[4]);
        is_verbose = atoi(argv[5]) == 1;
    }

    int seq1len, seq2len;
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
            
            // if (seq1len < seq2len){ // TODO
            //     BeamSankoffParser parser(beam_size, beam_size2, aln_weight, is_verbose);
            //     swap(seqs[0], seqs[1]);
            //     parser.parse(seqs);
            // }
            // else {
            if (beam_size == -2) {
                SankoffParser parser(beam_size, beam_size2, aln_weight, energy_diff, is_verbose);
                parser.parse(seqs);
            } else {
                BeamSankoffParser parser(beam_size, beam_size2, aln_weight, energy_diff, is_verbose);
                parser.parse(seqs);
            }
            // }
                
            gettimeofday(&parse_endtime, NULL);
            double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
            printf("seqs %d %d time: %f seconds.\n", seq1len, seq2len, parse_elapsed_time);
        }
    }

    return 0;
}
