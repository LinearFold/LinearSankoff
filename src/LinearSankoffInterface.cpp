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
    float aln_weight;
    int beam_size, LFbeam, LAbeam, LAwidth;
    bool use_astar, add_branch;
    float energy_diff;
    bool is_verbose=false;

    // weight, beamsize, LFbeam, LAbeam, use_astar, energy_diff, is_verbose
    if (argc > 1) {
        aln_weight = atof(argv[1]) * 100;
        beam_size = atoi(argv[2]);
        LFbeam = atoi(argv[3]);
        LAbeam = atoi(argv[4]);
        LAwidth = atoi(argv[5]);
        use_astar = atoi(argv[6]) == 1;
        add_branch = atoi(argv[7]) == 1;
        energy_diff = atof(argv[8]);
        is_verbose = atoi(argv[9]) == 1;
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
            

            // if (beam_size == -2) {
            //     SankoffParser parser(aln_weight, beam_size, LFbeam, LAbeam, energy_diff, is_verbose);
            //     parser.parse(seqs);
            // } else {

            try{ 
                BeamSankoffParser parser(aln_weight, beam_size, LFbeam, LAbeam, LAwidth, use_astar, add_branch, energy_diff, is_verbose);
                parser.parse(seqs);
            }  
            catch (const std::overflow_error& e)
            {
                std::cout << " a overflow exception was caught, with message '"
                        << e.what() << "'\n";
            } // this executes if f() throws std::overflow_error (same type rule)
            catch (const std::exception& e) // caught by reference to base
            {
                std::cout << " a standard exception was caught, with message '"
                        << e.what() << "'\n";
            }
            // }
            // }
                
            gettimeofday(&parse_endtime, NULL);
            double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
            printf("seqs %d %d time: %f seconds.\n", seq1len, seq2len, parse_elapsed_time);
        }
    }

    return 0;
}
