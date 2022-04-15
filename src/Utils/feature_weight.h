/*
 *feature_weight.h*
 the feature weight for the parser. This code is automatically generated from feature weights in other formats.

 author: Kai Zhao, Dezhong Deng
 edited by: 02/2018
*/


#ifndef FASTCKY_W
#define FASTCKY_W
extern double multi_base;
extern double multi_unpaired;
extern double multi_paired;
extern double external_unpaired;
extern double external_paired;
extern double base_pair[25];
extern double internal_1x1_nucleotides[25];
extern double helix_stacking[625];
extern double terminal_mismatch[625];
extern double bulge_0x1_nucleotides[5];
extern double helix_closing[25];
extern double dangle_left[125];
extern double dangle_right[125];
extern double internal_explicit[21];
extern double hairpin_length[31];
extern double bulge_length[31];
extern double internal_length[31];
extern double internal_symmetric_length[16];
extern double internal_asymmetry[29];
extern double hairpin_length_at_least[31];
extern double bulge_length_at_least[31];
extern double internal_length_at_least[31];
extern double internal_symmetric_length_at_least[16];
extern double internal_asymmetry_at_least[29];
#endif
