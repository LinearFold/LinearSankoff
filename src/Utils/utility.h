/*
 *utility.h*
 provides feature functions.

 author: Kai Zhao, Dezhong Deng
 edited by: 02/2018
*/

#ifndef FASTCKY_UTILITY_H
#define FASTCKY_UTILITY_H

#include <algorithm>
#include <cstring>
#include <assert.h>
#include <math.h> 
#include <vector>

#include "feature_weight.h"

#define INF 1000000007

#define NOTON 5 // NUM_OF_TYPE_OF_NUCS
#define NOTOND 25
#define NOTONT 125

#define EXPLICIT_MAX_LEN 4
#define SINGLE_MIN_LEN 0
#define SINGLE_MAX_LEN 30  // NOTE: *must* <= sizeof(char), otherwise modify State::TraceInfo accordingly

#define MAX_LOOP_LEN 20  // NOTE: *must* <= sizeof(char), otherwise modify State::TraceInfo accordingly

#define HAIRPIN_MAX_LEN 30
#define BULGE_MAX_LEN SINGLE_MAX_LEN
#define INTERNAL_MAX_LEN SINGLE_MAX_LEN
#define SYMMETRIC_MAX_LEN 15
#define ASYMMETRY_MAX_LEN 28

#define GET_ACGU_NUM(x) ((x=='A'? 0 : (x=='C'? 1 : (x=='G'? 2 : (x=='U'?3: 4)))))
#define GET_NUC(x) ((x==0? 'A' : (x==1 ? 'C' : (x==2 ? 'G' : (x==3 ? 'U': '.')))))
#define HELIX_STACKING_OLD(x, y, z, w) (_helix_stacking[GET_ACGU_NUM(x)][GET_ACGU_NUM(y)][GET_ACGU_NUM(z)][GET_ACGU_NUM(w)])

extern bool _allowed_pairs[NOTON][NOTON];
extern bool _helix_stacking[NOTON][NOTON][NOTON][NOTON];
extern double cache_single[SINGLE_MAX_LEN+1][SINGLE_MAX_LEN+1];

inline void initialize_cachesingle()
{
    memset(cache_single, 0, sizeof(cache_single));
    for (int l1 = SINGLE_MIN_LEN; l1 <= SINGLE_MAX_LEN; l1 ++)
        for (int l2 = SINGLE_MIN_LEN; l2 <= SINGLE_MAX_LEN; l2 ++)
        {
            if (l1 == 0 && l2 == 0)
                continue;

                // bulge
            else if (l1 == 0)
                cache_single[l1][l2] += bulge_length[l2];
            else if (l2 == 0)
                cache_single[l1][l2] += bulge_length[l1];
            else
            {

                // internal
                cache_single[l1][l2] += internal_length[std::min(l1+l2, INTERNAL_MAX_LEN)];

                // internal explicit
                if (l1 <= EXPLICIT_MAX_LEN && l2 <= EXPLICIT_MAX_LEN)
                    cache_single[l1][l2] +=
                            internal_explicit[l1<=l2 ? l1*EXPLICIT_MAX_LEN+l2 : l2*EXPLICIT_MAX_LEN+l1];

                // internal symmetry
                if (l1 == l2)
                    cache_single[l1][l2] += internal_symmetric_length[std::min(l1, SYMMETRIC_MAX_LEN)];

                else {  // internal asymmetry
                    int diff = l1 - l2; if (diff < 0) diff = -diff;
                    cache_single[l1][l2] += internal_asymmetry[std::min(diff, ASYMMETRY_MAX_LEN)];
                }
            }
        }
    return;
}

inline void initialize()
{
    _allowed_pairs[GET_ACGU_NUM('A')][GET_ACGU_NUM('U')] = true;
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('A')] = true;
    _allowed_pairs[GET_ACGU_NUM('C')][GET_ACGU_NUM('G')] = true;
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('C')] = true;
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('U')] = true;
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('G')] = true;

    HELIX_STACKING_OLD('A', 'U', 'A', 'U') = true;
    HELIX_STACKING_OLD('A', 'U', 'C', 'G') = true;
    HELIX_STACKING_OLD('A', 'U', 'G', 'C') = true;
    HELIX_STACKING_OLD('A', 'U', 'G', 'U') = true;
    HELIX_STACKING_OLD('A', 'U', 'U', 'A') = true;
    HELIX_STACKING_OLD('A', 'U', 'U', 'G') = true;
    HELIX_STACKING_OLD('C', 'G', 'A', 'U') = true;
    HELIX_STACKING_OLD('C', 'G', 'C', 'G') = true;
    HELIX_STACKING_OLD('C', 'G', 'G', 'C') = true;
    HELIX_STACKING_OLD('C', 'G', 'G', 'U') = true;
    HELIX_STACKING_OLD('C', 'G', 'U', 'G') = true;
    HELIX_STACKING_OLD('G', 'C', 'A', 'U') = true;
    HELIX_STACKING_OLD('G', 'C', 'C', 'G') = true;
    HELIX_STACKING_OLD('G', 'C', 'G', 'U') = true;
    HELIX_STACKING_OLD('G', 'C', 'U', 'G') = true;
    HELIX_STACKING_OLD('G', 'U', 'A', 'U') = true;
    HELIX_STACKING_OLD('G', 'U', 'G', 'U') = true;
    HELIX_STACKING_OLD('G', 'U', 'U', 'G') = true;
    HELIX_STACKING_OLD('U', 'A', 'A', 'U') = true;
    HELIX_STACKING_OLD('U', 'A', 'G', 'U') = true;
    HELIX_STACKING_OLD('U', 'G', 'G', 'U') = true;
}

// ------------- nucs based scores -------------

// parameters: nucs[i], nucs[j]
inline double base_pair_score(int nuci, int nucj) {
  return base_pair[nucj*NOTON + nuci];
}

// parameters: nucs[i], nucs[i+1], nucs[j-1], nucs[j]
inline double helix_stacking_score(int nuci, int nuci1, int nucj_1, int nucj) {
  return helix_stacking[nuci*NOTONT + nucj*NOTOND + nuci1*NOTON + nucj_1];
}

// parameters: nucs[i], nucs[j]
inline double helix_closing_score(int nuci, int nucj) {
    return helix_closing[nuci*NOTON + nucj];
}

// parameters: nucs[i], nucs[i+1], nucs[j-1], nucs[j]
inline double terminal_mismatch_score(int nuci, int nuci1, int nucj_1, int nucj) {
    return terminal_mismatch[nuci*NOTONT+nucj*NOTOND + nuci1*NOTON + nucj_1];
}



// parameter: nucs[i]
inline double bulge_nuc_score(int nuci) {
    return bulge_0x1_nucleotides[nuci];
}

// parameters: nucs[i], nucs[j]
inline double internal_nuc_score(int nuci, int nucj) {
  return internal_1x1_nucleotides[nuci*NOTON + nucj];
}

// parameters: nucs[i], nucs[i+1], nucs[j]
inline double dangle_left_score(int nuci, int nuci1, int nucj) {
    return dangle_left[nuci*NOTOND + nucj*NOTON + nuci1];
}

// parameters: nucs[i], nucs[j-1], nucs[j]
inline double dangle_right_score(int nuci, int nucj_1, int nucj) {
    return dangle_right[nuci*NOTOND + nucj*NOTON + nucj_1];
}



// ------------- length based scores -------------

inline double hairpin_score(int i, int j) {
    return hairpin_length[std::min(j-i-1, HAIRPIN_MAX_LEN)];
}

inline double internal_length_score(int l) {
    return internal_length[std::min(l, INTERNAL_MAX_LEN)];
}

inline double internal_explicit_score(int l1, int l2){
    int l1_ = std::min(l1, EXPLICIT_MAX_LEN);
    int l2_ = std::min(l2, EXPLICIT_MAX_LEN);
    return internal_explicit[l1_<=l2_ ? l1_*NOTON+l2_ : l2_*NOTON+l1_];
}

inline double internal_sym_score(int l) {
    return internal_symmetric_length[std::min(l, SYMMETRIC_MAX_LEN)];
}

inline double internal_asym_score(int l1, int l2)
{
    int diff = l1 - l2; if (diff < 0) diff = -diff;
    return internal_asymmetry[std::min(diff, ASYMMETRY_MAX_LEN)];
}

inline double bulge_length_score(int l){
    return bulge_length[std::min(l, BULGE_MAX_LEN)];
}

inline double hairpin_at_least_score(int l) {
    return hairpin_length_at_least[std::min(l, HAIRPIN_MAX_LEN)];
}

inline double buldge_length_at_least_score(int l) {
    return bulge_length_at_least[std::min(l, BULGE_MAX_LEN)];
}

inline double internal_length_at_least_score(int l) {
    return internal_length_at_least[std::min(l, INTERNAL_MAX_LEN)];
}


//-----------------------------------------------------
inline double score_junction_A(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    return helix_closing_score(nuci, nucj) +
            (i < len - 1 ? dangle_left_score(nuci, nuci1, nucj) : 0) +
            (j > 0 ? dangle_right_score(nuci, nucj_1, nucj) : 0);
}

inline double score_junction_B(int i, int j, int nuci, int nuci1, int nucj_1, int nucj) {
    return helix_closing_score(nuci, nucj) + terminal_mismatch_score(nuci, nuci1, nucj_1, nucj);
}

inline double score_hairpin_length(int len) {
  return hairpin_length[std::min(len, HAIRPIN_MAX_LEN)];
}

inline double score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj) {
    return hairpin_length[std::min(j-i-1, HAIRPIN_MAX_LEN)] +
            score_junction_B(i, j, nuci, nuci1, nucj_1, nucj);
}

inline double score_helix(int nuci, int nuci1, int nucj_1, int nucj) {
    return helix_stacking_score(nuci, nuci1, nucj_1, nucj) + base_pair_score(nuci1, nucj_1);
}

inline double score_single_nuc(int i, int j, int p, int q, int nucp_1, int nucq1) {
    int l1 = p-i-1, l2=j-q-1;
    if (l1==0 && l2==1) return bulge_nuc_score(nucq1);
    if (l1==1 && l2==0) return bulge_nuc_score(nucp_1);
    if (l1==1 && l2==1) return internal_nuc_score(nucp_1, nucq1);
    return 0;
}

inline double score_single(int i, int j, int p, int q, int len,
                           int nuci, int nuci1, int nucj_1, int nucj,
                           int nucp_1, int nucp, int nucq, int nucq1) {
    int l1 = p-i-1, l2=j-q-1;
    return cache_single[l1][l2] +
           base_pair_score(nucp, nucq) +
           score_junction_B(i, j, nuci, nuci1, nucj_1, nucj) +
           score_junction_B(q, p, nucq, nucq1, nucp_1, nucp) +
           score_single_nuc(i, j, p, q, nucp_1, nucq1);
}

// score_single without socre_junction_B
inline double score_single_without_junctionB(int i, int j, int p, int q,
                           int nucp_1, int nucp, int nucq, int nucq1) {
    int l1 = p-i-1, l2=j-q-1;
    return cache_single[l1][l2] +
           base_pair_score(nucp, nucq) +
           score_single_nuc(i, j, p, q, nucp_1, nucq1);
}

inline double score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    return score_junction_A(i, j, nuci, nuci1, nucj_1, nucj, len) +
           multi_paired + multi_base;
}

inline double score_multi_unpaired(int i, int j) {
    return (j-i+1) * multi_unpaired;
}

inline double score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len) {
    return score_junction_A(k, i, nuck, nuck1, nuci_1, nuci, len) +
           score_multi_unpaired(k+1, j) + base_pair_score(nuci, nuck) + multi_paired;
}

inline double score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len) {
    return score_junction_A(j, i, nucj, nucj1, nuci_1, nuci, len) +
           external_paired + base_pair_score(nuci, nucj);
}

inline double score_external_unpaired(int i, int j) {
    return (j-i+1) * external_unpaired;
}

// data type and 
// quick select 
// log space (borrowed from CONTRAfold)

#define NEG_INF -2e20 

#ifdef lv
  typedef float pf_type;
#else
  typedef double pf_type;
#endif

#ifdef lv
  typedef int value_type;
  #define VALUE_MIN_PF numeric_limits<double>::lowest()
  #define VALUE_MIN numeric_limits<int>::lowest()
#else
  typedef double value_type;
  #define VALUE_MIN_PF numeric_limits<double>::lowest()
  #define VALUE_MIN numeric_limits<double>::lowest()
#endif

#ifdef dynalign
  typedef int aln_value_type;
#else
  typedef float aln_value_type;
#endif
#define ALN_VALUE_MIN numeric_limits<aln_value_type>::lowest()
#define VALUE_FMIN numeric_limits<float>::lowest()

// quick select (data type: int)
inline unsigned long quickselect_partition(std::vector<std::pair<int, int>>& scores, unsigned long lower, unsigned long upper) {
    int pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}
// in-place quick-select
inline int quickselect(std::vector<std::pair<int, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

// quick select (data type: float)
inline unsigned long quickselect_partition(std::vector<std::pair<float, int>>& scores, unsigned long lower, unsigned long upper) {
    float pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}
// in-place quick-select
inline float quickselect(std::vector<std::pair<float, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


inline pf_type Fast_LogExpPlusOne(pf_type x){
    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.
    
    assert(pf_type(0.0000000000) <= x && x <= pf_type(11.8624794162) && "Argument out-of-range.");
    if (x < pf_type(3.3792499610))
    {
        if (x < pf_type(1.6320158198))
        {
            if (x < pf_type(0.6615367791))
                return ((pf_type(-0.0065591595)*x+pf_type(0.1276442762))*x+pf_type(0.4996554598))*x+pf_type(0.6931542306);
            return ((pf_type(-0.0155157557)*x+pf_type(0.1446775699))*x+pf_type(0.4882939746))*x+pf_type(0.6958092989);
        }
        if (x < pf_type(2.4912588184))
            return ((pf_type(-0.0128909247)*x+pf_type(0.1301028251))*x+pf_type(0.5150398748))*x+pf_type(0.6795585882);
        return ((pf_type(-0.0072142647)*x+pf_type(0.0877540853))*x+pf_type(0.6208708362))*x+pf_type(0.5909675829);
    }
    if (x < pf_type(5.7890710412))
    {
        if (x < pf_type(4.4261691294))
            return ((pf_type(-0.0031455354)*x+pf_type(0.0467229449))*x+pf_type(0.7592532310))*x+pf_type(0.4348794399);
        return ((pf_type(-0.0010110698)*x+pf_type(0.0185943421))*x+pf_type(0.8831730747))*x+pf_type(0.2523695427);
    }
    if (x < pf_type(7.8162726752))
        return ((pf_type(-0.0001962780)*x+pf_type(0.0046084408))*x+pf_type(0.9634431978))*x+pf_type(0.0983148903);
    return ((pf_type(-0.0000113994)*x+pf_type(0.0003734731))*x+pf_type(0.9959107193))*x+pf_type(0.0149855051);
}

inline void Fast_LogPlusEquals (pf_type &x, pf_type y)
{
    if (x < y) std::swap (x, y);
    if (y > pf_type(NEG_INF/2) && x-y < pf_type(11.8624794162))
        x = Fast_LogExpPlusOne(x-y) + y;
}

inline pf_type Fast_Exp(pf_type x)
{
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.
    
    if (x < pf_type(-2.4915033807))
    {
        if (x < pf_type(-5.8622823336))
        {
            if (x < pf_type(-9.91152))
                return pf_type(0);
            return ((pf_type(0.0000803850)*x+pf_type(0.0021627428))*x+pf_type(0.0194708555))*x+pf_type(0.0588080014);
        }
        if (x < pf_type(-3.8396630909))
            return ((pf_type(0.0013889414)*x+pf_type(0.0244676474))*x+pf_type(0.1471290604))*x+pf_type(0.3042757740);
        return ((pf_type(0.0072335607)*x+pf_type(0.0906002677))*x+pf_type(0.3983111356))*x+pf_type(0.6245959221);
    }
    if (x < pf_type(-0.6725053211))
    {
        if (x < pf_type(-1.4805375919))
            return ((pf_type(0.0232410351)*x+pf_type(0.2085645908))*x+pf_type(0.6906367911))*x+pf_type(0.8682322329);
        return ((pf_type(0.0573782771)*x+pf_type(0.3580258429))*x+pf_type(0.9121133217))*x+pf_type(0.9793091728);
    }
    if (x < pf_type(0))
        return ((pf_type(0.1199175927)*x+pf_type(0.4815668234))*x+pf_type(0.9975991939))*x+pf_type(0.9999505077);
    return (x > pf_type(46.052) ? pf_type(1e20) : expf(x));
}


#endif //FASTCKY_UTILITY_H
