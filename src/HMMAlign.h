/*
 *HMMAlign.h*
 header file for HMMAlign.cpp.

 author: Sizhen Li
 created by: 08/2021
*/

#ifndef HMM_ALIGN_H
#define HMM_ALIGN_H

#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <math.h> 

#include "Utils/utility.h"
// #include "LinearSankoff.h"

using namespace std;

#define GET_NUC_HMM(x) ((x==0? 'A' : (x==1 ? 'C' : (x==2 ? 'G' : (x==3 ? 'U': 'N')))))


// from LinearPartition
// inline float Fast_Exp(float x)
// {
//     // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
//     // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
//     // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
//     // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
//     // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
//     // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
//     // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
//     // 6 polynomials needed.
    
//     if (x < float(-2.4915033807))
//     {
//         if (x < float(-5.8622823336))
//         {
//             if (x < float(-9.91152))
//                 return float(0);
//             return ((float(0.0000803850)*x+float(0.0021627428))*x+float(0.0194708555))*x+float(0.0588080014);
//         }
//         if (x < float(-3.8396630909))
//             return ((float(0.0013889414)*x+float(0.0244676474))*x+float(0.1471290604))*x+float(0.3042757740);
//         return ((float(0.0072335607)*x+float(0.0906002677))*x+float(0.3983111356))*x+float(0.6245959221);
//     }
//     if (x < float(-0.6725053211))
//     {
//         if (x < float(-1.4805375919))
//             return ((float(0.0232410351)*x+float(0.2085645908))*x+float(0.6906367911))*x+float(0.8682322329);
//         return ((float(0.0573782771)*x+float(0.3580258429))*x+float(0.9121133217))*x+float(0.9793091728);
//     }
//     if (x < float(0))
//         return ((float(0.1199175927)*x+float(0.4815668234))*x+float(0.9975991939))*x+float(0.9999505077);
//     return (x > float(46.052) ? float(1e20) : expf(x));
// }

// ----------------------------------------
// log computation from RNAstructure
#define USE_XLOG_ZERO
// #define DBL_MAX __DBL_MAX__
#ifdef USE_XLOG_ZERO
    static const float LOG_OF_ZERO = VALUE_FMIN;
    #define IS_LOG_ZERO(x) (x<=LOG_OF_ZERO)
#endif
// ----------------------------------------


// ----------------------------------------
// HMM paramethers from RNAstructure
#define N_BINZ (10)
#define N_STATES (3)
#define N_OUTPUTS (27)

static float ML_emit_probs[N_OUTPUTS][N_STATES] = 
{
    {0.000000, 0.000000, 0.134009}, // AA 0
    {0.000000, 0.000000, 0.027164}, // AC
    {0.000000, 0.000000, 0.049659}, // AG
    {0.000000, 0.000000, 0.028825}, // AU
    {0.211509, 0.000000, 0.000000}, // A.
    {0.000000, 0.000000, 0.027164}, // CA 5
    {0.000000, 0.000000, 0.140242}, // CC
    {0.000000, 0.000000, 0.037862}, // CG
    {0.000000, 0.000000, 0.047735}, // CU
    {0.257349, 0.000000, 0.000000}, // C. 
    {0.000000, 0.000000, 0.049659}, // GA 10
    {0.000000, 0.000000, 0.037862}, // GC
    {0.000000, 0.000000, 0.178863}, // GG
    {0.000000, 0.000000, 0.032351}, // GU
    {0.271398, 0.000000, 0.000000}, // G.
    {0.000000, 0.000000, 0.028825}, // UA 15 
    {0.000000, 0.000000, 0.047735}, // UC
    {0.000000, 0.000000, 0.032351}, // UG
    {0.000000, 0.000000, 0.099694}, // UU
    {0.259744, 0.000000, 0.000000}, // U.
    {0.000000, 0.211509, 0.000000}, // .A 20 
    {0.000000, 0.257349, 0.000000}, // .C 
    {0.000000, 0.271398, 0.000000}, // .G
    {0.000000, 0.259744, 0.000000}, // .U
    {0.000000, 0.000000, 0.000000}, // ..
    {0.000000, 0.000000, 1.000000}, // START 25
    {0.000000, 0.000000, 1.000000}  // END
};

static float ML_trans_probs[N_STATES][N_STATES] = 
{
    {0.666439, 0.041319, 0.292242}, // INS1
    {0.041319, 0.666439, 0.292242}, // INS2
    {0.022666, 0.022666, 0.954668}  // ALIGN
};
// ----------------------------------------

// ----------------------------------------
// log computation from RNAstructure
inline float xlog(const float& value) noexcept {
    #ifdef USE_XLOG_ZERO
    if (value == 0.0) return LOG_OF_ZERO; // the log function itself will throw an exception if value < 0.
    #endif
    return log(value); // the log function itself will throw an exception if value < 0.
}
inline float xlog_mul(const float& log1, const float& log2) {
    #ifdef USE_XLOG_ZERO
    return (log1 <= LOG_OF_ZERO || log2 <= LOG_OF_ZERO) ? LOG_OF_ZERO : log1+log2;
    #else
    return log1+log2;
    #endif
}
// Computes log(exp(a)+exp(b))
inline float xlog_sum(const float& a, const float& b) {
	// Derivation:   log(exp(a)-exp(b))  
	//             = log(exp(a) * (1-exp(b-a)) )
	//             = a+log(1-exp(b-a))
	//             = a+log1p(-exp(b-a))
    #ifdef USE_XLOG_ZERO
	if(IS_LOG_ZERO(a)) return b;
    if(IS_LOG_ZERO(b)) return a;
	#endif
	// Note: The test of a>b is important when A or B is greater than MAX_DOUBLE.
	// e.g. A=1E+500 and B=1 --- exp(a-b) will overflow, while exp(b-a) will not.
	// As long as a>=b, the value of exp(b-a) will always be in the range (0,1], 
	// so it will not overflow (even when one or both numbers are greater than 
	// MAX_DOUBLE in the linear scale).
	// In fact, the result will always be max(a,b)+Q where Q is in the range [0..ln(2)].
	return a>b ? a+log1p(exp(b-a)): b+log1p(exp(a-b));
	/* 
	// -----Possible alternatives (todo: test for speed)-----
    // -----alt.1----- (temporary vars for min,max)
	const double mx = std::max(a,b), mn=std::min(a,b);
    return mx+log1pexp(mn-mx);
	// -----alt.2----- (implicit compiler-defined temporary vars for min,max)
    return std::max(a,b)+log1pexp(std::min(a,b)-std::max(a,b));
	// -----alt.3----- (max and -abs)
	return std::max(a,b)+log1pexp(-std::abs(a-b));
	*/
}
// Returns 0 if log1 is 0 no matter what log2 is.
inline double xlog_div(const double& log1, const double& log2) {
	#ifdef USE_XLOG_ZERO
	if(log1 <= LOG_OF_ZERO) return LOG_OF_ZERO;
	if(log2 <= LOG_OF_ZERO) throw std::runtime_error("Division by xlog zero-value (in " __FILE__ ")" );
	#endif
	return log1-log2;
}
// ----------------------------------------

enum HMMManner {
  HMMMANNER_NONE = 0,        // 0: empty
  MANNER_INS1 = 1,           // 1: INS1
  MANNER_INS2 = 2,           // 2: INS2
  MANNER_ALN = 3             // 3: ALN
};

struct AlnScore {
    float prob;
    HMMManner manner;

    AlnScore(): prob(xlog(0)), manner(HMMMANNER_NONE) {};

     void set(float score_, HMMManner manner_) {
        prob = score_; 
        manner = manner_;
    }
};

struct AlignState {
    int i;
    int k;
    float ml;
    float beta;
    HMMManner start_manner;
    // HMMManner manner;
    HMMManner pre;
    int step;

    // AlignState(): ml(xlog(0)), manner(HMMMANNER_NONE), start_manner(HMMMANNER_NONE), pre(HMMMANNER_NONE), step(0), i(0), k(0) {};
    AlignState(): ml(xlog(0)), beta(xlog(0)), pre(HMMMANNER_NONE), start_manner(HMMMANNER_NONE), step(0), i(0), k(0) {};

    void set(float score_, HMMManner pre_manner_) {
        ml = score_; 
        pre = pre_manner_;
    }

    void set(float score_, HMMManner pre_manner_, int step_, int i_, int k_) {
        ml = score_; 
        pre = pre_manner_;
        step = step_;
        i = i_;
        k = k_;
    }

    void set(float score_, HMMManner pre_manner_, int step_, int i_, int k_, HMMManner start_manner_) {
        ml = score_; 
        pre = pre_manner_;
        step = step_;
        i = i_;
        k = k_;
        start_manner = start_manner_;
    }

    void set_backward(float score_) {
        beta = score_; 
    }
};

struct AlignState3 {
    int i;
    int k;
    // float ml;
    // int step;

    AlignState alnobj;
    AlignState ins1obj;
    AlignState ins2obj;

    AlignState3(): alnobj(), ins1obj(), ins2obj(), i(0), k(0) {};
};

class BeamAlign{
public:
    int beam; // , max_len;
    int start1, start2; // end1, end2;
    vector<int> seq1, seq2;
    int seq1_len, seq2_len;
    int max_range;
    bool eval=false;
    bool verbose=false;

    BeamAlign(int beam_size=100,
              bool is_eval=false,
              bool is_verbose=false);

    vector<int> up_bounds;
    vector<int> low_bounds;
    vector<int> max_j1;
    vector<int> min_j1;

    float** emission_probs;
	float** trans_probs;

    unordered_map<int, AlignState> *bestINS1, *bestINS2, *bestALN;

    // local alignment scores
    vector<unordered_map<int, AlignState3> > local_scores;
    float *****all_local_scores;
    float ******left_local_scores;
    float ******right_local_scores;

    void set(int beam_size, vector<int> &seq1_nuc_types, vector<int> &seq2_nuc_types);
    void clear(bool local_path=true);

    float viterbi_path(bool newpara);
    float evalulate(bool newpara, const std::set<pair<int, int>> &align_pairs);
    void viterbi_path_all_locals();
    float viterbi_path_local(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose=false);
    float viterbi_path_local_left(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose=false);
    float viterbi_path_local_right(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool limited=true, bool verbose=false);
    void traceback2(int i1, int j1, int i2, int j2, vector<char> &aln1, vector<char> &aln2, HMMManner endmanner);

    float get_trans_emit_prob0(int prev_state, int current_state, int i, int k, bool new_pars=false);
    float get_trans_emit_prob1(int current_state, int next_state, int i, int k, bool new_pars=false);
private:
    float* fam_hmm_pars;
	float* fam_thresholds;
    float forward_score;

    vector<unordered_map<int, AlnScore> > aln_env;

    // pre-computed
    void viterbi_path_local22(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose=false);
    void viterbi_path_local_left22(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose=false);
    void viterbi_path_local_right22(int i1, int j1, int i2, int j2, HMMManner s1, HMMManner s2, bool verbose=false);

    void prepare(int j1, int j2);
    bool update_if_better(AlignState &state, float newscore, HMMManner pre_manner, int step, int i, int k, HMMManner start_manner=HMMMANNER_NONE); 
    bool update_if_better_backward(AlignState &state, float newscore); 
    void update(AlignState &state, float newscore, HMMManner pre_manner, int step, int i, int k); 
    
    // traceback
    float get_aln_similarity(char* &seq1_aln_line, char* &seq2_aln_line, char gap_symbol='-');
    void set_parameters_by_sim(float similarity);
    int get_bin_index(float similarity, int n_bins);
    void load_init_params();

    // beam prune
    void beam_prune(std::unordered_map<int, AlignState> &beamstep);
    vector<pair<float, int> > scores;

    float get_trans_emit_prob(int prev_state, int current_state, int i_1, int k_1, int i, int k, HMMManner s1=HMMMANNER_NONE, HMMManner s2=HMMMANNER_NONE, bool new_pars=false);
    float get_trans_emit_prob(int prev_state, int current_state, int i, int k, bool new_pars=false);
    float get_trans_emit_prob_right(int prev_state, int current_state, int i_1, int k_1, int i, int k, HMMManner s1=HMMMANNER_NONE, HMMManner s2=HMMMANNER_NONE, bool new_pars=false);
    float get_trans_emit_prob_left(int current_state, int next_state, int i, int k, HMMManner s1=HMMMANNER_NONE, HMMManner s2=HMMMANNER_NONE, bool new_pars=false);
    
    void traceback(vector<char> &aln1, vector<char> &aln2, HMMManner endmanner);

    void viterbi_backward(bool newpara);
    void forward();
	void backward();
    void cal_align_prob(float threshold);
};

#endif // HMM_ALIGN_H