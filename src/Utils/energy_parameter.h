/*
 *energy_parameter.h*
 feature values from ViennaRNA.

 author: Dezhong Deng, He Zhang
 edited by: 02/2018
*/

#ifndef VIE_INF
// #define VIE_INF 999999999
#define VIE_INF 10000000 // to be the same as in vienna
#endif
#ifndef NBPAIRS
#define NBPAIRS 7 // NP CG GC GU UG AU UA NN
#endif

// nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U
// TODO: unify

#define SPECIAL_HP
//int special_hp = 1;

using namespace std;

extern double lxc37;
extern int ML_intern37;
extern int ML_closing37;
extern int ML_BASE37;
extern int MAX_NINIO;
extern int ninio37;
extern int TerminalAU37;  // lhuang: outermost pair is AU or GU; also used in tetra_loop triloop

// double lxc37=107.856;
// int ML_intern37=-90;
// int ML_closing37=930;
// int ML_BASE37=0;
// int MAX_NINIO=300;
// int ninio37=60;
// int TerminalAU37=50;  // lhuang: outermost pair is AU or GU; also used in tetra_loop triloop

extern char Triloops[241];
extern int Triloop37[2];

extern char Tetraloops[281];
extern int Tetraloop37[16];
extern char Hexaloops[361];
extern int Hexaloop37[4];

extern int stack37[NBPAIRS+1][NBPAIRS+1];

// lhuang                 0          1           2        3      4      5      6      7      8      9     10     
extern int hairpin37[31];
extern int bulge37[31];
extern int internal_loop37[31];

// lhuang: terminal mismatch for internal loop
extern int mismatchI37[NBPAIRS+1][5][5];
 
// lhuang: terminal mismatch for hairpin
extern int mismatchH37[NBPAIRS+1][5][5];

extern int mismatchM37[NBPAIRS+1][5][5];

extern int mismatch1nI37[NBPAIRS+1][5][5]; 
extern int mismatch23I37[NBPAIRS+1][5][5];

// lhuang: terminal mismatch for external loop
extern int mismatchExt37[NBPAIRS+1][5][5];

/* dangle5 */
extern int dangle5_37[NBPAIRS+1][5];

/* dangle3 */
extern int dangle3_37[NBPAIRS+1][5];
