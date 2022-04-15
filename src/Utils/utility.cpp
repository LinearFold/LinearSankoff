#include "utility.h"

bool _allowed_pairs[NOTON][NOTON];
bool _helix_stacking[NOTON][NOTON][NOTON][NOTON];
double cache_single[SINGLE_MAX_LEN+1][SINGLE_MAX_LEN+1];