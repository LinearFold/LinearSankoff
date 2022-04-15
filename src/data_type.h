/*
 *data_type.h*
 header file 

 author: Sizhen Li
 created by: 08/2021
*/

#ifndef DATA_TYPE_H
#define DATA_TYPE_H

using namespace std;

#ifdef dynalign
  typedef int aln_value_type;
#else
  typedef float aln_value_type;
#endif

#ifdef lv
  typedef int value_type;
#define VALUE_MIN numeric_limits<int>::lowest()
#else
  typedef double value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#endif

#define ALN_VALUE_MIN numeric_limits<aln_value_type>::lowest()
#define VALUE_FMIN numeric_limits<float>::lowest()

#endif // DATA_TYPE_H
