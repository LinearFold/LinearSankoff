/*
 *time_mem.h*
 header file for time_mem.cpp.

 author: Sizhen Li
 created by: 08/2021
*/

#ifndef CHECK_MEMORY_H
#define CHECK_MEMORY_H

#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

// memory check
typedef struct {
   unsigned virtualMem;
   unsigned VmPeak;
} processMem_t;

int memoryparseLine(char *line);

processMem_t GetProcessMemory();

#endif
