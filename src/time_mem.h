/*
 *time_mem.h*
 header file for time_mem.cpp.

 author: Sizhen Li
 created by: 08/2021
*/

#ifndef TIME_MEMORY_H
#define TIME_MEMORY_H

// memory check
typedef struct {
   unsigned virtualMem;
   unsigned VmPeak;
} processMem_t;

int memoryparseLine(char *line);

processMem_t GetProcessMemory();

#endif
