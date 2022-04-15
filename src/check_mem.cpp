#include "check_mem.h"

using namespace std;

// memory check
int memoryparseLine(char *line) {
  int i = strlen(line);
  const char *p = line;
  while (*p < '0' || *p > '9') p++;
  line[i - 3] = '\0';
  i = atoi(p);
  return i;
}

processMem_t GetProcessMemory() {
  FILE *file = fopen("/proc/self/status", "r");
  char line[128];
  processMem_t processMem;
  int flag = 0;

  while (fgets(line, 128, file) != NULL) {
    if (strncmp(line, "VmSize:", 7) == 0) {
      processMem.virtualMem = memoryparseLine(line);
      flag++;
    }

    if (strncmp(line, "VmPeak:", 7) == 0) {
      processMem.VmPeak = memoryparseLine(line);
      flag++;
    }

    if (flag == 2) break;
  }
  fclose(file);
  return processMem;
}