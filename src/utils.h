#ifndef TIMOTHY_UTILS_H
#define TIMOTHY_UTILS_H
#include<stdint.h>
#include<string.h>
void stopRun(int, const char*, const char*, int);
void swap_Nbyte(char *data, uint64_t n, int m);
void checkEndiannes();
void getLocalRankAndSize(int rank, int size, int32_t * lrank,	int32_t *lsize);
size_t getMemoryPerProcess(int32_t lsize);
#endif //TIMOTHY_UTILS_H

