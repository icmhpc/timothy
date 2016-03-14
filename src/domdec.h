#ifndef TIMOTHY_DOMDEC_H
#define TIMOTHY_DOMDEC_H
#include <mpi.h>
void decompositionInit(int argc, char **argv);
void decompositionExecute(uint64_t total_number_of_cells);
void decompositionFinalize();

#endif /* ifndef TIMOTHY_DOMDEC_H */
