#ifndef TIMOTHY_DOMDEC_H
#define TIMOTHY_DOMDEC_H
#include <mpi.h>
void decompositionInit(int argc, char **argv, MPI_Comm Comm);
void decompositionExecute();
void decompositionFinalize();

#endif /* ifndef TIMOTHY_DOMDEC_H */
