#ifndef TIMOTHY_DOMDEC_H
#define TIMOTHY_DOMDEC_H
#include <mpi.h>
void decompositionInit(int argc, char **argv);
void decompositionExecute();
void decompositionFinalize();

#endif /* ifndef TIMOTHY_DOMDEC_H */
