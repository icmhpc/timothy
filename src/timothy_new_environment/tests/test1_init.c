//
// Created by Grzegorz Bokota on 12.06.2016.
//
#include <stdlib.h>
#include "../environment.h"
struct gridData grid;
struct settings simsetup;
struct state simstate;

void prepareTestEnvironment(struct settings *set) {
  size_t i;
  set->environments =
      calloc(set->numberOfEnvironments, sizeof(struct environment));
  for (i = 0; i < set->numberOfEnvironments; i++) {
    set->environments[i].diffusion_coefficient = 1.82e-5;
    set->environments[i].lambda_delay = 0.0;
    set->environments[i].boundary_condition = ((double)i + 1.0) * 0.1575e-6;
    set->environments[i].initial_condition_mean = ((double)i + 1.0) * 0.1575e-6;
    set->environments[i].initial_condition_variance = 0.0;
  }
}

extern void computeGridSize(const struct settings *simsetup,
                     const struct state *simstate, struct gridData *grid);

extern void allocateGrid(const struct settings *simsetup, const struct state *simstate,
                  struct gridData *grid);

int main(int argc, char * argv[]){
  simstate.MPIrank = 0;
  simstate.MPIsize = 27;
  simstate.MPIdim[0] = 3;
  simstate.MPIdim[1] = 3;
  simstate.MPIdim[2] = 3;
  simsetup.numberOfEnvironments = 5;
  simsetup.gfH = 128.0;
  simsetup.maxCellsPerProc = 300000;
  simsetup.size_x = 1024;
  simsetup.size_y = 1024;
  simsetup.size_z = 1024;
  simsetup.dimension = 3;
  simstate.MPIcoords = (int **)malloc(simstate.MPIsize * sizeof(int *));
  for (int i = 0; i < simstate.MPIsize; i++) {
    simstate.MPIcoords[i] = (int *)malloc(3 * sizeof(int));
    simstate.MPIcoords[i][0] = i%3;
    simstate.MPIcoords[i][1] = (i/3)%3;
    simstate.MPIcoords[i][2] = (i/9)%3;
  }

  prepareTestEnvironment(&simsetup);
  computeGridSize(&simsetup,&simstate, &grid);
  allocateGrid(&simsetup, &simstate, &grid);

  return 0;

}