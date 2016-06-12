#ifndef TIMOTHY_ENVIRONMENT_H
#define TIMOTHY_ENVIRONMENT_H
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <mpi.h>
#include "global_structures.h"
#ifdef __MIC__
#define MIC_ATTR __attribute__((target(mic)))
#else
#define MIC_ATTR
#endif

/* NEW */


/* NEW */
struct gridData{
  struct int64Vector3d global_size;
  doubleVector3d lower_corner,upper_corner;
  struct intVector3d local_size;
  double resolution;
  double voxel_volume;
  struct int64Vector3d *lower_indices,*upper_indices;
  doubleVector3d *buffer;
};
/* NEW */

void envCalculate(const struct settings *simsetup, struct state *simstate,
                  struct gridData *grid);

void prepareEnvironment(const struct settings *simsetup, struct state *simstate,
                        struct gridData *grid);


#endif //TIMOTHY_ENVIRONMENT_H
