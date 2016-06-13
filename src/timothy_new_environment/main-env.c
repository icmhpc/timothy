//
// Created by Grzegorz Bokota on 10.06.2016.
//
#include<mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "fields.h"
#include "environment.h"

void prepareTestEnvironment(struct settings *set) {
  size_t i;
  set->environments =
      calloc(set->numberOfEnvironments, sizeof(struct environment));
  for (i = 0; i < set->numberOfEnvironments; i++) {
    set->environments[i].diffusion_coefficient = 1.82e-5;
    set->environments[i].lambda_delay = 0.0;
    set->environments[i].boundary_condition = ((double)i + 1.0) * 0.1575e-6;
    set->environments[i].initial_condition_mean = ((double)i + 1.0) * 0.1575e-6;
    set->environments[i].initial_condition_variance = 0.1575e-6;
  }
}


/* zmienne globalne z timothy-ego */
struct gridData grid;
struct settings simsetup;
struct state simstate;

extern double **envField;

// int nnutrients;

double minim_arr(double * arr, int size){
  double res = arr[0];
  for (int i = 1; i < size ; ++i) {
    if(arr[i] < res)
      res = arr[i];
  }
  return res;
};

double maxim_arr(double * arr, int size){
  double res = arr[0];
  for (int i = 1; i < size ; ++i) {
    if(arr[i] > res)
      res = arr[i];
  }
  return res;
};

int main(int argc, char *argv[]) {

  int i;
  double t0, t1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &simstate.MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &simstate.MPIrank);
  if (simstate.MPIrank == 0) printf("rank: %d\n", simstate.MPIrank);
  if (simstate.MPIrank == 0) printf("size: %d\n", simstate.MPIsize);

  simsetup.gfH = 128.0;
  simsetup.maxCellsPerProc = 300000;
  simsetup.size_x = 1024;
  simsetup.size_y = 1024;
  simsetup.size_z = 1024;
  simsetup.dimension = 3;

  simstate.step = 1;

  simsetup.numberOfEnvironments = atoi(argv[1]);
  simsetup.gfH = atof(argv[2]);

  /* to jest robione w init.c */
  int periods[3];
  int reorder;
  MPI_Dims_create(simstate.MPIsize, simsetup.dimension, simstate.MPIdim);
  simstate.MPIreverseCords = calloc(simstate.MPIsize, sizeof(int));
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, simsetup.dimension, simstate.MPIdim, periods,
                  reorder, &(simstate.MPI_CART_COMM));
  simstate.MPIcoords = (int **) malloc(simstate.MPIsize * sizeof(int *));
  for (i = 0; i < simstate.MPIsize; i++) {
    simstate.MPIcoords[i] = (int *) malloc(3 * sizeof(int));
    MPI_Cart_coords(simstate.MPI_CART_COMM, i, simsetup.dimension,
                    simstate.MPIcoords[i]);
    struct intVector3d pos;
    pos.x = simstate.MPIcoords[i][0];
    pos.y = simstate.MPIcoords[i][1];
    pos.z = simstate.MPIcoords[i][2];
    *getProcesNum(&simstate, pos) = i;

  }
  /* koniec */
  if (simstate.MPIrank == 0) printf("rank2: %d\n", simstate.MPIrank);

  prepareTestEnvironment(&simsetup);
  if (simstate.MPIrank == 0) printf("rank3: %d\n", simstate.MPIrank);
  prepareEnvironment(&simsetup, &simstate, &grid);
  if (simstate.MPIrank == 0) printf("rank4: %d\n", simstate.MPIrank);
  int grid_size = grid.local_size.x * grid.local_size.y * grid.local_size.z;
  double * min_arr = calloc(simsetup.numberOfEnvironments, sizeof(double));
  double * max_arr = calloc(simsetup.numberOfEnvironments, sizeof(double));
  double * min_arr2 = calloc(simsetup.numberOfEnvironments, sizeof(double));
  double * max_arr2 = calloc(simsetup.numberOfEnvironments, sizeof(double));
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for (int i = 0; i < 1000; i++) {
    envCalculate(&simsetup, &simstate, &grid);
    for (size_t j = 0; j < simsetup.numberOfEnvironments; ++j) {
      min_arr[j] = minim_arr(envField[j], grid_size) * 1e8;
      max_arr[j] = maxim_arr(envField[j], grid_size) * 1e8;
    }
    if (i%10 == 0) {
      MPI_Allreduce(min_arr, min_arr2, simsetup.numberOfEnvironments, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(max_arr, max_arr2, simsetup.numberOfEnvironments, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (simstate.MPIrank == 0) {
        printf("Stat:");
        for (size_t j = 0; j < simsetup.numberOfEnvironments; ++j) {
          printf(" |%f, %f, %f ", min_arr2[j], max_arr2[j], max_arr2[j] - min_arr2[j]);
        }
        printf("\n");
      }
    }

  }
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (simstate.MPIrank == 0)
    printf("TIME: %f\n", t1 - t0);
  freeFieldGradient();

  MPI_Finalize();
}
