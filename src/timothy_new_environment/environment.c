/* **************************************************************************
 * This file is part of Timothy
 *
 * Copyright (c) 2014/15 Maciej Cytowski
 * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * *************************************************************************/

#include <float.h>
#include <inttypes.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "environment.h"
#include "fields.h"
#include "mpi.h"

/* zmienne lokalne w pliku environment */
HYPRE_SStructGrid envGrid;
HYPRE_SStructGraph envGraph;
HYPRE_SStructStencil *envStencil;
HYPRE_SStructMatrix A;
HYPRE_SStructVector b;
HYPRE_SStructVector x;
HYPRE_ParCSRMatrix parA;
HYPRE_ParVector parb;
HYPRE_ParVector parx;
HYPRE_Solver envSolver;
HYPRE_Solver envPrecond;
int envObjectType;
HYPRE_Int envLower[3], envUpper[3];
HYPRE_Int bcLower[3];
HYPRE_Int bcUpper[3];
double *dt;
int envIter;
double *envZ;
HYPRE_SStructVariable *envVarTypes;
int numberOfIters;
/* koniec */

/* do polaczenia z fields */
double *fieldDt;
double **fieldAddr;
double **envField;
double *fieldConsumption; /* units - mol (cell s)^-1 */
double *fieldProduction;  /* units - mol (cell s)^-1 */
/* koniec */

/* tych zmiennych nie odnalazlem w strukturach */
double csize = 0.2735;
double csizeInUnits = 10.0;
int bvsim = 0;
/* koniec */

/*!
 * This function computes the sizes of the grid.
 */
void computeGridSize(const struct settings *simsetup,
                     const struct state *simstate, struct gridData *grid) {
  doubleVector3d globalGridSize;

  grid->resolution = csize * simsetup->gfH;
  grid->voxel_volume =
      pow((grid->resolution / csize) * csizeInUnits * 0.0001, 3);
  grid->lower_corner.x = 0.0;
  grid->lower_corner.y = 0.0;
  grid->lower_corner.z = 0.0;
  grid->upper_corner.x = (double)simsetup->size_x - 1;
  grid->upper_corner.y = (double)simsetup->size_y - 1;
  grid->upper_corner.z = (double)simsetup->size_z - 1;

  globalGridSize.x = grid->upper_corner.x - grid->lower_corner.x + 1;
  globalGridSize.y = grid->upper_corner.y - grid->lower_corner.y + 1;
  globalGridSize.z = grid->upper_corner.z - grid->lower_corner.z + 1;

  grid->global_size.x = (int64_t)((globalGridSize.x + 1) / grid->resolution);
  grid->global_size.y = (int64_t)((globalGridSize.y + 1) / grid->resolution);
  if (simsetup->dimension == 3)
    grid->global_size.z = (int64_t)((globalGridSize.z + 1) / grid->resolution);
  else
    grid->global_size.z = 0;

  grid->global_size.x =
      grid->global_size.x +
      (simstate->MPIdim[0] - grid->global_size.x % simstate->MPIdim[0]);
  grid->global_size.y =
      grid->global_size.y +
      (simstate->MPIdim[1] - grid->global_size.y % simstate->MPIdim[1]);
  if (simsetup->dimension == 3)
    grid->global_size.z =
        grid->global_size.z +
        (simstate->MPIdim[2] - grid->global_size.z % simstate->MPIdim[2]);

  grid->local_size.x = grid->global_size.x / simstate->MPIdim[0];
  grid->local_size.y = grid->global_size.y / simstate->MPIdim[1];
  if (simsetup->dimension == 3)
    grid->local_size.z = grid->global_size.z / simstate->MPIdim[2];
  else
    grid->local_size.z = 1;

  if (simstate->MPIrank == 0)
    printf("Grid size:%" PRId64 "x%" PRId64 "x%" PRId64 "\n",
           grid->global_size.x, grid->global_size.y, grid->global_size.z);
}

/*!
 * This function allocates grid arrays.
 */
void allocateGrid(const struct settings *simsetup, const struct state *simstate,
                  struct gridData *grid) {
  int i, j, k;
#define grid_node(i, j, k)                                                     \
  (grid->buffer[grid->local_size.y * grid->local_size.z * i +                  \
                grid->local_size.z * j + k])
  grid->lower_indices = (struct int64Vector3d *)malloc(
      simstate->MPIsize * sizeof(struct int64Vector3d));
  grid->upper_indices = (struct int64Vector3d *)malloc(
      simstate->MPIsize * sizeof(struct int64Vector3d));

  for (i = 0; i < simstate->MPIsize; i++) {
    grid->lower_indices[i].x = grid->local_size.x * simstate->MPIcoords[i][0];
    grid->lower_indices[i].y = grid->local_size.y * simstate->MPIcoords[i][1];
    if (simsetup->dimension == 3)
      grid->lower_indices[i].z = grid->local_size.z * simstate->MPIcoords[i][2];
    else
      grid->lower_indices[i].z = 0;
    grid->upper_indices[i].x =
        grid->lower_indices[i].x + grid->local_size.x - 1;
    grid->upper_indices[i].y =
        grid->lower_indices[i].y + grid->local_size.y - 1;
    if (simsetup->dimension == 3)
      grid->upper_indices[i].z =
          grid->lower_indices[i].z + grid->local_size.z - 1;
    else
      grid->upper_indices[i].z = 0;
  }

  if (!(grid->buffer = (doubleVector3d *)calloc((size_t)grid->local_size.x *
                                                    grid->local_size.y *
                                                    grid->local_size.z,
                                                sizeof(doubleVector3d))))
    exit(1);
  // stopRun(106, "gridBuffer", __FILE__, __LINE__);

  for (i = 0; i < grid->local_size.x; i++)
    for (j = 0; j < grid->local_size.y; j++)
      for (k = 0; k < grid->local_size.z; k++) {
        grid_node(i, j, k).x =
            grid->lower_corner.x +
            grid->resolution * (grid->lower_indices[simstate->MPIrank].x + i);
        grid_node(i, j, k).y =
            grid->lower_corner.y +
            grid->resolution * (grid->lower_indices[simstate->MPIrank].y + j);
        grid_node(i, j, k).z =
            grid->lower_corner.z +
            grid->resolution * (grid->lower_indices[simstate->MPIrank].z + k);
      }
  printf("Grid allocation finished.\n");
#undef grid_node
}

/*!
 * This function sets boundary conditions for domain faces.
 */
void envSetBoundary(int coord, int boundary) {
  if (coord == 0 && boundary == -1) {
    bcLower[0] = envLower[0];
    bcUpper[0] = envLower[0];
    bcLower[1] = envLower[1];
    bcUpper[1] = envUpper[1];
    bcLower[2] = envLower[2];
    bcUpper[2] = envUpper[2];
  }
  if (coord == 0 && boundary == 1) {
    bcLower[0] = envUpper[0];
    bcUpper[0] = envUpper[0];
    bcLower[1] = envLower[1];
    bcUpper[1] = envUpper[1];
    bcLower[2] = envLower[2];
    bcUpper[2] = envUpper[2];
  }
  if (coord == 1 && boundary == -1) {
    bcLower[0] = envLower[0];
    bcUpper[0] = envUpper[0];
    bcLower[1] = envLower[1];
    bcUpper[1] = envLower[1];
    bcLower[2] = envLower[2];
    bcUpper[2] = envUpper[2];
  }
  if (coord == 1 && boundary == 1) {
    bcLower[0] = envLower[0];
    bcUpper[0] = envUpper[0];
    bcLower[1] = envUpper[1];
    bcUpper[1] = envUpper[1];
    bcLower[2] = envLower[2];
    bcUpper[2] = envUpper[2];
  }
  if (coord == 2 && boundary == -1) {
    bcLower[0] = envLower[0];
    bcUpper[0] = envUpper[0];
    bcLower[1] = envLower[1];
    bcUpper[1] = envUpper[1];
    bcLower[2] = envLower[2];
    bcUpper[2] = envLower[2];
  }
  if (coord == 2 && boundary == 1) {
    bcLower[0] = envLower[0];
    bcUpper[0] = envUpper[0];
    bcLower[1] = envLower[1];
    bcUpper[1] = envUpper[1];
    bcLower[2] = envUpper[2];
    bcUpper[2] = envUpper[2];
  }
}

void envInit(size_t numberOfEnvironments) {
  envZ = (double *)malloc(numberOfEnvironments * sizeof(double));
  envStencil =
      (HYPRE_SStructStencil *)malloc(numberOfEnvironments * sizeof(HYPRE_SStructStencil));
  envVarTypes = (HYPRE_SStructVariable *)malloc(numberOfEnvironments *
                                                sizeof(HYPRE_SStructVariable));
}

/*!
 * This function initializes grid, stencil and matrix for a given envical field.
 */
void envInitSystem(const struct state *simstate, const struct settings *set,
                   struct gridData *grid) {
  size_t i;
  int j;//, k, c;
  int entry;
  size_t var;
  HYPRE_Int offsets[7][3] = {{0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                             {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};

  double gridResolutionInUnits; /* grid resolution in centimeters */

  numberOfIters = 1;
  envIter = 0;
  for (var = 0; var < set->numberOfEnvironments; var++) {
    envVarTypes[var] = HYPRE_SSTRUCT_VARIABLE_NODE;
    dt[var] = fieldDt[var];
  }

  /* 1. INIT GRID */

  /* create an empty 3D grid object */
  HYPRE_SStructGridCreate(simstate->MPI_CART_COMM, 3, 1, &envGrid);

  /* set this process box */
  envLower[0] = grid->lower_indices[simstate->MPIrank].x;
  envLower[1] = grid->lower_indices[simstate->MPIrank].y;
  envLower[2] = grid->lower_indices[simstate->MPIrank].z;

  envUpper[0] = grid->upper_indices[simstate->MPIrank].x;
  envUpper[1] = grid->upper_indices[simstate->MPIrank].y;
  envUpper[2] = grid->upper_indices[simstate->MPIrank].z;

  /* add a new box to the grid */
  HYPRE_SStructGridSetExtents(envGrid, 0, envLower, envUpper);

  HYPRE_SStructGridSetVariables(
      envGrid, 0, (HYPRE_Int)set->numberOfEnvironments, envVarTypes);
  HYPRE_SStructGridAssemble(envGrid);

  //  2. INIT STENCIL
  // HYPRE_SStructStencilCreate(3, 7, &envStencil);
  for (var = 0; var < set->numberOfEnvironments; var++) {
    HYPRE_SStructStencilCreate(3, 7, &envStencil[var]);
    for (entry = 0; entry < 7; entry++)
      HYPRE_SStructStencilSetEntry(envStencil[var], entry, offsets[entry], var);
  }

  // 3. SET UP THE GRAPH
  // assumption - all stencils are the same
  envObjectType = HYPRE_PARCSR;
  HYPRE_SStructGraphCreate(simstate->MPI_CART_COMM, envGrid, &envGraph);
  HYPRE_SStructGraphSetObjectType(envGraph, envObjectType);
  for (var = 0; var < set->numberOfEnvironments; var++)
    HYPRE_SStructGraphSetStencil(envGraph, 0, var, envStencil[var]);
  HYPRE_SStructGraphAssemble(envGraph);

  // 4. SET UP MATRIX
  long long nentries = 7;
  size_t nvalues;
  double *values;
  HYPRE_Int stencil_indices[7];

  nvalues = (size_t)nentries * grid->local_size.x * grid->local_size.y *
            grid->local_size.z;
  // create an empty matrix object
  HYPRE_SStructMatrixCreate(simstate->MPI_CART_COMM, envGraph, &A);
  HYPRE_SStructMatrixSetObjectType(A, envObjectType);
  // indicate that the matrix coefficients are ready to be set
  HYPRE_SStructMatrixInitialize(A);

  values = calloc(nvalues, sizeof(double));

  for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

  gridResolutionInUnits = (grid->resolution / csize) * csizeInUnits * 0.0001;

  for (var = 0; var < set->numberOfEnvironments; var++) {

    envZ[var] = set->environments[var].diffusion_coefficient * dt[var] /
                (gridResolutionInUnits * gridResolutionInUnits);

    // set the standard stencil at each grid point,
    //  we will fix the boundaries later
    for (i = 0; i < nvalues; i += nentries) {
      values[i] =
          1 + 6.0 * envZ[var] + set->environments[var].lambda_delay * dt[var];
      for (j = 1; j < nentries; j++)
        values[i + j] = -envZ[var];
    }

    HYPRE_SStructMatrixSetBoxValues(A, 0, envLower, envUpper, var, nentries,
                                    stencil_indices, values);
  }

  free(values);
}

/*!
 * This function computes cell production/consumption function based on
 * the interpolated cell density field.
 */
void envCellPC(const struct state *simstate, double *envPC, HYPRE_Int var,
               struct gridData *grid) {
  int ch __attribute__((unused)), i, j, k;

  if (simstate->step == 0)
    return;

  size_t idx = 0;
  for (k = 0; k < grid->local_size.z; k++)
    for (j = 0; j < grid->local_size.y; j++)
      for (i = 0; i < grid->local_size.x; i++, idx++) {
        envPC[idx] = -fieldConsumption[var] * 0.1 * dt[var];
        /*envPC[idx] = -fieldConsumption[nch] * tissueField[gridSize.z *
        gridSize.y * i + gridSize.z * j + k] * dt[nch];
        if(bvsim) envPC[idx]+=fieldProduction[nch] * vesselField[gridSize.z *
        gridSize.y * i + gridSize.z * j +k] * dt[nch] ;	*/
        //*(cellVolume/boxVolume);//*(1.0/cellVolume);//*dt[nch];//*dt[nch];
      }
}

/*!
 * This function initializes boundary conditions for a given envical field.
 */
void envInitBC(const struct state *simstate, const struct settings *set,
               struct gridData *grid) {
  HYPRE_Int i, j, k;
  int mi;
  HYPRE_Int var;
  int nentries = 1;
  HYPRE_Int stencil_indices[1];
  long long nvalues =
      grid->local_size.x * grid->local_size.y * grid->local_size.z;
  double *values, *bvalues;
  double *envPC;

  envPC = (double *)calloc(nvalues, sizeof(double));
  values = calloc(nvalues, sizeof(double));
  bvalues = calloc(nvalues, sizeof(double));

  /* 5. SETUP STRUCT VECTORS FOR B AND X */

  /* create an empty vector object */
  HYPRE_SStructVectorCreate(simstate->MPI_CART_COMM, envGrid, &b);
  HYPRE_SStructVectorCreate(simstate->MPI_CART_COMM, envGrid, &x);

  /* as with the matrix, set the appropriate object type for the vectors */
  HYPRE_SStructVectorSetObjectType(b, envObjectType);
  HYPRE_SStructVectorSetObjectType(x, envObjectType);

  /* indicate that the vector coefficients are ready to be set */
  HYPRE_SStructVectorInitialize(b);
  HYPRE_SStructVectorInitialize(x);

  for (var = 0; var < (int) set->numberOfEnvironments; var++) {

    envCellPC(simstate, envPC, var, grid);

    /* set the values */
    mi = 0;
    for (k = envLower[2]; k <= envUpper[2]; k++)
      for (j = envLower[1]; j <= envUpper[1]; j++)
        for (i = envLower[0]; i <= envUpper[0]; i++) {
          values[mi] = set->environments[var].initial_condition_mean;
          mi++;
        }

    HYPRE_SStructVectorSetBoxValues(b, 0, envLower, envUpper, var, values);

    mi = 0;
    for (k = envLower[2]; k <= envUpper[2]; k++)
      for (j = envLower[1]; j <= envUpper[1]; j++)
        for (i = envLower[0]; i <= envUpper[0]; i++) {
          values[mi] = set->environments[var].initial_condition_mean;
          mi++;
        }

    HYPRE_SStructVectorSetBoxValues(x, 0, envLower, envUpper, var, values);

    /* incorporate boundary conditions; Dirichlet on 6 faces */

    for (i = 0; i < nvalues; i++)
      values[i] = envZ[var];
    for (i = 0; i < nvalues; i++)
      bvalues[i] = envZ[var] * set->environments[var].boundary_condition;

    if (simstate->MPIcoords[simstate->MPIrank][0] == 0) {
      nvalues = nentries * grid->local_size.y * grid->local_size.z;
      envSetBoundary(0, -1);
      stencil_indices[0] = 1;
      HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var, nentries,
                                        stencil_indices, values);
      HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var, bvalues);
    }
    if (simstate->MPIcoords[simstate->MPIrank][0] == simstate->MPIdim[0] - 1) {
      nvalues = nentries * grid->local_size.y * grid->local_size.z;
      envSetBoundary(0, 1);
      stencil_indices[0] = 2;
      HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var, nentries,
                                        stencil_indices, values);
      HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var, bvalues);
    }
    if (simstate->MPIcoords[simstate->MPIrank][1] == 0) {
      nvalues = nentries * grid->local_size.x * grid->local_size.z;
      envSetBoundary(1, -1);
      stencil_indices[0] = 3;
      HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var, nentries,
                                        stencil_indices, values);
      HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var, bvalues);
    }
    if (simstate->MPIcoords[simstate->MPIrank][1] == simstate->MPIdim[1] - 1) {
      nvalues = nentries * grid->local_size.x * grid->local_size.z;
      envSetBoundary(1, 1);
      stencil_indices[0] = 4;
      HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var, nentries,
                                        stencil_indices, values);
      HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var, bvalues);
    }
    if (simstate->MPIcoords[simstate->MPIrank][2] == 0) {
      nvalues = nentries * grid->local_size.x * grid->local_size.y;
      envSetBoundary(2, -1);
      stencil_indices[0] = 5;
      HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var, nentries,
                                        stencil_indices, values);
      HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var, bvalues);
    }
    if (simstate->MPIcoords[simstate->MPIrank][2] == simstate->MPIdim[2] - 1) {
      nvalues = nentries * grid->local_size.x * grid->local_size.y;
      envSetBoundary(2, 1);
      stencil_indices[0] = 6;
      HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var, nentries,
                                        stencil_indices, values);
      HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var, bvalues);
    }

    /* add production consumption function to the right side */
    HYPRE_SStructVectorAddToBoxValues(b, 0, envLower, envUpper, var, envPC);
  }

  free(envPC);
  free(values);
  free(bvalues);
  /* stdout brought back */
}

/*!
 * This function initializes Hypre for solving a given envical field.
 */
void envInitSolver(struct state *simstate) {

  HYPRE_SStructMatrixAssemble(A);
  /* This is a collective call finalizing the vector assembly.
     The vector is now ``ready to be used'' */
  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);

  HYPRE_SStructMatrixGetObject(A, (void **)&parA);
  HYPRE_SStructVectorGetObject(b, (void **)&parb);
  HYPRE_SStructVectorGetObject(x, (void **)&parx);

  HYPRE_ParCSRPCGCreate(simstate->MPI_CART_COMM, &envSolver);
  HYPRE_ParCSRPCGSetTol(envSolver, 1.0e-12);
  HYPRE_ParCSRPCGSetPrintLevel(envSolver, 2);
  HYPRE_ParCSRPCGSetMaxIter(envSolver, 50);

  HYPRE_BoomerAMGCreate(&envPrecond);
  HYPRE_BoomerAMGSetMaxIter(envPrecond, 1);
  HYPRE_BoomerAMGSetTol(envPrecond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(envPrecond, 2);
  HYPRE_BoomerAMGSetCoarsenType(envPrecond, 6);
  HYPRE_BoomerAMGSetRelaxType(envPrecond, 6);
  HYPRE_BoomerAMGSetNumSweeps(envPrecond, 1);

  HYPRE_ParCSRPCGSetPrecond(envSolver, HYPRE_BoomerAMGSolve,
                            HYPRE_BoomerAMGSetup, envPrecond);
  HYPRE_ParCSRPCGSetup(envSolver, parA, parb, parx);
}

/*!
 * This is a driving function for solving next time step
 * of a given envical field.
 */
void envSolve(const struct state *simstate, const struct settings *simsetup,
              struct gridData *grid) {
  size_t i;//, j, k;
  //int idx;
  HYPRE_Int var;
  double *values;
  int stepIter = 0;
  size_t nvalues =
      (size_t)grid->local_size.x * grid->local_size.y * grid->local_size.z;
  double *envPC;
  if (simstate->MPIrank == 0) {
    printf("Solving field.");
    fflush(stdout);
  }

  values = (double *)calloc(nvalues, sizeof(double));
  envPC = (double *)calloc(nvalues, sizeof(double));

  while (stepIter < numberOfIters) {
    if (envIter > 0) {
      /* update right hand side */
      for (var = 0; var < (int) simsetup->numberOfEnvironments; var++) {

        envCellPC(simstate, envPC, var, grid);
        HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper, var, values);
        HYPRE_SStructVectorSetBoxValues(b, 0, envLower, envUpper, var, values);
        for (i = 0; i < nvalues; i++)
          values[i] =
              envZ[var] * simsetup->environments[var].boundary_condition;
        if (simstate->MPIcoords[simstate->MPIrank][0] == 0) {
          envSetBoundary(0, -1);
          HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                            values);
        }
        if (simstate->MPIcoords[simstate->MPIrank][0] ==
            simstate->MPIdim[0] - 1) {
          envSetBoundary(0, 1);
          HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                            values);
        }
        if (simstate->MPIcoords[simstate->MPIrank][1] == 0) {
          envSetBoundary(1, -1);
          HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                            values);
        }
        if (simstate->MPIcoords[simstate->MPIrank][1] ==
            simstate->MPIdim[1] - 1) {
          envSetBoundary(1, 1);
          HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                            values);
        }
        if (simstate->MPIcoords[simstate->MPIrank][2] == 0) {
          envSetBoundary(2, -1);
          HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                            values);
        }
        if (simstate->MPIcoords[simstate->MPIrank][2] ==
            simstate->MPIdim[2] - 1) {
          envSetBoundary(2, 1);
          HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                            values);
        }
        HYPRE_SStructVectorAddToBoxValues(b, 0, envLower, envUpper, var, envPC);
        HYPRE_SStructVectorAssemble(b);
        HYPRE_SStructVectorAssemble(x);
      }
    }

    HYPRE_ParCSRPCGSolve(envSolver, parA, parb, parx);

    for (var = 0; var < (int) simsetup->numberOfEnvironments; var++) {
      HYPRE_SStructVectorGather(x);
      HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper, var, values);
      /*         idx = 0;
         for (k = 0; k < gridSize.z; k++)
           for (j = 0; j < gridSize.y; j++)
             for (i = 0; i < gridSize.x; i++, idx++) {
               //envField[nch][gridSize.y * gridSize.z * i + gridSize.z * j +
               //               k] = values[idx];
               printf("[%d,%d,%d] %.12f\n",i,j,k,values[idx]);
             }
     */
    }

    /* copy solution to field buffer */
    /*   HYPRE_SStructVectorGather(x);
        HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper, 0,
                                        values);
        idx = 0;
        for (k = 0; k < gridSize.z; k++)
          for (j = 0; j < gridSize.y; j++)
            for (i = 0; i < gridSize.x; i++, idx++) {
              envField[nch][gridSize.y * gridSize.z * i + gridSize.z * j +
                             k] = values[idx];
            }
    */
    envIter++;
    stepIter++;
  }

  free(values);
  free(envPC);
}

void allocateFields(const struct settings *simsetup, struct gridData *grid) {

  int i;
  size_t  nf;

  dt = (double *)malloc(simsetup->numberOfEnvironments * sizeof(double));
  fieldDt = (double *)malloc(simsetup->numberOfEnvironments * sizeof(double));
  fieldAddr =
      (double **)malloc(simsetup->numberOfEnvironments * sizeof(double *));
  envField =
      (double **)malloc(simsetup->numberOfEnvironments * sizeof(double *));

  fieldConsumption =
      (double *)malloc(simsetup->numberOfEnvironments * sizeof(double));
  fieldProduction =
      (double *)malloc(simsetup->numberOfEnvironments * sizeof(double));
  for (nf = 0; nf < simsetup->numberOfEnvironments; nf++) {
    fieldAddr[nf] = (double *)calloc((size_t) grid->local_size.x * grid->local_size.y *
                                         grid->local_size.z,
                                     sizeof(double));
    envField[nf] = (double *)fieldAddr[nf];
    for (i = 0;
         i < grid->local_size.x * grid->local_size.y * grid->local_size.z; i++)
      envField[nf][i] = 0.0;
  }

  printf("Fields allocation finished\n");
}

/* ustawia testowe parametry dla fields */
void initFields(const struct settings *simsetup, struct gridData *grid) {

  int i, j, k;
  size_t var;

  for (var = 0; var < simsetup->numberOfEnvironments; var++) {

    fieldConsumption[var] = 8.3e-17;
    fieldProduction[var] = 8.3e-1;
    fieldDt[var] = 4000;
    for (k = 0; k < grid->local_size.z; k++)
      for (j = 0; j < grid->local_size.y; j++)
        for (i = 0; i < grid->local_size.x; i++) {
          envField[var][grid->local_size.y * grid->local_size.z * i +
                        grid->local_size.z * j + k] =
              simsetup->environments[var].initial_condition_mean;
        }
  }
}

/* to jest glowna funkcja "biblioteczna" */
void prepareEnvironment(const struct settings *simsetup, struct state *simstate,
                        struct gridData *grid) {
  //double t0, t1;
  computeGridSize(simsetup, simstate, grid);
  allocateGrid(simsetup, simstate, grid);
  allocateFields(simsetup, grid);
  initFields(simsetup, grid);
  envInit(simsetup->numberOfEnvironments);
  envInitSystem(simstate, simsetup, grid);
  envInitBC(simstate, simsetup, grid);
  envInitSolver(simstate);
  allocateFieldGradient(simsetup, simstate, grid);

}

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

void envCalculate (const struct settings *simsetup, struct state *simstate,
                   struct gridData *grid) {
  initFieldHaloExchange(simstate, simsetup, grid);
  envSolve(simstate, simsetup, grid);
}

/* zmienne globalne z timothy-ego */
struct gridData grid;
struct environment *nutrient;
struct settings simsetup;
struct state simstate;
// int nnutrients;

int main(int argc, char *argv[]) {

  int i;
  double t0, t1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &simstate.MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &simstate.MPIrank);

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
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, simsetup.dimension, simstate.MPIdim, periods,
                  reorder, &(simstate.MPI_CART_COMM));
  simstate.MPIcoords = (int **)malloc(simstate.MPIsize * sizeof(int *));
  for (i = 0; i < simstate.MPIsize; i++) {
    simstate.MPIcoords[i] = (int *)malloc(3 * sizeof(int));
    MPI_Cart_coords(simstate.MPI_CART_COMM, i, simsetup.dimension,
                    simstate.MPIcoords[i]);
  }
  /* koniec */
  prepareTestEnvironment(&simsetup);
  prepareEnvironment(&simsetup, &simstate, &grid);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  envCalculate(&simsetup, &simstate, &grid);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (simstate.MPIrank == 0)
    printf("TIME: %f\n", t1 - t0);
  freeFieldGradient();

  MPI_Finalize();
}
