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

/* B0 definicje typów */
struct doubleVector3d {
  double x;
  double y;
  double z;
};

typedef struct _doubleVector3d {
  double x;
  double y;
  double z;
} doubleVector3d;

struct int64Vector3d {
  int64_t x;
  int64_t y;
  int64_t z;
};
/* koniec B0 */

int nnutrients;

/* B1 zmienne wejściowe Timothy */
float gfH=128.0;
double *fieldConsumption; /* units - mol (cell s)^-1 */
double *fieldProduction;  /* units - mol (cell s)^-1 */
double envLambda = 0.25;
double csize=0.2735;
double csizeInUnits=10.0;
int maxCellsPerProc=300000;
int step=1;
int bvsim=0;
int nx=1024;
int ny=1024;
int nz=1024;
int sdim=3;
/* koniec B1 */

struct environment *nutrient;

/* B2 zmienne MPI */
MPI_Comm MPI_CART_COMM;
int **MPIcoords;
int MPIrank,MPIsize;
int MPIdim[3];
/* koniec B2 */

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

double *fieldDt;

doubleVector3d lowerGridCorner,upperGridCorner;
struct int64Vector3d gridSize;
struct int64Vector3d *gridStartIdx,*gridEndIdx;
doubleVector3d *gridBuffer;
double gridResolution;
double boxVolume;
doubleVector3d globalGridSize;
int64_t gridI,gridJ,gridK;
#define grid(i,j,k) (gridBuffer[gridSize.y*gridSize.z*i+gridSize.z*j+k])

double **fieldAddr;
double **envField;

/*!
 * This function computes the sizes of the grid.
 */
void computeGridSize()
{
  gridResolution = csize * gfH;
  boxVolume = pow((gridResolution / csize) * csizeInUnits * 0.0001, 3);
  lowerGridCorner.x = 0.0;
  lowerGridCorner.y = 0.0;
  lowerGridCorner.z = 0.0;
  upperGridCorner.x = (double) nx - 1;
  upperGridCorner.y = (double) ny - 1;
  upperGridCorner.z = (double) nz - 1;


  globalGridSize.x = upperGridCorner.x - lowerGridCorner.x + 1;
  globalGridSize.y = upperGridCorner.y - lowerGridCorner.y + 1;
  globalGridSize.z = upperGridCorner.z - lowerGridCorner.z + 1;

  gridI = (int64_t) ((globalGridSize.x + 1) / gridResolution);
  gridJ = (int64_t) ((globalGridSize.y + 1) / gridResolution);
  if (sdim == 3)
    gridK = (int64_t) ((globalGridSize.z + 1) / gridResolution);
  else
    gridK = 0;

  gridI = gridI + (MPIdim[0] - gridI % MPIdim[0]);
  gridJ = gridJ + (MPIdim[1] - gridJ % MPIdim[1]);
  if (sdim == 3)
    gridK = gridK + (MPIdim[2] - gridK % MPIdim[2]);

  gridSize.x = gridI / MPIdim[0];
  gridSize.y = gridJ / MPIdim[1];
  if (sdim == 3)
    gridSize.z = gridK / MPIdim[2];
  else
    gridSize.z = 1;

  if (MPIrank == 0)
    printf("Grid size:%" PRId64 "x%" PRId64 "x%" PRId64 "\n", gridI, gridJ,
           gridK);
}

/*!
 * This function allocates grid arrays.
 */
void allocateGrid()
{
  int i, j, k;

  gridStartIdx =
      (struct int64Vector3d *) malloc(MPIsize *
                                      sizeof(struct int64Vector3d));
  gridEndIdx =
      (struct int64Vector3d *) malloc(MPIsize *
                                      sizeof(struct int64Vector3d));

  for (i = 0; i < MPIsize; i++) {
    gridStartIdx[i].x = gridSize.x * MPIcoords[i][0];
    gridStartIdx[i].y = gridSize.y * MPIcoords[i][1];
    if (sdim == 3)
      gridStartIdx[i].z = gridSize.z * MPIcoords[i][2];
    else
      gridStartIdx[i].z = 0;
    gridEndIdx[i].x = gridStartIdx[i].x + gridSize.x - 1;
    gridEndIdx[i].y = gridStartIdx[i].y + gridSize.y - 1;
    if (sdim == 3)
      gridEndIdx[i].z = gridStartIdx[i].z + gridSize.z - 1;
    else
      gridEndIdx[i].z = 0;
  }

  if (!
      (gridBuffer =
       (doubleVector3d *) calloc(gridSize.x * gridSize.y *
                                        gridSize.z,
                                        sizeof(doubleVector3d))))
    exit(1);
    //stopRun(106, "gridBuffer", __FILE__, __LINE__);

  for (i = 0; i < gridSize.x; i++)
    for (j = 0; j < gridSize.y; j++)
      for (k = 0; k < gridSize.z; k++) {
        grid(i, j, k).x =
            lowerGridCorner.x + gridResolution * (gridStartIdx[MPIrank].x +
                                                  i);
        grid(i, j, k).y =
            lowerGridCorner.y + gridResolution * (gridStartIdx[MPIrank].y +
                                                  j);
        grid(i, j, k).z =
            lowerGridCorner.z + gridResolution * (gridStartIdx[MPIrank].z +
                                                  k);
      }
  printf("Grid allocation finished.\n");
}

/*!
 * This function sets boundary conditions for domain faces.
 */
void envSetBoundary(int coord, int boundary)
{
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

void envInit() {
  envZ=(double*)malloc(nnutrients*sizeof(double));
  envStencil=(HYPRE_SStructStencil*)malloc(nnutrients*sizeof(HYPRE_SStructStencil));
  envVarTypes=(HYPRE_SStructVariable*)malloc(nnutrients*sizeof(HYPRE_SStructVariable));
}

/*!
 * This function initializes grid, stencil and matrix for a given envical field.
 */
void envInitSystem()
{
  int i, j, k, c;
  int entry;
  int var;
  HYPRE_Int offsets[7][3] = {
    {0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, 
    {0, 1, 0}, { 0, 0,-1}, {0, 0, 1}
  };

  double gridResolutionInUnits;	/* grid resolution in centimeters */

  numberOfIters = 1;
  envIter=0;
  for(var=0;var<nnutrients;var++) {
    envVarTypes[var]=HYPRE_SSTRUCT_VARIABLE_NODE;
    dt[var]=fieldDt[var];
  }

  /* 1. INIT GRID */

  /* create an empty 3D grid object */
  HYPRE_SStructGridCreate(MPI_CART_COMM, 3, 1, &envGrid);

  /* set this process box */
  envLower[0] = gridStartIdx[MPIrank].x;
  envLower[1] = gridStartIdx[MPIrank].y;
  envLower[2] = gridStartIdx[MPIrank].z;

  envUpper[0] = gridEndIdx[MPIrank].x;
  envUpper[1] = gridEndIdx[MPIrank].y;
  envUpper[2] = gridEndIdx[MPIrank].z;

  /* add a new box to the grid */
  HYPRE_SStructGridSetExtents(envGrid, 0, envLower, envUpper);

  HYPRE_SStructGridSetVariables(envGrid, 0, nnutrients, envVarTypes);
  HYPRE_SStructGridAssemble(envGrid);

  //  2. INIT STENCIL
  //HYPRE_SStructStencilCreate(3, 7, &envStencil);
  for(var = 0; var < nnutrients; var++) {
    HYPRE_SStructStencilCreate(3, 7, &envStencil[var]);
    for(entry = 0; entry < 7; entry++)
      HYPRE_SStructStencilSetEntry(envStencil[var], entry, offsets[entry],var);

  }

  // 3. SET UP THE GRAPH
  // assumption - all stencils are the same
  envObjectType = HYPRE_PARCSR;
  HYPRE_SStructGraphCreate(MPI_CART_COMM, envGrid, &envGraph);
  HYPRE_SStructGraphSetObjectType(envGraph, envObjectType);
  for(var=0;var<nnutrients;var++) 
    HYPRE_SStructGraphSetStencil(envGraph, 0, var, envStencil[var]);
  HYPRE_SStructGraphAssemble(envGraph);

  // 4. SET UP MATRIX
  long long nentries = 7;
  long long nvalues;
  double *values;
  HYPRE_Int stencil_indices[7];

  nvalues = nentries * gridSize.x * gridSize.y * gridSize.z;
  // create an empty matrix object
  HYPRE_SStructMatrixCreate(MPI_CART_COMM, envGraph, &A);
  HYPRE_SStructMatrixSetObjectType(A, envObjectType);
  // indicate that the matrix coefficients are ready to be set
  HYPRE_SStructMatrixInitialize(A);

  values = calloc(nvalues, sizeof(double));

  for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

  gridResolutionInUnits = (gridResolution / csize) * csizeInUnits * 0.0001;

  for(var=0;var<nnutrients;var++) {

    envZ[var] =
      nutrient[var].diffusion_coefficient * dt[var] / (gridResolutionInUnits *
                                      gridResolutionInUnits);

    // set the standard stencil at each grid point,
    //  we will fix the boundaries later
    for (i = 0; i < nvalues; i += nentries) {
      values[i] = 1 + 6.0 * envZ[var] + nutrient[var].lambda_delay * dt[var];
      for (j = 1; j < nentries; j++)
        values[i + j] = -envZ[var];
    }

    HYPRE_SStructMatrixSetBoxValues(A, 0, envLower, envUpper, var,
                                    nentries, stencil_indices, values);
  }

  free(values);
}

/*!
 * This function computes cell production/consumption function based on
 * the interpolated cell density field.
 */
void envCellPC(double *envPC, int var)
{
  int ch __attribute__((unused)) , i, j, k;

  if (step == 0)
    return;

  int idx = 0;
  for (k = 0; k < gridSize.z; k++)
    for (j = 0; j < gridSize.y; j++)
      for (i = 0; i < gridSize.x; i++, idx++) {
        envPC[idx] = -fieldConsumption[var] * 0.1 * dt[var];
        /*envPC[idx] = -fieldConsumption[nch] * tissueField[gridSize.z * gridSize.y * i + gridSize.z * j + k] * dt[nch];
        if(bvsim) envPC[idx]+=fieldProduction[nch] * vesselField[gridSize.z * gridSize.y * i + gridSize.z * j +k] * dt[nch] ;	*/ 
        //*(cellVolume/boxVolume);//*(1.0/cellVolume);//*dt[nch];//*dt[nch]; 
      }
}

/*!
 * This function initializes boundary conditions for a given envical field.
 */
void envInitBC()
{
  int i, j, k;
  int mi;
  int var;
  int nentries = 1;
  HYPRE_Int stencil_indices[1];
  long long nvalues = gridSize.x * gridSize.y * gridSize.z;
  double *values, *bvalues;
  double *envPC;

  envPC = (double *) calloc(nvalues, sizeof(double));
  values = calloc(nvalues, sizeof(double));
  bvalues = calloc(nvalues, sizeof(double));


  /* 5. SETUP STRUCT VECTORS FOR B AND X */

  /* create an empty vector object */
  HYPRE_SStructVectorCreate(MPI_CART_COMM, envGrid, &b);
  HYPRE_SStructVectorCreate(MPI_CART_COMM, envGrid, &x);

  /* as with the matrix, set the appropriate object type for the vectors */
  HYPRE_SStructVectorSetObjectType(b, envObjectType);
  HYPRE_SStructVectorSetObjectType(x, envObjectType);

  /* indicate that the vector coefficients are ready to be set */
  HYPRE_SStructVectorInitialize(b);
  HYPRE_SStructVectorInitialize(x);

  for(var=0;var<nnutrients;var++) {

  envCellPC(envPC,var);

  /* set the values */
  mi = 0;
  for (k = envLower[2]; k <= envUpper[2]; k++)
    for (j = envLower[1]; j <= envUpper[1]; j++)
      for (i = envLower[0]; i <= envUpper[0]; i++) {
        values[mi] = nutrient[var].initial_condition_mean;
        mi++;
      }

  HYPRE_SStructVectorSetBoxValues(b, 0, envLower, envUpper, var,
                                  values);

  mi = 0;
  for (k = envLower[2]; k <= envUpper[2]; k++)
    for (j = envLower[1]; j <= envUpper[1]; j++)
      for (i = envLower[0]; i <= envUpper[0]; i++) {
        values[mi] = nutrient[var].initial_condition_mean;
        mi++;
      }

  HYPRE_SStructVectorSetBoxValues(x, 0, envLower, envUpper, var,
                                  values);

  /* incorporate boundary conditions; Dirichlet on 6 faces */

  for (i = 0; i < nvalues; i++)
    values[i] = envZ[var];
  for (i = 0; i < nvalues; i++)
    bvalues[i] = envZ[var] * nutrient[var].boundary_condition;

  if (MPIcoords[MPIrank][0] == 0) {
    nvalues = nentries * gridSize.y * gridSize.z;
    envSetBoundary(0, -1);
    stencil_indices[0] = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                      nentries, stencil_indices, values);
    HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                      bvalues);
  }
  if (MPIcoords[MPIrank][0] == MPIdim[0] - 1) {
    nvalues = nentries * gridSize.y * gridSize.z;
    envSetBoundary(0, 1);
    stencil_indices[0] = 2;
    HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                      nentries, stencil_indices, values);
    HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                      bvalues);
  }
  if (MPIcoords[MPIrank][1] == 0) {
    nvalues = nentries * gridSize.x * gridSize.z;
    envSetBoundary(1, -1);
    stencil_indices[0] = 3;
    HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                      nentries, stencil_indices, values);
    HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                      bvalues);
  }
  if (MPIcoords[MPIrank][1] == MPIdim[1] - 1) {
    nvalues = nentries * gridSize.x * gridSize.z;
    envSetBoundary(1, 1);
    stencil_indices[0] = 4;
    HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                      nentries, stencil_indices, values);
    HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                      bvalues);
  }
  if (MPIcoords[MPIrank][2] == 0) {
    nvalues = nentries * gridSize.x * gridSize.y;
    envSetBoundary(2, -1);
    stencil_indices[0] = 5;
    HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                      nentries, stencil_indices, values);
    HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                      bvalues);
  }
  if (MPIcoords[MPIrank][2] == MPIdim[2] - 1) {
    nvalues = nentries * gridSize.x * gridSize.y;
    envSetBoundary(2, 1);
    stencil_indices[0] = 6;
    HYPRE_SStructMatrixAddToBoxValues(A, 0, bcLower, bcUpper, var,
                                      nentries, stencil_indices, values);
    HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper, var,
                                      bvalues);
  }

  /* add production consumption function to the right side */
  HYPRE_SStructVectorAddToBoxValues(b, 0, envLower, envUpper, var,
                                    envPC);

  }

  free(envPC);
  free(values);
  free(bvalues);
  /* stdout brought back */
}

/*!
 * This function initializes Hypre for solving a given envical field.
 */
void envInitSolver()
{

  HYPRE_SStructMatrixAssemble(A);
  /* This is a collective call finalizing the vector assembly.
     The vector is now ``ready to be used'' */
  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);

  HYPRE_SStructMatrixGetObject(A, (void **) &parA);
  HYPRE_SStructVectorGetObject(b, (void **) &parb);
  HYPRE_SStructVectorGetObject(x, (void **) &parx);

  HYPRE_ParCSRPCGCreate(MPI_CART_COMM, &envSolver);
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
  HYPRE_ParCSRPCGSetup(envSolver, parA, parb,
                       parx);

}

/*!
 * This is a driving function for solving next time step
 * of a given envical field.
 */
void envSolve()
{
  int i, j, k;
  int idx;
  int var;
  double *values;
  int stepIter = 0;
  long long nvalues = gridSize.x * gridSize.y * gridSize.z;
  double *envPC;
  if (MPIrank == 0 ) {
    printf("Solving field.");
    fflush(stdout);
  }

  values = (double *) calloc(nvalues, sizeof(double));
  envPC = (double *) calloc(nvalues, sizeof(double));

  while (stepIter < numberOfIters) {
    if (envIter > 0) {
      /* update right hand side */
      for(var=0;var<nnutrients;var++) {

      envCellPC(envPC,var);
      HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper,
                                      var, values);
      HYPRE_SStructVectorSetBoxValues(b, 0, envLower, envUpper,
                                      var, values);
      for (i = 0; i < nvalues; i++)
        values[i] = envZ[var] * nutrient[var].boundary_condition;
      if (MPIcoords[MPIrank][0] == 0) {
        envSetBoundary(0, -1);
        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                          var, values);
      }
      if (MPIcoords[MPIrank][0] == MPIdim[0] - 1) {
        envSetBoundary(0, 1);
        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                          var, values);
      }
      if (MPIcoords[MPIrank][1] == 0) {
        envSetBoundary(1, -1);
        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                          var, values);
      }
      if (MPIcoords[MPIrank][1] == MPIdim[1] - 1) {
        envSetBoundary(1, 1);
        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                          var, values);
      }
      if (MPIcoords[MPIrank][2] == 0) {
        envSetBoundary(2, -1);
        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                          var, values);
      }
      if (MPIcoords[MPIrank][2] == MPIdim[2] - 1) {
        envSetBoundary(2, 1);
        HYPRE_SStructVectorAddToBoxValues(b, 0, bcLower, bcUpper,
                                          var, values);
      }
      HYPRE_SStructVectorAddToBoxValues(b, 0, envLower,
                                        envUpper, var, envPC);
      HYPRE_SStructVectorAssemble(b);
      HYPRE_SStructVectorAssemble(x);
    }

    }

    HYPRE_ParCSRPCGSolve(envSolver, parA, parb,
                         parx);

    for(var=0;var<nnutrients;var++) {
      HYPRE_SStructVectorGather(x);
      HYPRE_SStructVectorGetBoxValues(x, 0, envLower, envUpper, var,
                                    values);
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

void allocateFields() {

  int i,nf;

  dt=(double*)malloc(nnutrients*sizeof(double));
  fieldDt=(double*)malloc(nnutrients*sizeof(double));
  fieldAddr=(double**)malloc(nnutrients*sizeof(double*));
  envField=(double**)malloc(nnutrients*sizeof(double*));

  fieldConsumption=(double*)malloc(nnutrients*sizeof(double));
  fieldProduction=(double*)malloc(nnutrients*sizeof(double));
  nutrient=(struct environment*)malloc(nnutrients*sizeof(struct environment));

  for(nf=0;nf<nnutrients;nf++) {  
    fieldAddr[nf] =
        (double *) calloc(gridSize.x * gridSize.y * gridSize.z,
                          sizeof(double));
    envField[nf] = (double *) fieldAddr[nf];
    for (i = 0; i < gridSize.x * gridSize.y * gridSize.z; i++)
      envField[nf][i] = 0.0;
  }

  printf("Fields allocation finished\n");
}

/* ustawia testowe parametry dla fields */
void initFields() {

  int i,j,k;
  int var;

  for(var=0;var<nnutrients;var++) {
    nutrient[var].diffusion_coefficient=1.82e-5;
    nutrient[var].lambda_delay=0.0;
    nutrient[var].boundary_condition=((double)var+1.0)*0.1575e-6; 
    nutrient[var].initial_condition_mean=((double)var+1.0)*0.1575e-6; 
    nutrient[var].initial_condition_variance=0.0; 
    fieldConsumption[var]=8.3e-17; 
    fieldProduction[var]=8.3e-1;  
    fieldDt[var]=4000;
    for (k = 0; k < gridSize.z; k++)
      for (j = 0; j < gridSize.y; j++)
        for (i = 0; i < gridSize.x; i++) {
          envField[var][gridSize.y * gridSize.z * i + gridSize.z * j + k] =
            nutrient[var].initial_condition_mean;
        }
  }

}

void initEnvironment(struct state * st, struct settings *set){

}


int main(int argc,char* argv[]) {

  int i;
  int periods[3];
  int reorder;
  double t0,t1;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  nnutrients=atoi(argv[1]);
  gfH=atof(argv[2]);

  MPI_Dims_create(MPIsize, sdim, MPIdim);
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, sdim, MPIdim, periods, reorder,
                  &MPI_CART_COMM);
  MPIcoords = (int **) malloc(MPIsize * sizeof(int *));
  for (i = 0; i < MPIsize; i++) {
    MPIcoords[i] = (int *) malloc(3 * sizeof(int));
    MPI_Cart_coords(MPI_CART_COMM, i, sdim, MPIcoords[i]);
  }

  computeGridSize();
  allocateGrid();
  allocateFields();
  initFields();
  envInit();
  
  MPI_Barrier(MPI_COMM_WORLD);
  t0=MPI_Wtime();
  envInitSystem();
  envInitBC();
  envInitSolver();
  envSolve();
  MPI_Barrier(MPI_COMM_WORLD);
  t1=MPI_Wtime();
  if(MPIrank==0) printf("TIME: %f\n",t1-t0);

  MPI_Finalize();

}

