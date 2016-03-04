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

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <inttypes.h>
#include <sprng.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include "global.h"
#include "inline.h"
#include "fields.h"
#include "utils.h"
#include "cells.h"

/*! \file cells.c
 *  \brief contains functions which control current states and evolution of cells
 */

unsigned char *celld;

/*!
 * This function checks whether the cell p is outside the computational box.
 */
static inline int outsideTheBox(int p)
{
  double x, y, z, r;

  x = cellsData.cells[p].x;
  y = cellsData.cells[p].y;
  z = cellsData.cells[p].z;
  r = cellsData.cells[p].size;

  if (x - r < 0 || x + r > (double) nx)
    return 1;
  if (y - r < 0 || y + r > (double) ny)
    return 1;
  if (sdim == 3 && (z - r < 0 || z + r > (double) nz))
    return 1;

  return 0;
}

/*!
 * This function allocates tables responsible for carrying
 * informations about cells, their current state and evolution.
 */
void cellsAllocate()
{

  int f;
  int64_t cellsActualSize;

  if (sdim == 2)
    csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 2.0);
  if (sdim == 3)
    csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 3.0);

  maxCellsPerProc = 1.5 * maxCells / MPIsize;

  cellsActualSize = maxCellsPerProc * sizeof(struct cellData);

  localID = 0;

#ifdef __MIC__
  if (!(cells = (struct cellData *) _mm_malloc(cellsActualSize,64)))
    stopRun(106, "cells", __FILE__, __LINE__);
#else
  if (!(cellsData.cells = (struct cellData *) malloc(cellsActualSize)))
    stopRun(106, "cells", __FILE__, __LINE__);
#endif

#ifdef __MIC__
  if (!(velocity = (struct doubleVector3d *) _mm_malloc(maxCellsPerProc *
                   sizeof(struct doubleVector3d),64)))
    stopRun(106, "velocity", __FILE__, __LINE__);
#else
  if (!(velocity = (struct doubleVector3d *) malloc(maxCellsPerProc *
                   sizeof(struct doubleVector3d))))
    stopRun(106, "velocity", __FILE__, __LINE__);
#endif

  cellsData.cellFields = (double **) malloc((NFIELDS+NCHEM*3) * sizeof(double *));
  for (f = 0; f < NFIELDS+NCHEM*3; f++)
    if (!
        (cellsData.cellFields[f] =
           (double *) malloc(maxCellsPerProc * sizeof(double))))
      stopRun(106, "cellFields", __FILE__, __LINE__);

  if (!(tlnc = (int64_t *) calloc(MPIsize, sizeof(int64_t))))
    stopRun(106, "tlnc", __FILE__, __LINE__);

  /* cell size */
  if (sdim == 2)
    csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 2.0);
  if (sdim == 3)
    csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 3.0);
  h=3.0*csize;
  simTime=0.0;
  h2 = h * h;
  h3 = h2 * h;
  h4 = h3 * h;


}

/*!
 * This function initializes counters of cells in various cell phases.
 */
void cellsCycleInit()
{
  /* global numbers of cells */
  g0nc = nc;
  g1nc = 0;
  snc = 0;
  g2nc = 0;
  mnc = 0;
  cnc = 0;
  /* local numbers of cells */
  lg0nc = lnc;
  lg1nc = 0;
  lsnc = 0;
  lg2nc = 0;
  lmnc = 0;
  lcnc = 0;
  lnnc = 0;
  /* number of cancer cells */
  cancer = 0;
  /* cell size */
  /*  if (sdim == 2)
      csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 2.0);
    if (sdim == 3)
      csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 3.0);
    h=3.0*csize;
    simTime=0.0;
    h2 = h * h;
    h3 = h2 * h;
    h4 = h3 * h;
  */
}

/*!
 * This function initializes cell data.
 * Locations of cells in space are generated randomly.
 */
int cellsRandomInit()
{

  int i;

  /* uniform distribution */
  if (!strcmp(rng, "UNB")) {
    double D=1.0;

    if (sdim == 2)
      csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 2.0);
    if (sdim == 3)
      csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 3.0);
    if (sdim == 2)
      D = csize * pow(8.0 * nc, 1.0 / 2.0);
    if (sdim == 3)
      D = csize * pow(8.0 * nc, 1.0 / 3.0);

    h = 3.0 * csize;
    simTime = 0;

    for (i = 0; i < lnc; i++) {
      cellsData.cells[i].x = D * (sprng(stream) * 2 - 1);
      cellsData.cells[i].y = D * (sprng(stream) * 2 - 1);
      if (sdim == 3)
        cellsData.cells[i].z = D * (sprng(stream) * 2 - 1);
      else
        cellsData.cells[i].z = 0.0;

      cellsData.cells[i].x += nx / 2;
      cellsData.cells[i].y += nx / 2;
      if (sdim == 3)
        cellsData.cells[i].z += nx / 2;
      else
        cellsData.cells[i].z = 0.0;

      cellsData.cells[i].size = pow(2.0, -(1.0 / 3.0)) * csize;
      cellsData.cells[i].gid =
        (unsigned long long int) MPIrank *(unsigned long long int)
        maxCellsPerProc + (unsigned long long int) i;
      cellsData.cells[i].v = 0.0;
      cellsData.cells[i].density = 0.0;
      cellsData.cells[i].h = h;
      cellsData.cells[i].young = (float) (2100.0 + sprng(stream) * 100.0);
      cellsData.cells[i].halo = 0;
      cellsData.cells[i].phase = 0;
      cellsData.cells[i].g1 = (float) (g1 * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].g2 = (float) (g2 * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].s = (float) (s * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].m = (float) (m * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].phasetime = 0.0;
      cellsData.cells[i].age = 0;
      cellsData.cells[i].death = 0;
      cellsData.cells[i].tumor = 0;
      cellsData.cells[i].ctype = 0;
      cellsData.cells[i].scstage = 0;
      localID++;
    }
    nscinst[0]+=lnc;
  }
  /* normal distribution (Box-Muller transform) */
  if (!strcmp(rng, "BM")) {
    double x1, x2, x3;
    double z1, z2, z3;
    double r1, r2;
    double l;
    double D=1.0;
    if (sdim == 2)
      csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 2.0);
    if (sdim == 3)
      csize = (nx / 2) / pow(8.0 * maxCells, 1.0 / 3.0);
    if (sdim == 2)
      D = csize * pow(8.0 * nc, 1.0 / 2.0);
    if (sdim == 3)
      D = csize * pow(8.0 * nc, 1.0 / 3.0);

    h = 3.0 * csize;
    simTime = 0;
    for (i = 0; i < lnc-lvc-lbnc; i++) {

      r2 = 1.1;

      while (r2 >= 1.0) {
        r1 = 1.1;
        while (r1 == 0 || r1 >= 1.0) {
          x1 = sprng(stream) * 2 - 1;
          x2 = sprng(stream) * 2 - 1;
          x3 = sprng(stream) * 2 - 1;
          r1 = x1 * x1 + x2 * x2 + x3 * x3;
        }
        l = sqrt(-2 * log(r1) / r1);
        z1 = x1 * l;
        z2 = x2 * l;
        z3 = x3 * l;

        r2 = z1 * z1 + z2 * z2 + z3 * z3;
      }

      if(bvsim) {
        cellsData.cells[i].x=cellsData.cells[middleCellIdx].x-0.8;
        cellsData.cells[i].y=cellsData.cells[middleCellIdx].y-0.8;
        cellsData.cells[i].z=cellsData.cells[middleCellIdx].z+0.01;
      } else {
        cellsData.cells[i].x = z1 * D + nx / 2;
        cellsData.cells[i].y = z2 * D + nx / 2;
        if (sdim == 3)
          cellsData.cells[i].z = z3 * D + nx / 2;
        else
          cellsData.cells[i].z = 0.0;
      }

      cellsData.cells[i].size = pow(2.0, -(1.0 / 3.0)) * csize;
      cellsData.cells[i].gid =
        (unsigned long long int) MPIrank *(unsigned long long int)
        maxCellsPerProc + (unsigned long long int) i;
      cellsData.cells[i].v = 0.0;
      cellsData.cells[i].density = 0.0;
      cellsData.cells[i].h = h;
      cellsData.cells[i].young = (float) (2100.0 + sprng(stream) * 100.0);
      cellsData.cells[i].halo = 0;
      cellsData.cells[i].phase = 0;
      cellsData.cells[i].g1 = (float) (g1 * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].g2 = (float) (g2 * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].s = (float) (s * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].m = (float) (m * (1 + (sprng(stream) * 2 - 1) * v));
      cellsData.cells[i].phasetime = 0.0;
      cellsData.cells[i].tumor = 0;
      cellsData.cells[i].age = 0;
      cellsData.cells[i].death = 0;
      cellsData.cells[i].ctype = 0;
      cellsData.cells[i].scstage = 0;
      localID++;
    }
  }

  nscinst[0]+=lnc;

  /* powers of h are calculated only once here */
  h2 = h * h;
  h3 = h2 * h;
  h4 = h3 * h;

  return 0;
}

/*!
 * This function implements mitosis of cells.
 */
void mitosis(int c)
{

  double sc;
  double shift[3];

  if (lnc + 1 > maxCellsPerProc)
    stopRun(109, NULL, __FILE__, __LINE__);

  /* stem cells counters */
  if(cellsData.cells[c].scstage==nscstages-1) {
    if(sprng(stream) > sctprob[cellsData.cells[c].scstage]) {
      celld[c]=1;
      localbc+=1;
    }
    return;
  }

  if(cellsData.cells[c].scstage<nscstages-1) {
    if(sprng(stream)>sctprob[cellsData.cells[c].scstage]) {
      cellsData.cells[lnc].scstage=cellsData.cells[c].scstage+1;
      nscinst[cellsData.cells[lnc].scstage]+=1;
    } else {
      cellsData.cells[lnc].scstage=cellsData.cells[c].scstage;
      nscinst[cellsData.cells[lnc].scstage]+=1;
    }
    if(sprng(stream)>sctprob[cellsData.cells[c].scstage]) {
      nscinst[cellsData.cells[c].scstage]-=1;
      cellsData.cells[c].scstage=cellsData.cells[c].scstage+1;
      nscinst[cellsData.cells[c].scstage]+=1;
    }
    /*if(sprng(stream)>sctprob[cells[c].scstage]) {
      cells[lnc].scstage=cells[c].scstage+1;
    } else {
      cells[lnc].scstage=cells[c].scstage;
    }*/
    //nscinst[cells[lnc].scstage]+=1;
  }

  sc = sqrt(velocity[c].x * velocity[c].x + velocity[c].y * velocity[c].y +
            velocity[c].z * velocity[c].z);

  /* daughter cells are shifted away from the center of parent cell */
  if (sc > 0 && mitrand == 0) {	/* direction of shift related to velocity vector */
    sc = cellsData.cells[c].size / (2 * sc);
    shift[0] = sc * velocity[c].x;
    shift[1] = sc * velocity[c].y;
    if (sdim == 3)
      shift[2] = sc * velocity[c].z;
    else
      shift[2] = 0.0;
  } else {			/* direction of shift chosen randomly */
    int accept = 0;
    while (accept == 0) {
      shift[0] = sprng(stream) * 2.0 - 1.0;
      shift[1] = sprng(stream) * 2.0 - 1.0;
      if (sdim == 3)
        shift[2] = sprng(stream) * 2.0 - 1.0;
      else
        shift[2] = 0.0;
      sc = sqrt(pow(shift[0], 2) + pow(shift[1], 2) + pow(shift[2], 2));
      if (sc == 0)
        continue;
      sc = cellsData.cells[c].size / (2 * sc);
      shift[0] = sc * shift[0];
      shift[1] = sc * shift[1];
      shift[2] = sc * shift[2];
      accept = 1;
    }
  }
  /* 1st daughter cell position, size, type and age */
  cellsData.cells[lnc].x = cellsData.cells[c].x + shift[0];
  cellsData.cells[lnc].y = cellsData.cells[c].y + shift[1];
  cellsData.cells[lnc].z = cellsData.cells[c].z + shift[2];
  cellsData.cells[lnc].size = pow(2.0, -(1.0 / 3.0)) * cellsData.cells[c].size;;
  cellsData.cells[lnc].tumor = cellsData.cells[c].tumor;
  cellsData.cells[lnc].age = cellsData.cells[c].age + 1;

  /* 2nd daughter cell position, size, type and age */
  cellsData.cells[c].x -= shift[0];
  cellsData.cells[c].y -= shift[1];
  cellsData.cells[c].z -= shift[2];
  cellsData.cells[c].size = cellsData.cells[lnc].size;;
  cellsData.cells[c].age += 1;

  /* 2nd daughter cell cycle phases lenghts */
  if (cellsData.cells[c].tumor == 1) {
    cellsData.cells[c].g1 = (float) (cg1 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[c].g2 = (float) (cg2 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[c].s = (float) (cs * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[c].m = (float) (cm * (1 + (sprng(stream) * 2 - 1) * v));
  } else {
    cellsData.cells[c].g1 = (float) (g1 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[c].g2 = (float) (g2 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[c].s = (float) (s * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[c].m = (float) (m * (1 + (sprng(stream) * 2 - 1) * v));
  }
  /* 1st daughter cell global ID */
  cellsData.cells[lnc].gid =
    (unsigned long long int) MPIrank *(unsigned long long int)
    maxCellsPerProc + (unsigned long long int) lnc;

  /* 1st daughter cell parameters */
  cellsData.cells[lnc].v = 0.0;
  cellsData.cells[lnc].density = cellsData.cells[c].density;
  cellsData.cells[lnc].h = h;
  cellsData.cells[lnc].young = (float) (2100.0 + sprng(stream) * 100.0);
  cellsData.cells[lnc].halo = 0;
  cellsData.cells[lnc].phase = 1;
  cellsData.cells[lnc].death = 0;
  cellsData.cells[lnc].phasetime = 0.0;
  /* 1st daughter cell cycle phases lenghts */
  if (cellsData.cells[lnc].tumor == 1) {
    cellsData.cells[lnc].g1 = (float) (cg1 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[lnc].g2 = (float) (cg2 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[lnc].s = (float) (cs * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[lnc].m = (float) (cm * (1 + (sprng(stream) * 2 - 1) * v));
  } else {
    cellsData.cells[lnc].g1 = (float) (g1 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[lnc].g2 = (float) (g2 * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[lnc].s = (float) (s * (1 + (sprng(stream) * 2 - 1) * v));
    cellsData.cells[lnc].m = (float) (m * (1 + (sprng(stream) * 2 - 1) * v));
  }

  /* update local cell counters */
  if (cellsData.cells[lnc].tumor == 1)
    lcnc += 1;
  lnc = lnc + 1;
  lg1nc += 1;
  /* increment local ID */
  localID++;

}

/*!
 * This function finds locates cell closest to the center of mass of the system
 * and marks this cell as a cancer cell.
 */
void markMiddleCancerCell()
{
  int c;
  int middle = 0;
  double dist;
  struct {
    double val;
    int rank;
  } lmdist, gmdist;
  double center[3];
  double gcenter[3];

  /* each process computes its local center of mass */
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;
  for (c = 0; c < lnc; c++) {
    center[0] += cellsData.cells[c].x / nc;
    center[1] += cellsData.cells[c].y / nc;
    center[2] += cellsData.cells[c].z / nc;
  }

  /* MPI Reduce operation computes global center of mass */
  MPI_Allreduce(center, gcenter, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* intialization */
  lmdist.rank = MPIrank;
  lmdist.val = INT_MAX;

  /* each process finds local cell closest to the global center of mass */
  for (c = 0; c < lnc; c++) {
    dist =
      sqrt((cellsData.cells[c].x - gcenter[0]) * (cellsData.cells[c].x - gcenter[0]) +
           (cellsData.cells[c].y - gcenter[1]) * (cellsData.cells[c].y - gcenter[1]) +
           (cellsData.cells[c].z - gcenter[2]) * (cellsData.cells[c].z - gcenter[2]));
    if (dist < lmdist.val) {
      lmdist.val = dist;
      middle = c;
    }
  }

  /* MPI_Allreduce locates the cell closest to the global center of mass */
  MPI_Allreduce(&lmdist, &gmdist, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                MPI_COMM_WORLD);
  /* mark the found cell as cancer one */
  if (MPIrank == gmdist.rank) {
    cellsData.cells[middle].tumor = 1;
    cellsData.cells[middle].phase = 1;
    lg0nc--;
    lg1nc++;
    lcnc++;
  }

  /* indicate that there is a cancer cell in the system */
  cancer = 1;
}

/*!
 * This function dealocates all tables allocated during initialization of cell data
 */
void cellsCleanup()
{
  int f;
  free(tlnc);
#ifdef __MIC__
  _mm_free(cells);
#else
  free(cellsData.cells);
#endif
  for (f = 0; f < NFIELDS; f++)
    free(cellsData.cellFields[f]);
  free(cellsData.cellFields);
#ifdef __MIC__
  _mm_free(velocity);
#else
  free(velocity);
#endif
}

/*!
 * This function removes a dead cell from the simulation.
 */
void cellsDeath(int lnc_old)
{
  int c, pos;

  pos = 0;
  for (c = 0; c < lnc; c++) {
    /* shift cells after dead cell removal */
    if (c >= lnc_old) {
      cellsData.cells[pos] = cellsData.cells[c];
      pos++;
      continue;
    }
    if (c != pos && celld[c] == 0)
      cellsData.cells[pos] = cellsData.cells[c];
    if (celld[c] == 0)
      pos++;
    /* update cell counters */
    if (celld[c] == 1) {
      switch (cellsData.cells[c].phase) {
      case 0:
        lg0nc--;
        break;
      case 1:
        lg1nc--;
        break;
      case 2:
        lsnc--;
        break;
      case 3:
        lg2nc--;
        break;
      case 4:
        lmnc--;
        break;
      }
      if (cellsData.cells[c].tumor == 1)
        lcnc--;
    }
  }
  lnc -= rsum;
}

/*!
 * This function updates cell counters.
 */
void updateCellCounters()
{
  MPI_Allgather(&lnc, 1, MPI_INT64_T, tlnc, 1, MPI_INT64_T,
                MPI_COMM_WORLD);
  MPI_Allreduce(localCellCount, totalCellCount, numberOfCounts,
                MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(nscinst, gnscinst, nscstages,
                MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localbc, &globalbc, 1,
                MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
}

void updateChemotaxis()
{
  int c;
  for(c=0; c<lnc; c++) {


  }
}


/*!
 * This function updates cells' positions.
 */
void updateCellPositions()
{
  int c;
#ifdef DEBUG
  if (MPIrank == 0 && !(step % statOutStep)) {
    printf(" Cells movement...");
    fflush(stdout);
  }
#endif
  if ((statistics.mindist >= 0.95 * 2.0 * pow(2.0, -(1.0 / 3.0)) * csize
       && simStart == 0) || (nc == 1 && simStart == 0)) {
    simStart = 1;
    if (MPIrank == 0)
      printf("\nSimulation started.\n");
  }

  /* move cells */
  for (c = 0; c < lnc; c++) {
    if(cellsData.cells[c].ctype==1) continue;
    cellsData.cells[c].x += velocity[c].x ;
    cellsData.cells[c].y += velocity[c].y ;
    cellsData.cells[c].z += velocity[c].z ;
    /* random movement */
    double alpha=0.000001;
    cellsData.cells[c].x += alpha*(2*sprng(stream)-1);
    cellsData.cells[c].y += alpha*(2*sprng(stream)-1);
    cellsData.cells[c].z += alpha*(2*sprng(stream)-1);

    // Mark cells that are out of the box and need to be removed
    //if(outside_the_box(c)) { celld[c]=1; rsum++; }
  }
#ifdef DEBUG
  if (MPIrank == 0 && !(step % statOutStep))
    printf("done\n");
#endif
}

/*!
 * This function updates cells' cycle phases.
 */
int updateCellCycles()
{

  int c;
  double eps, epsCancer;
  int lncAtThisStep;

  eps = densityCriticalLevel1;
  epsCancer = densityCriticalLevel2;

  lncAtThisStep = lnc;

  for (c = 0; c < lncAtThisStep; c++) {

    if (outsideTheBox(c)) {
      celld[c] = 1;
      rsum++;
      continue;
    }

    if (celld[c])
      continue;

    if (simStart) {

      if (cellsData.cells[c].phase != 0
          && ((cellsData.cells[c].tumor == 0 && cellsData.cells[c].density <= eps)
              || (cellsData.cells[c].tumor == 1 && cellsData.cells[c].density <= epsCancer)))
        cellsData.cells[c].phasetime += gfDt / 3600.0;

      switch (cellsData.cells[c].phase) {

      case 0:			/* G0 phase */
        if (gfields && oxygen
            && cellsData.cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
          cellsData.cells[c].phase = 5;
          cellsData.cells[c].phasetime = 0;
          lg0nc--;
          lnnc++;
          break;
        }
        /* transition to G1 phase */
        if ((cellsData.cells[c].tumor == 0 && cellsData.cells[c].density <= eps) ||	/* enough space for healthy cell */
            (cellsData.cells[c].tumor == 1 && cellsData.cells[c].density <= epsCancer) ||	/* enough space for tumor cell */
            nc == 1 ||		/* only single cell in the simulation */
            (gfields && oxygen && cellsData.cellFields[OXYG][c] >= fieldCriticalLevel1[OXYG])) {	/* sufficient level of oxygen */
          cellsData.cells[c].phase = 1;
          lg0nc--;
          lg1nc++;
          break;
        }
        break;
      case 1:			/* G1 phase */
        /* transition to G0 or Necrotic phase */
        if ((cellsData.cells[c].tumor == 0 && cellsData.cells[c].density > eps) ||	/* too crowdy for healthy cell */
            (cellsData.cells[c].tumor == 1 && cellsData.cells[c].density > epsCancer) ||	/* too crowdy for tumor cell */
            (gfields && oxygen && cellsData.cellFields[OXYG][c] < fieldCriticalLevel1[OXYG])) {	/* too low oxygen level */
          if (gfields && oxygen && cellsData.cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {	/* transition to Necrotic phase */
            cellsData.cells[c].phase = 5;
            cellsData.cells[c].phasetime = 0;
            lg1nc--;
            lnnc++;
          } else {		/* transition to G0 phase */
            cellsData.cells[c].phase = 0;
            lg1nc--;
            lg0nc++;
          }
          break;
        }
        /* cells grow in phase G1 */
        if (cellsData.cells[c].size < csize) {
          cellsData.cells[c].size +=
            (csize -
             pow(2.0,
                 -(1.0 / 3.0)) * csize) * (gfDt) / (3600.0 *
                    cellsData.cells[c].g1);
        }
        if (cellsData.cells[c].size > csize)
          cellsData.cells[c].size = csize;
        if (cellsData.cells[c].phasetime >= cellsData.cells[c].g1) {
          int death;
          cellsData.cells[c].phase = 2;
          cellsData.cells[c].phasetime = 0;
          lg1nc--;
          lsnc++;
          if (cellsData.cells[c].tumor == 0) {
            death = (sprng(stream) < rd ? 1 : 0);
            if (death) {
              celld[c] = 1;
              rsum++;
            }
          }
        }
        break;
      case 2:			/* S phase */
        if (gfields && oxygen
            && cellsData.cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
          cellsData.cells[c].phase = 5;
          cellsData.cells[c].phasetime = 0;
          lsnc--;
          lnnc++;
          break;
        }
        if (cellsData.cells[c].phasetime >= cellsData.cells[c].s) {
          cellsData.cells[c].phase = 3;
          cellsData.cells[c].phasetime = 0;
          lsnc--;
          lg2nc++;
          break;
        }
        break;
      case 3:			/* G2 phase */
        if (gfields && oxygen
            && cellsData.cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
          cellsData.cells[c].phase = 5;
          cellsData.cells[c].phasetime = 0;
          lg2nc--;
          lnnc++;
          break;
        }
        if (cellsData.cells[c].phasetime >= cellsData.cells[c].g2) {
          int death;
          cellsData.cells[c].phase = 4;
          cellsData.cells[c].phasetime = 0;
          lg2nc--;
          lmnc++;
          if (cellsData.cells[c].tumor == 0) {
            death = (sprng(stream) < rd ? 1 : 0);
            if (death) {
              celld[c] = 1;
              rsum++;
            }
          }
          break;
        }
        break;
      case 4:			/* M phase */
        if (gfields && oxygen
            && cellsData.cellFields[OXYG][c] < fieldCriticalLevel2[OXYG]) {
          cellsData.cells[c].phase = 5;
          cellsData.cells[c].phasetime = 0;
          lmnc--;
          lnnc++;

        } else if (cellsData.cells[c].phasetime >= cellsData.cells[c].m) {
          mitosis(c);
          cellsData.cells[c].phase = 1;
          cellsData.cells[c].phasetime = 0;
          lmnc--;
          lg1nc++;
        }
        break;
      }				// switch
    }				// if
  }				// for loop

  /* update global number of cells */
  MPI_Allreduce(&lnc, &nc, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

  return 0;
}

/*!
 * This function fills the scalarOutput field of each cell.
 * It can be modified to output any float scalars that
 * user would like to analyze or visualize after simulation.
 * This field is printed to the output VTK files.
 */
void additionalScalarField()
{
  int c;
  for (c = 0; c < lnc; c++) {
    if (cellsData.cells[c].tumor == 1)
      cellsData.cells[c].scalarField = 8.0;
    else
      cellsData.cells[c].scalarField = cellsData.cells[c].density;
  }
}

/*!
 * This function drives the whole cell cycle update.
 */
void updateCellStates()
{
  int lnc_old;
  /* number of local cells might change during the update */
  lnc_old = lnc;
  celld = (unsigned char *) calloc(lnc_old, sizeof(unsigned char));
  rsum = 0;

  updateCellCycles();
  if (nhs > 0 && nc > nhs && tgs == 1 && cancer == 0)
    markMiddleCancerCell();
  if (nhs > 0 && nc > nhs)
    cellsDeath(lnc_old);
  updateCellCounters();
  additionalScalarField();
  free(celld);
}

void constructCell(struct cellData * cell, int type_num, const struct cellTypeData * type){
  cell->g1 = (float) (type->g1 * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->g2 = (float) (type->g2 * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->s = (float) (type->s * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->m = (float) (type->m * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->ctype = type_num;
  cell->phasetime = 0;
  cell->phase = G1_phase;
  cell->young = (float) (2100.0 + sprng(stream) * 100.0);
}

void mitosis2(struct cellsInfo *ci, uint64_t cell_pos){
  struct cellData * cell_base;
  struct cellData * cell_new;
  cell_base = &ci->cells[cell_pos];
  if (ci->localCellCount.number_of_cells + 1 > maxCellsPerProc)
    stopRun(109, NULL, __FILE__, __LINE__);
#pragma omp critical (local_number_of_cells)
  {
    cell_new = &ci->cells[ci->localCellCount.number_of_cells];
    ci->localCellCount.number_of_cells++;
  }
  double new_size = pow(2.0, -(1.0 / 3.0)) * cell_base->size;
  constructCell(cell_new, cell_base->ctype, &ci->cellTypes[cell_base->ctype]);
  constructCell(cell_base, cell_base->ctype, &ci->cellTypes[cell_base->ctype]);
  cell_new->size = new_size;
  cell_base->size = new_size;
  cell_new->age = ++cell_base->age;
}

void updateCellCycles2(struct cellsInfo * ci, struct settings * s){
  uint64_t number_of_cells = ci->localCellCount.number_of_cells;
  struct cellData * cells = ci->cells;
  bool to_g0_phase;
  for(uint64_t i = 0; i < number_of_cells; i++){
    if (cells[i].phase == Necrotic_phase || cells[i].phase == const_cell)
      continue;
    bool alive = true;
    for (size_t j = 0; j < s->numberOfEnvironments; j++){
      if (ci->cellFields[i][j] < s->environments[j].critical_level_2){
        alive = false;
        break;
      }
    }
    if (!alive){
      switch (cells[i].phase){
        case G0_phase:
#pragma omp atomic
          ci->localCellCount.g0_phase_number_of_cells--;
          break;
        case G1_phase:
#pragma omp atomic
          ci->localCellCount.g1_phase_number_of_cells--;
          break;
        case S_phase:
#pragma omp atomic
          ci->localCellCount.s_phase_number_of_cells--;
          break;
        case G2_phase:
#pragma omp atomic
          ci->localCellCount.g2_phase_number_of_cells--;
          break;
        case M_phase:
#pragma omp atomic
          ci->localCellCount.m_phase_number_of_cells--;
          break;
        default: /*catch in first if in for loop */ ;
      }
#pragma omp atomic
      ci->localCellCount.necrotic_phase_number_of_cells++;
      cells[i].phase = Necrotic_phase;
    }
    switch (cells[i].phase){
      case G0_phase:
        if (cells[i].density <= ci->cellTypes[cells[i].ctype].eps ){
          cells[i].phase = G1_phase;
          //cells[i].phasetime = 0;
#pragma omp atomic
          ci->localCellCount.g0_phase_number_of_cells--;
#pragma omp atomic
          ci->localCellCount.g1_phase_number_of_cells++;
        };
        break;
      case G1_phase:
        to_g0_phase = false;
        if (cells[i].density > ci->cellTypes[cells[i].ctype].eps)
          to_g0_phase = true;
        for (size_t j = 0; j < s->numberOfEnvironments && (!to_g0_phase); j++){
          if (ci->cellFields[i][j] < s->environments[j].critical_level_1){
            to_g0_phase = true;
          }
        }
        if (to_g0_phase){
          cells[i].phase = G0_phase;
          //cells[i].phasetime = 0;
#pragma omp atomic
          ci->localCellCount.g0_phase_number_of_cells++;
#pragma omp atomic
          ci->localCellCount.g1_phase_number_of_cells--;
          break;
        }
        double max_size = ci->cellTypes[cells[i].ctype].max_size;
        if (cells[i].size < max_size){
          cells[i].size +=  (max_size - pow(2.0, -(1.0 / 3.0)) * max_size) *
                  (s->gloabal_fields_time_delta) / (3600.0 * cells[i].g1);
          //TODO This should depend of phase length and biomass income
        }
        if (cells[i].size > max_size){
          cells[i].size = max_size;
        }
        if (cells[i].phasetime >= cells[i].g1){
          cells[i].phase = S_phase;
          cells[i].phasetime = 0;
#pragma omp atomic
          ci->localCellCount.g1_phase_number_of_cells--;
#pragma omp atomic
          ci->localCellCount.s_phase_number_of_cells++;
          if (ci->cellTypes[cells[i].ctype].enable_random_death){
            if (sprng(stream) < rd){
              //FIXME I make assumption that in case of necrotic cells code should check that is enough place to destroy cell
              cells[i].phase = Necrotic_phase;
#pragma omp atomic
              ci->localCellCount.s_phase_number_of_cells--;
#pragma omp atomic
              ci->localCellCount.necrotic_phase_number_of_cells++;
            }

          }
        }
        break;
      case S_phase:
        if (cells[i].phasetime >= cells[i].s){
          cells[i].phase = G2_phase;
          cells[i].phasetime = 0;
#pragma omp atomic
          ci->localCellCount.s_phase_number_of_cells--;
#pragma omp atomic
          ci->localCellCount.g2_phase_number_of_cells++;
        }
        break;
      case G2_phase:
        if (cells[i].phasetime >= cells[i].g2){
          cells[i].phase = M_phase;
          cells[i].phasetime = 0;
#pragma omp atomic
          ci->localCellCount.g2_phase_number_of_cells--;
#pragma omp atomic
          ci->localCellCount.m_phase_number_of_cells++;
          if (ci->cellTypes[cells[i].ctype].enable_random_death){
            if (sprng(stream) < rd){
              //FIXME I make assumption that in case of necrotic cells code should check that is enough place to destroy cell
              cells[i].phase = Necrotic_phase;
#pragma omp atomic
              ci->localCellCount.m_phase_number_of_cells--;
#pragma omp atomic
              ci->localCellCount.necrotic_phase_number_of_cells++;
            }

          }
        }break;
      case M_phase:
        if (cells[i].phasetime >= cells[i].m){

        }
        break;
      default: /*catch in first if in for loop */ ;
    }
  }

}

void initiateCellsInfo(struct cellsInfo *ci, const struct settings *s){
  ci->localTypeCellCount = (uint64_t *)
          calloc(s->numberOfCellTypes, sizeof(uint64_t));
  if (!ci->localTypeCellCount){
    stopRun(106, "cellsInfo::localTypeCellCount", __FILE__, __LINE__);
  }
  ci->totalTypeCellCount = (uint64_t *)
          calloc(s->numberOfCellTypes, sizeof(uint64_t));
  if (!ci->totalTypeCellCount){
    stopRun(106, "cellsInfo::totalTypeCellCount", __FILE__, __LINE__);
  }
  // Memory pointed by ci->localCellCount and ci->totalCellCount is set to 0 by calloc
  ci->cells = (struct cellData *) aligned_alloc(64, s->maxCellsPerProc * sizeof(struct cellData) );
  if (!ci->cells){
    stopRun(106, "cellsInfo::cells", __FILE__, __LINE__);
  }
  size_t cellFieldsSize = s->numberOfEnvironments * s->maxCellsPerProc * sizeof(double);
  double * fields = (double *) aligned_alloc(64, cellFieldsSize);
  if (!fields){
    stopRun(106, "initiateCellsInfo::fields", __FILE__, __LINE__);
  }
  ci->cellFields = (double **) aligned_alloc(64, s->maxCellsPerProc * sizeof(double *));
  if (!ci->cellFields){
    stopRun(106, "cellsInfo::cellFieldsSize", __FILE__, __LINE__);
  }
  for (size_t i = 0; i < s->maxCellsPerProc; i++){
    ci->cellFields[i] = &fields[s->numberOfEnvironments * i];
  }
  memset(ci->cellFields, 0, cellFieldsSize);
  ci->velocity = (struct doubleVector3d *) aligned_alloc(64, sizeof(struct doubleVector3d) * s->maxCellsPerProc);
  if (!ci->cells){
    stopRun(106, "cellsInfo::velocity", __FILE__, __LINE__);
  }
  memset(ci->velocity, 0, sizeof(struct doubleVector3d) * s->maxCellsPerProc);
  ci->cellTypes = s->cellTypes;
  ci->cellTypeNumberDict = s->cellTypeNumberDict;
}

void freeCellsInfo(struct cellsInfo *ci){
  free(ci->localTypeCellCount);
  free(ci->totalTypeCellCount);
  free(ci->cells);
  free(ci->cellFields[0]);
  free(ci->cellFields);
  free(ci->velocity);
}

void changeCellType(const struct cellsInfo *ci, struct cellData * cell, int ctype){
  uint64_t * counter;
  counter = &ci->localTypeCellCount[cell->ctype];
#pragma omp atomic
  (*counter)--;
  cell->ctype = ctype;
  counter = &ci->localTypeCellCount[cell->ctype];
#pragma omp atomic
  (*counter)++;
}

void changeCellTypeByName(const struct cellsInfo *ci, struct cellData * cell, char * name){
  changeCellType(ci, cell, get_value(ci->cellTypeNumberDict, name));
}



