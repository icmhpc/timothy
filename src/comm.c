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

#include <stdlib.h>

#include "global.h"
#include "utils.h"

/*! \file comm.c
 *  \brief contains communication functions
 */

/* data type for IDs of exported cells */
struct expData {
  int cell;
  int proc;
};

MPI_Request *reqSend;
MPI_Request *reqRecv;

int64_t *sendOffset;
int64_t *recvOffset;

struct expData *expList;

int *recvCount;
int *sendCount;

#define MAX_EXPORTED_PER_PROC 2*maxCellsPerProc

/*!
 * This function is a comparison function used
 * for sorting export list table.
 */
int comm_compare_exp_list(const void *a, const void *b)
{
  return ((struct expData *) a)->proc - ((struct expData *) b)->proc;
}

/*!
 * This function uses Zoltan's library function Zoltan_LB_Box_Assign
 * to find possible intersections of cells' neighbourhoods
 * and other processes' geometries.
 */
void createExportList()
{

  int i, p;
  int procs[MPIsize];
  int numprocs;

  numExp = 0;
  numImp = 0;
  if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
    return;

  expList =
    (struct expData *) malloc(sizeof(struct expData) *
                              MAX_EXPORTED_PER_PROC);
  recvCount = (int *) calloc(MPIsize, sizeof(int));
  sendCount = (int *) calloc(MPIsize, sizeof(int));
  sendOffset = (int64_t *) calloc(MPIsize, sizeof(int64_t));
  recvOffset = (int64_t *) calloc(MPIsize, sizeof(int64_t));

  /* loop over local cells */
  /*#pragma omp parallel for private(procs) */
  for (p = 0; p < lnc; p++) {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double r;

    cellsData.cells[p].halo = 0;

    if (nc < MPIsize*MIN_CELLS_PER_PROC)
      continue;

    r = h * 1.5;

    /* compute neighbourhood box */
    xmin = cellsData.cells[p].x - r;
    xmax = cellsData.cells[p].x + r;
    ymin = cellsData.cells[p].y - r;
    ymax = cellsData.cells[p].y + r;
    if (sdim == 3) {
      zmin = cellsData.cells[p].z - r;
      zmax = cellsData.cells[p].z + r;
    } else {
      zmin = 0.0;
      zmax = 0.0;
    }

    /* look for possible neighbours */
    Zoltan_LB_Box_Assign(ztn, xmin, ymin, zmin, xmax, ymax, zmax, procs,
                         &numprocs);

    /* loop over receivers */
    for (i = 0; i < numprocs; i++) {
      if (procs[i] == MPIrank || tlnc[procs[i]] == 0)
        continue;
      expList[numExp].cell = p;
      expList[numExp].proc = procs[i];
      cellsData.cells[p].halo = MPIrank + 1;
      sendCount[procs[i]]++;
      numExp++;
      /* upps! too many refugees */
      if (numExp >= MAX_EXPORTED_PER_PROC)
        stopRun(110, NULL, __FILE__, __LINE__);
    }
  }

  /* sort export list with respect to process number */
  qsort(expList, numExp, sizeof(struct expData), comm_compare_exp_list);

  /* distribute the information on transfer sizes between each process */
  MPI_Alltoall(sendCount, 1, MPI_INT, recvCount, 1, MPI_INT,
               MPI_COMM_WORLD);

  /* compute send offsets */
  sendOffset[0] = 0;
  for (i = 1; i < MPIsize; i++)
    sendOffset[i] = sendOffset[i - 1] + sendCount[i - 1];

  /* compute receive offsets */
  recvOffset[0] = 0;
  for (i = 1; i < MPIsize; i++)
    recvOffset[i] = recvOffset[i - 1] + recvCount[i - 1];

  /* count cells to be imported */
  for (i = 0; i < MPIsize; i++)
    numImp += recvCount[i];

}

/*!
 * This function deallocates all communication buffers and
 * auxiliary tables.
 */
void commCleanup()
{
  if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
    return;

#ifdef __MIC__
  _mm_free(recvData);
  _mm_free(recvDensPotData);
#else
  free(recvData);
  free(recvDensPotData);
#endif

  free(expList);

  free(recvCount);
  free(sendCount);

  free(sendOffset);
  free(recvOffset);

}

/*!
 * This function initiate sending and receiving cells' data between processes.
 */
void cellsExchangeInit()
{
  int i;

  if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
    return;

  /* allocate communication buffers */
  sendData = (struct partData *) malloc(numExp * sizeof(struct partData));
#ifdef __MIC__
  recvData = (struct partData *) _mm_malloc(numImp * sizeof(struct partData),64);
#else
  recvData = (struct partData *) malloc(numImp * sizeof(struct partData));
#endif
  reqSend = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);
  reqRecv = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);

  /* create reduced particle data buffer for exporting */
  for (i = 0; i < numExp; i++) {
    sendData[i].x = cellsData.cells[expList[i].cell].x;
    sendData[i].y = cellsData.cells[expList[i].cell].y;
    sendData[i].z = cellsData.cells[expList[i].cell].z;
    sendData[i].size = cellsData.cells[expList[i].cell].size;
    sendData[i].h = h;
    sendData[i].young = (double) cellsData.cells[expList[i].cell].young;
    sendData[i].ctype = cellsData.cells[expList[i].cell].ctype;
  }

  /* send cells - asynchronous MPI call */
  for (i = 0; i < MPIsize; i++) {
    if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    MPI_Isend(&sendData[sendOffset[i]],
              sendCount[i] * sizeof(struct partData), MPI_BYTE, i, MPIrank,
              MPI_COMM_WORLD, &reqSend[i]);
  }

  /* receive cells - asynchronous MPI call */
  for (i = 0; i < MPIsize; i++) {
    if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    MPI_Irecv(&recvData[recvOffset[i]],
              recvCount[i] * sizeof(struct partData), MPI_BYTE, i, i,
              MPI_COMM_WORLD, &reqRecv[i]);
  }

}

/*!
 * This function waits for cells' data exchange completion.
 */
void cellsExchangeWait()
{
  int i;
  MPI_Status status;

  if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
    return;

  /* wait for send completion */
  for (i = 0; i < MPIsize; i++) {
    if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    if (MPI_Wait(&reqSend[i], &status) != MPI_SUCCESS)
      stopRun(103, "reqSend", __FILE__, __LINE__);
  }

  /* wait for receive completion */
  for (i = 0; i < MPIsize; i++) {
    if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    if (MPI_Wait(&reqRecv[i], &status) != MPI_SUCCESS)
      stopRun(103, "reqRecv", __FILE__, __LINE__);
  }

  /* some of the buffers can be deallocated here */
  free(sendData);
  free(reqSend);
  free(reqRecv);
}

/*!
 * This function initiate sending and receiving density
 * and potential values between processes.
 */
void densPotExchangeInit()
{
  int i;

  if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
    return;

  /* allocate communication buffers */
  sendDensPotData =
    (struct densPotData *) malloc(numExp * sizeof(struct densPotData));
#ifdef __MIC__
  recvDensPotData =
    (struct densPotData *) _mm_malloc(numImp * sizeof(struct densPotData),64);
#else
  recvDensPotData =
    (struct densPotData *) malloc(numImp * sizeof(struct densPotData));
#endif
  reqSend = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);
  reqRecv = (MPI_Request *) malloc(sizeof(MPI_Request) * MPIsize);

  /* create density and potential buffer for exporting */
  for (i = 0; i < numExp; i++) {
    sendDensPotData[i].v = cellsData.cells[expList[i].cell].v;
    sendDensPotData[i].density = cellsData.cells[expList[i].cell].density;
  }

  /* send data - asynchronous MPI call */
  for (i = 0; i < MPIsize; i++) {
    if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    MPI_Isend(&sendDensPotData[sendOffset[i]],
              sendCount[i] * sizeof(struct densPotData), MPI_BYTE, i,
              MPIrank, MPI_COMM_WORLD, &reqSend[i]);
  }

  /* receive data - asynchronous MPI call */
  for (i = 0; i < MPIsize; i++) {
    if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    MPI_Irecv(&recvDensPotData[recvOffset[i]],
              recvCount[i] * sizeof(struct densPotData), MPI_BYTE, i, i,
              MPI_COMM_WORLD, &reqRecv[i]);
  }

}

/*!
 * This function waits for density and potential data exchange completion.
 */
void densPotExchangeWait()
{
  int i;
  MPI_Status status;

  if (nc < MPIsize*MIN_CELLS_PER_PROC || MPIsize == 1)
    return;

  // Wait for send completion
  for (i = 0; i < MPIsize; i++) {
    if (sendCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    if (MPI_Wait(&reqSend[i], &status) != MPI_SUCCESS)
      stopRun(103, "sending", __FILE__, __LINE__);
  }

  // Wait for receive completion
  for (i = 0; i < MPIsize; i++) {
    if (recvCount[i] == 0 || tlnc[i] == 0 || lnc == 0)
      continue;
    if (MPI_Wait(&reqRecv[i], &status) != MPI_SUCCESS)
      stopRun(103, "receiving", __FILE__, __LINE__);
  }

  /* some of the buffers can be deallocated */
  free(sendDensPotData);
  free(reqSend);
  free(reqRecv);
}
