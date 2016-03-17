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
#include <string.h>

#include "global.h"
#include "domdec.h"
#include "utils.h"

/*! \file domdec.c
 *  \brief contains domain decomposition functions
 */

/* arrays and variables used by Zoltan library */
int changes;			/* 1 if partitioning was changed, 0 otherwise */
int numGidEntries;		/* number of integers used for a global ID */
int numLidEntries;		/* number of integers used for a local ID */
int numImport;			/* number of objects to be sent to me */
ZOLTAN_ID_PTR importGlobalGids;	/* global IDs of objects to be sent to me */
ZOLTAN_ID_PTR importLocalGids;	/* local IDs of objects to be sent to me */
int *importProcs;		/* process rank for source of each incoming object */
int *importToPart;		/* new partition for each incoming object */
int numExport;			/* number of objects I must send to other processes */
ZOLTAN_ID_PTR exportGlobalGids;	/* global IDs of the objects I must send */
ZOLTAN_ID_PTR exportLocalGids;	/* local IDs of the objects I must send */
int *exportProcs;		/* process to which I send each of the objects */
int *exportToPart;		/* partition to which each object will belong */

struct Zoltan_Struct *ztn;


/*!
 * Zoltan callback function. This function returns the dimension (2D or 3D) of the system.
 */
int ztnReturnDimension(void *data __attribute__((unused)), int *ierr __attribute__((unused)))
{
  if (sdim == 3)
    return 3;
  if (sdim == 2)
    return 2;
  return 0;
}

/*!
 * Zoltan callback function. This function returns the spatial coordinates of the cell identified by its global and local id.
 */
void ztnReturnCoords(void *data, int numGidEntries __attribute__((unused)), int numLidEntries __attribute__ ((unused)),
                     ZOLTAN_ID_PTR globalId __attribute__((unused)), ZOLTAN_ID_PTR localId,
                     double *geomVec, int *ierr __attribute__((unused)))
{
  struct cellData * cells = (struct cellData *) data;
  if (sdim == 3) {
    geomVec[0] = cells[localId[0]].x;
    geomVec[1] = cells[localId[0]].y;
    geomVec[2] = cells[localId[0]].z;
  }
  if (sdim == 2) {
    geomVec[0] = cells[localId[0]].x;
    geomVec[1] = cells[localId[0]].y;
  }
}

/*!
 * Zoltan callback function. This function returns the number of cells assigned to this process.
 */
int ztnReturnNumNode(void *data __attribute__((unused)), int *ierr __attribute__((unused)))
{
  uint64_t * localCellCount = (uint64_t *) data;
  return (int) *localCellCount;
}

/*!
 * Zoltan callback function. This function fills the tables of global ids, local ids and weights for all cells assigned to this process.
 */
void ztnReturnOwnedNodes(void *data, int numGIdEntries, int numLIdEntries,
                         ZOLTAN_ID_PTR globalIds, ZOLTAN_ID_PTR localIds,
                         int wgtDim __attribute__((unused)), float *objWgts, int *ierr __attribute__((unused)))
{
  struct cellData * cells = ((struct cellsInfo *) data)->cells;
  uint64_t i;
  uint64_t local_number_of_cells = ((struct cellsInfo *) data)->localCellCount.number_of_cells;
  for (i = 0; i < local_number_of_cells; i++) { //TODO maybe use wgtDim instead of lnc?
    globalIds[i * numGIdEntries] = cells[i].gid;
    localIds[i * numLIdEntries] = i;
    objWgts[i] = 1.0;
    //if(nc==1 || step==0) objWgts[i]=1.0;
    //else objWgts[i]=cells[i].density;
  }
}

/*!
 * Zoltan callback function. This function returns the size of a data structure used for keeping the data of a single cell.
 */
int ztnReturnParticleDataSize(void *data __attribute__((unused)), int numGIdEntries __attribute__((unused)),
                              int numLIdEntries __attribute__((unused)), ZOLTAN_ID_PTR globalId __attribute__((unused)),
                              ZOLTAN_ID_PTR localId __attribute__((unused)), int *ierr __attribute__((unused)))
{
  return sizeof(struct cellData);
}

/*!
 * Zoltan callback function. This function packs data into a send buffer before migration.
 */
void ztnPack(void *data, int numGIdEntries __attribute__((unused)), int numLIdEntries __attribute__((unused)),
             ZOLTAN_ID_PTR globalId __attribute__((unused)), ZOLTAN_ID_PTR localId, int dest __attribute__((unused)),
             int size __attribute__((unused)), char *buf, int *ierr __attribute__((unused)))
{
  struct cellData *c = (struct cellData *) data;
  memcpy(buf, &(c[localId[0]]), sizeof(struct cellData));
  c[(int) (*localId)].gid = -1;	/* mark local particle as exported */
}

/*!
 * Zoltan callback function. This function is executed before migration of data between processes.
 */
void ztnPre(void *data __attribute__((unused)), int numGIdEntries __attribute__((unused)), int numLIdEntries __attribute__((unused)),
            int numImport __attribute__((unused)), ZOLTAN_ID_PTR importGlobalIds __attribute__((unused)),
            ZOLTAN_ID_PTR importLocalIds __attribute__((unused)), int *importProcs __attribute__((unused)),
            int *importToPart __attribute__((unused)), int numExport __attribute__((unused)),
            ZOLTAN_ID_PTR exportGlobalIds __attribute__((unused)), ZOLTAN_ID_PTR exportLocalIds __attribute__((unused)),
            int *exportProcs __attribute__((unused)), int *exportToPart __attribute__((unused)), int *ierr __attribute__((unused)))
{
  /* any pre communication operations should go here */
}

/*!
 * Zoltan callback function. This function is executed after packing of send buffer and unpacking of receive buffer during migration.
 */
void ztnMid(void *data, int numGIdEntries __attribute__((unused)), int numLIdEntries __attribute__((unused)),
            int numImport __attribute__((unused)), ZOLTAN_ID_PTR importGlobalIds __attribute__((unused)),
            ZOLTAN_ID_PTR importLocalIds __attribute__((unused)), int *importProcs __attribute__((unused)),
            int *importToPart __attribute__((unused)), int numExport,
            ZOLTAN_ID_PTR exportGlobalIds __attribute__((unused)), ZOLTAN_ID_PTR exportLocalIds __attribute__((unused)),
            int *exportProcs __attribute__((unused)), int *exportToPart __attribute__((unused)), int *ierr __attribute__((unused)))
{
  uint64_t pos, i;
  //struct ztnMidData * d = (struct ztnMidData * ) data;
  struct cellData *c = ((struct cellsInfo *) data)->cells;
  uint64_t * localCellCount = &((struct cellsInfo *) data)->localCellCount.number_of_cells;
  pos = 0;
  for (i = 0; i < *localCellCount; i++) {
    if (i != pos && c[i].gid != (unsigned long) -1) { //TODO why compare unsigned long with -1?
      c[pos] = c[i];
    }
    if (c[i].gid != (unsigned long) -1)
      pos++;
  }
  *localCellCount = *localCellCount - numExport;
}

/*!
 * Zoltan callback function. This function is executed after migration of data between processes.
 */
void ztnPost(void *data, int numGIdEntries __attribute__((unused)), int numLIdEntries __attribute__((unused)),
             int numImport __attribute__((unused)), ZOLTAN_ID_PTR importGlobalIds __attribute__((unused)),
             ZOLTAN_ID_PTR importLocalIds __attribute__((unused)), int *importProcs __attribute__((unused)),
             int *importToPart __attribute__((unused)), int numExport __attribute__((unused)),
             ZOLTAN_ID_PTR exportGlobalIds __attribute__((unused)), ZOLTAN_ID_PTR exportLocalIds __attribute__((unused)),
             int *exportProcs __attribute__((unused)), int *exportToPart __attribute__((unused)), int *ierr __attribute__((unused)))
{
  struct cellsInfo * ci = (struct cellsInfo *) data;
  /* any post communication operations should go here */
  /* gather number of cells from each process */
  MPI_Allgather(&ci->localCellCount.number_of_cells, 1, MPI_UINT64_T, ci->numberOfCellsInEachProcess, 1, MPI_UINT64_T,
                MPI_COMM_WORLD);
}

/*!
 * Zoltan callback function. This function unpacks data from the receive buffer.
 */
void ztnUnpack(void *data, int numGIdEntries __attribute__((unused)), ZOLTAN_ID_PTR globalId __attribute__((unused)),
               int size __attribute__((unused)), char *buf, int *ierr __attribute__((unused)))
{
  struct cellsInfo * ci = (struct cellsInfo * ) data;
  //struct cellData *c = d->c;
  //int64_t * localCellCount = d->localCellCount;
  //struct cellData *c = (struct cellData *) data;
  memcpy(&ci->cells[ci->localCellCount.number_of_cells], buf, sizeof(struct cellData));
  ci->localCellCount.number_of_cells++;
}

/*!
 * This function initializes the Zoltan library.
 * It is called at the beginning of the simulation.
 */
void decompositionInit(int argc, char **argv)
{
  int rc;
  float version;

  rc = Zoltan_Initialize(argc, argv, &version);
  if (rc != ZOLTAN_OK)
    stopRun(112, NULL, __FILE__, __LINE__);

  if (State.MPIrank == 0)
    printf("Zoltan Version %.3f. Initialized.\n", version);

  ztn = Zoltan_Create(MPI_COMM_WORLD);



  Zoltan_Set_Param(ztn, "IMBALANCE_TOL", "1.4");
  Zoltan_Set_Param(ztn, "LB_METHOD", "HSFC");	/* Hilbert Space-Filling Curve Partitioning */
  Zoltan_Set_Param(ztn, "NUM_GID_ENTRIES", "1");	/* global ID is 1 integer */
  Zoltan_Set_Param(ztn, "NUM_LID_ENTRIES", "1");	/* local ID is 1 integer */
  Zoltan_Set_Param(ztn, "OBJ_WEIGHT_DIM", "1");	/* we use object weights */
  Zoltan_Set_Param(ztn, "DEBUG_LEVEL", "0");	/* quiet mode; no output unless an error or warning is produced */
  Zoltan_Set_Param(ztn, "KEEP_CUTS", "1");	/* save the cuts for later use */
  Zoltan_Set_Param(ztn, "AUTO_MIGRATE", "1");	/* use the auto migration mechanism */

  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_GEOM_FN_TYPE,
                (void (*)()) ztnReturnDimension, cellsData.cells);
  Zoltan_Set_Fn(ztn, ZOLTAN_GEOM_FN_TYPE, (void (*)()) ztnReturnCoords,
                cellsData.cells);
  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) ztnReturnNumNode,
                &cellsData.localCellCount.number_of_cells);
  Zoltan_Set_Fn(ztn, ZOLTAN_OBJ_LIST_FN_TYPE,
                (void (*)()) ztnReturnOwnedNodes, &cellsData);
  Zoltan_Set_Fn(ztn, ZOLTAN_OBJ_SIZE_FN_TYPE,
                (void (*)()) ztnReturnParticleDataSize, cellsData.cells);
  Zoltan_Set_Fn(ztn, ZOLTAN_PACK_OBJ_FN_TYPE, (void (*)()) ztnPack, cellsData.cells);
  Zoltan_Set_Fn(ztn, ZOLTAN_UNPACK_OBJ_FN_TYPE, (void (*)()) ztnUnpack,
                &cellsData);
  Zoltan_Set_Fn(ztn, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) ztnPre,
                cellsData.cells);
  Zoltan_Set_Fn(ztn, ZOLTAN_MID_MIGRATE_PP_FN_TYPE, (void (*)()) ztnMid,
                &cellsData);
  Zoltan_Set_Fn(ztn, ZOLTAN_POST_MIGRATE_PP_FN_TYPE, (void (*)()) ztnPost,
                &cellsData);

}

/*!
 * This function calls the Zoltan's domain decomposition and migration functions.
 * It is called at the beginning of each simulation step.
 */
void decompositionExecute(uint64_t total_number_of_cells)
{
  int rc;

  if (total_number_of_cells < State.uMPIsize*MIN_CELLS_PER_PROC)
    return;

  rc = Zoltan_LB_Partition(ztn,	/* input (all remaining fields are output) */
                           &changes,	/* 1 if partitioning was changed, 0 otherwise */
                           &numGidEntries,	/* number of integers used for a global ID */
                           &numLidEntries,	/* number of integers used for a local ID */
                           &numImport,	/* number of objects to be sent to me */
                           &importGlobalGids,	/* global IDs of objects to be sent to me */
                           &importLocalGids,	/* local IDs of objects to be sent to me */
                           &importProcs,	/* process rank for source of each incoming object */
                           &importToPart,	/* new partition for each incoming object */
                           &numExport,	/* number of objects I must send to other processes */
                           &exportGlobalGids,	/* global IDs of the objects I must send */
                           &exportLocalGids,	/* local IDs of the objects I must send */
                           &exportProcs,	/* process to which I send each of the objects */
                           &exportToPart);	/* partition to which each object will belong */

  if (rc != ZOLTAN_OK)
    stopRun(112, NULL, __FILE__, __LINE__);

  /* free the arrays allocated by Zoltan_LB_Partiotion */
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs,
                      &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs,
                      &exportToPart);
}

/*!
 * This function deactivates the Zoltan library.
 * It is called at the end of the simulation.
 */
void decompositionFinalize()
{
  Zoltan_Destroy(&ztn);
}
