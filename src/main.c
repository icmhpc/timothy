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

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include <stdbool.h>

#include "global.h"
#include "io.h"
#include "init.h"
#include "domdec.h"
#include "comm.h"
#include "octree.h"
#include "compute.h"
#include "stats.h"
#include "cells.h"
#include "random.h"
/*! \file main.c
 *  \brief contains the main simulation loop
 */

/*!
 * This function intializes MPI, calls Timothy initialization and allocation functions.
 * It also contains the main simulation loop where all important simulation steps are called.
*/

struct state State;
struct settings mainSettings;
struct cellsInfo cellsData;
struct statisticsData MIC_ATTR statistics;
octNode MIC_ATTR *octTree;


double **cellFields;
float simTime;          /* time of the simulation */

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &State.MPIsize);
  State.uMPIsize = (unsigned int) State.MPIsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &State.MPIrank);

  State.OMPthreads = omp_get_max_threads();

  simulationInit(argc, argv);

  for (State.step = 0; State.step < mainSettings.number_of_steps; State.step++) {

    if (!(State.step % mainSettings.statOutStep))
      printStepNum();

    decompositionExecute(cellsData.totalCellCount.number_of_cells);
    octBuild(&cellsData);
    createExportList(&cellsData);
    computeStep(&cellsData);
    if (!(State.step % mainSettings.statOutStep))
      statisticsPrint(&cellsData, &statistics);

    if (State.simStart)
      simTime += mainSettings.secondsPerStep / 3600.0;	/* biological process time in hours */


    if (!(State.step % mainSettings.vtkOutStep)) {
      if (mainSettings.vtkout)
        ioWriteStepVTK(State.step);
      if (mainSettings.povout)
        ioWriteStepPovRay(State.step, 0);
      if (mainSettings.vnfout)
        ioWriteFields(State.step);
    }

    updateCellPositions(&cellsData, &statistics);
    updateCellStates(&cellsData, &mainSettings);
    commCleanup(cellsData.totalCellCount.number_of_cells);
    octFree(&cellsData);

    if (!(State.step % mainSettings.rstOutStep))
      saveRstFile();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  decompositionFinalize();
  randomStreamFree();
  freeCellsInfo(&cellsData);

  if (State.MPIrank == 0)
    printf("\nEnd of simulation run.\n");

  MPI_Finalize();

  return 0;
}
