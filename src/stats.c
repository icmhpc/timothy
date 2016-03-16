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

#include "global.h"

/*! \file stats.c
 *  \brief contains functions computing and printing simulation statisctical information
 */

/*!
 * This function computes and prints out density statistics.
 */
void statisticsDensity(struct cellsInfo *ci, struct statisticsData *st)
{
  uint64_t c;
  double sum = 0.0, mean;
  double globalMaxDens;
  double globalMinDens;

  st->maxdens = 0;
  st->mindens = 1024;

  for (c = 0; c < ci->localCellCount.number_of_cells; c++) {
    sum += ci->cells[c].density;
    st->maxdens =
      (ci->cells[c].density >
       st->maxdens ? ci->cells[c].density : st->maxdens);
    st->mindens =
      (ci->cells[c].density <
       st->mindens ? ci->cells[c].density : st->mindens);
  }

  MPI_Allreduce(&sum, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&st->maxdens, &globalMaxDens, 1, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&st->mindens, &globalMinDens, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  st->maxdens = globalMaxDens;
  st->mindens = globalMinDens;

  mean /= ci->totalCellCount.number_of_cells;
  st->densavg = mean;
  sum = 0.0;

  for (c = 0; c < ci->localCellCount.number_of_cells; c++)
    sum += (ci->cells[c].density - mean) * (ci->cells[c].density - mean);

  MPI_Allreduce(&sum, &(st->densdev), 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  st->densdev /= ci->totalCellCount.number_of_cells;
  st->densdev = sqrt(st->densdev);

  if (State.MPIrank == 0)
    printf("%12s%10.4lf%10.4lf%10.4lf%10.4lf\n", "Density    ",
           globalMinDens, globalMaxDens, mean, st->densdev);
}

/*!
 * This function computes and prints out distance statistics.
 */
void statisticsDistance(struct statisticsData *st)
{
  double globalMinDist;

  MPI_Allreduce(&st->mindist, &globalMinDist, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);

  st->mindist = globalMinDist;

  if (State.MPIrank == 0)
    printf("%12s%10.4lf%10s%10s%10s\n", "Distance   ", globalMinDist, "-",
           "-", "-");
}

/*!
 * This function computes and prints out velocity statistics.
 */
void statisticsVelocity(struct statisticsData *st)
{
  MPI_Reduce(&st->minvel, &globalMinVel, 1, MPI_DOUBLE, MPI_MIN, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&st->maxvel, &globalMaxVel, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  if (State.MPIrank == 0)
    printf("%12s%10.4lf%10.4lf%10s%10s\n", "Velocity   ", globalMinVel,
           globalMaxVel, "-", "-");
}

/*!
 * This function computes and prints out cell size statistics.
 */
void statisticsSize(struct statisticsData *st)
{
  double globalMaxSize, globalMinSize;

  MPI_Reduce(&st->maxsize, &globalMaxSize, 1, MPI_DOUBLE, MPI_MAX,
             0, MPI_COMM_WORLD);
  MPI_Reduce(&st->minsize, &globalMinSize, 1, MPI_DOUBLE, MPI_MIN,
             0, MPI_COMM_WORLD);

  if (State.MPIrank == 0)
    printf("%12s%10.4lf%10.4lf%10s%10s\n", "Size       ", globalMinSize,
           globalMaxSize, "-", "-");
}

/*!
 * This function computes and prints out phases statistics.
 */
void statisticsPhases(struct cellsInfo *ci)
{
  if (State.MPIrank == 0) {
    printf("%12s%12s%12s%12s%12s%12s\n", "Cell phase ", "G0", "G1", "S",
           "G2", "M");
    printf("%12s%12" PRId64 "%12" PRId64 "%12" PRId64 "%12" PRId64 "%12"
           PRId64 "\n", "N. of cells", ci->totalCellCount.g0_phase_number_of_cells,
           ci->totalCellCount.g1_phase_number_of_cells, ci->totalCellCount.s_phase_number_of_cells,
           ci->totalCellCount.g2_phase_number_of_cells, ci->totalCellCount.m_phase_number_of_cells);
    for (size_t i=0; i < ci->numberOfCellTypes; i++){
      printf("%16s%12" PRId64 "\n", ci->cellTypes[i].name, ci->totalTypeCellCount[i]);
    }
    printf("%16s%12" PRId64 "\n", "Necrotic cells ", ci->totalCellCount.necrotic_phase_number_of_cells);

  }
}

/*!
 * This is a driving function for computing and printing out statistics.
 */
void statisticsPrint(struct cellsInfo *ci, struct statisticsData *st)
{
  if (State.MPIrank == 0)
    printf("%12s%10s%10s%10s%10s\n", "", "Min", "Max", "Avg", "Dev");
  statisticsDensity(ci, st);
  statisticsDistance(st);
  /*statisticsVelocity(); */
  statisticsSize(st);
  if (State.MPIrank == 0)
    printf("\n");
  statisticsPhases(ci);
}
