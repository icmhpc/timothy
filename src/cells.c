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
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include "dicts.h"
#include "global.h"
#include "fields.h"
#include "utils.h"
#include "cells.h"

/*! \file cells.c
 *  \brief contains functions which control current states and evolution of cells
 */
extern int *stream;

/*!
 * This function checks whether the cell p is outside the computational box.
 */
static inline bool outsideTheBox(struct cellData *cell, const struct settings *s)
{
  double x, y, z, r;

  x = cell->x;
  y = cell->y;
  z = cell->z;
  r = cell->size;

  if (x - r < 0 || x + r > (double) s->size_x)
    return true;
  if (y - r < 0 || y + r > (double) s->size_y)
    return true;
  if (s->numOfDimensions == 3 && (z - r < 0 || z + r > (double) s->size_z))
    return true;

  return false;
}




/* move to files with function which can be choose by settings */

/*!
 * This function finds locates cell (of type from type) closest to the center of mass of the system
 * and marks this cell as a given to_type.
 */
void markMiddleCell(struct cellsInfo *ci, int to_type, int from_type)
{
  uint64_t c;
  int middle = 0;
  double dist;
  struct {
    double val;
    int rank;
  } local_min_dist, global_min_dist;
  double center[3];
  double global_center[3];

  uint64_t local_number_off_cells, total_number_off_cells;
  local_number_off_cells = ci->localCellCount.number_of_cells;
  total_number_off_cells = ci->totalCellCount.number_of_cells;

  /* each process computes its local center of mass */
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;
  for (c = 0; c < local_number_off_cells; c++) {
    center[0] += ci->cells[c].x / total_number_off_cells;
    center[1] += ci->cells[c].y / total_number_off_cells;
    center[2] += ci->cells[c].z / total_number_off_cells;
  }

  /* MPI Reduce operation computes global center of mass */
  MPI_Allreduce(center, global_center, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* initialization */
  local_min_dist.rank = State.MPIrank;
  local_min_dist.val = INT_MAX;

  /* each process finds local cell closest to the global center of mass */
  for (c = 0; c < local_number_off_cells; c++) {
    /* check that cell can be changed */
    if (ci->cells[c].cell_type == from_type && ci->cells[c].phase != Necrotic_phase &&
            ci->cells[c].phase != const_cell && ci->cells[c].phase != cell_to_destroy) {
      dist = (ci->cells[c].x - global_center[0]) * (ci->cells[c].x - global_center[0]) +
             (ci->cells[c].y - global_center[1]) * (ci->cells[c].y - global_center[1]) +
             (ci->cells[c].z - global_center[2]) * (ci->cells[c].z - global_center[2]);
      if (dist < local_min_dist.val) {
        local_min_dist.val = dist;
        middle = (int) c;
      }
    }
  }

  /* MPI_Allreduce locates the cell closest to the global center of mass */
  MPI_Allreduce(&local_min_dist, &global_min_dist, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                MPI_COMM_WORLD);
  /* mark the found cell as cancer one */
  if (State.MPIrank == global_min_dist.rank) {
    struct cellData * cell  = &ci->cells[middle];
    switch (cell->phase){
      case G0_phase: ci->localCellCount.g0_phase_number_of_cells--; break;
      case G1_phase: ci->localCellCount.g1_phase_number_of_cells--; break;
      case S_phase: ci->localCellCount.s_phase_number_of_cells--; break;
      case G2_phase: ci->localCellCount.g2_phase_number_of_cells--; break;
      case M_phase: ci->localCellCount.m_phase_number_of_cells--; break;
      case Necrotic_phase: ci->localCellCount.necrotic_phase_number_of_cells--;break;
      case const_cell: ci->localCellCount.g0_phase_number_of_cells--; break;
      case cell_to_destroy:break;
    }
    ci->localCellCount.g1_phase_number_of_cells++;
    ci->localTypeCellCount[from_type]--;
    ci->localTypeCellCount[to_type]++;
    ci->cells[middle].cell_type = to_type;
    ci->cells[middle].phase = G1_phase;
  }
}





/*!
 * This function updates cells' positions.
 */
void updateCellPositions(struct cellsInfo *ci, struct statisticsData *st)
{
  uint64_t c;
#ifdef DEBUG
  if (MPIrank == 0 && !(step % statOutStep)) {
    printf(" Cells movement...");
    fflush(stdout);
  }
#endif
  if ((st->mindist >= 0.95 * 2.0 * pow(2.0, -(1.0 / 3.0)) * csize
       && State.simStart == 0) || (ci->totalCellCount.number_of_cells == 1 && State.simStart == 0)) {
    State.simStart = 1;
    if (State.MPIrank == 0)
      printf("\nSimulation started.\n");
  }

  /* move cells */
  uint64_t local_number_of_cells = ci->localCellCount.number_of_cells;
  for (c = 0; c < local_number_of_cells; c++) {
    if(ci->cells[c].cell_type==1) continue;
    ci->cells[c].x += ci->velocity[c].x ;
    ci->cells[c].y += ci->velocity[c].y ;
    ci->cells[c].z += ci->velocity[c].z ;
    /* random movement */
    double alpha=0.000001; //FIXME setings?
    ci->cells[c].x += alpha*(2*sprng(stream)-1);
    ci->cells[c].y += alpha*(2*sprng(stream)-1);
    ci->cells[c].z += alpha*(2*sprng(stream)-1);

    // Mark cells that are out of the box and need to be removed
    //if(outside_the_box(c)) { celld[c]=1; rsum++; }
  }
#ifdef DEBUG
  if (MPIrank == 0 && !(step % statOutStep))
    printf("done\n");
#endif
}


/*!
 * This function fills the scalarOutput field of each cell.
 * It can be modified to output any float scalars that
 * user would like to analyze or visualize after simulation.
 * This field is printed to the output VTK files.
 */
/*void additionalScalarField()
{
  int c;
  for (c = 0; c < lnc; c++) {
    if (cellsData.cells[c].tumor == 1)
      cellsData.cells[c].scalarField = 8.0;
    else
      cellsData.cells[c].scalarField = cellsData.cells[c].density;
  }
}*/




void updateCellCounters(struct cellsInfo *ci) {
  MPI_Allgather(&ci->localCellCount.number_of_cells, 1, MPI_UINT64_T,
                ci->numberOfCellsInEachProcess, 1, MPI_UINT64_T,
                MPI_COMM_WORLD);
  MPI_Allreduce(&ci->localCellCount, &ci->totalCellCount, sizeof(struct cellCountInfo)/sizeof(uint64_t),
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ci->localTypeCellCount, &ci->totalTypeCellCount, (int) ci->numberOfCellTypes,
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
}


/*!
 * This function removes a dead cell from the simulation.
 */
void cellsDeath(struct cellsInfo *ci){
  struct cellData * it;
  struct cellData * reverse_it;
  double minimal_density = ci->minimal_density_to_keep_necrotic_cell;
  it = ci->cells;
  reverse_it = &ci->cells[ci->localCellCount.number_of_cells - 1];
  while (it < reverse_it){
    while (reverse_it->phase == cell_to_destroy ||
            (reverse_it->phase == Necrotic_phase && reverse_it->density < minimal_density)){
      if (reverse_it->phase == Necrotic_phase) {
        ci->localCellCount.necrotic_phase_number_of_cells--;
      }
      ci->localTypeCellCount[reverse_it->cell_type]--;
      ci->localCellCount.number_of_cells--;
      reverse_it--;
    }

    if (it->phase == cell_to_destroy ||
              (it< reverse_it && Necrotic_phase == it->cell_type && it->density < minimal_density)){
      if (it->phase == Necrotic_phase) {
        ci->localCellCount.necrotic_phase_number_of_cells--;
      }
      ci->localTypeCellCount[it->cell_type]--;
      ci->localCellCount.number_of_cells--;
      *it = *reverse_it;
    }
    it++;
  }
}


/*!
 * This function drives the whole cell cycle update.
 */
void updateCellStates(struct cellsInfo *ci, const struct settings *s){
  updateCellCycles(ci, s);
  if (s->enable_step_transformation && s->step_transformation != NULL){
    (*s->step_transformation)(ci, s);
  }
  cellsDeath(ci);
  updateCellCounters(ci);
  /* no cancer mark and additional scalr fields*/

}

void constructCell(struct cellData * cell, int type_num, const struct cellTypeData * type){
  cell->g1 = (float) (type->g1 * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->g2 = (float) (type->g2 * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->s = (float) (type->s * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->m = (float) (type->m * (1 + (sprng(stream) * 2 - 1) * type->v));
  cell->cell_type = type_num;
  cell->phasetime = 0;
  cell->phase = G1_phase;
  cell->young = (float) (2100.0 + sprng(stream) * 100.0);
}

/*!
 * This function implements mitosis of cells.
 */
void mitosis(struct cellsInfo *ci, uint64_t cell_pos, const struct settings *s){
  struct cellData * cell_base;
  struct cellData * cell_new;
  struct doubleVector3d shift;
  double velocity_norm, mother_cell_radius;
  if (ci->localCellCount.number_of_cells <= cell_pos){
    stopRun(666, "Wrong position of mitosis cell", __FILE__, __LINE__);
  }
  velocity_norm = sqrt(ci->velocity[cell_pos].x * ci->velocity[cell_pos].x + ci->velocity[cell_pos].y * ci->velocity[cell_pos].y +
                    ci->velocity[cell_pos].z * ci->velocity[cell_pos].z);
  cell_base = &ci->cells[cell_pos];
  cell_base->age++;
  mother_cell_radius = cell_base->size;
  if (velocity_norm > 0 && !s->mitosis_random_direction){
    shift.x = ci->velocity[cell_pos].x * (mother_cell_radius/(2 * velocity_norm));
    shift.y = ci->velocity[cell_pos].y * (mother_cell_radius/(2 * velocity_norm));
    shift.z = ci->velocity[cell_pos].z * (mother_cell_radius/(2 * velocity_norm));
  } else {
    velocity_norm = 0.0;
    while (velocity_norm == 0.0) {
      shift.x = sprng(stream) * 2.0 - 1.0;
      shift.y = sprng(stream) * 2.0 - 1.0;
      if (sdim == 3)
        shift.z = sprng(stream) * 2.0 - 1.0;
      else
        shift.z = 0.0;
      velocity_norm = sqrt(pow(shift.x, 2) + pow(shift.y, 2) + pow(shift.z, 2));
    }
    shift.x = shift.x * (mother_cell_radius / (2 * velocity_norm));
    shift.y = shift.y * (mother_cell_radius / (2 * velocity_norm));
    shift.z = shift.z * (mother_cell_radius / (2 * velocity_norm));
  }
  if (ci->localCellCount.number_of_cells + 1 > s->maxCellsPerProc)
    stopRun(109, NULL, __FILE__, __LINE__);
#pragma omp critical (local_number_of_cells)
  {
    cell_new = &ci->cells[ci->localCellCount.number_of_cells];
    ci->localCellCount.number_of_cells++;
  }
  *cell_new = *cell_base;
  double new_size = pow(2.0, -(1.0 / s->dimension)) * cell_base->size;
#pragma omp atomic
  ci->localCellCount.m_phase_number_of_cells--;
#pragma omp atomic
  ci->localCellCount.g1_phase_number_of_cells += 2;
#pragma omp atomic
  ci->localCellCount.number_of_cells++;
#pragma omp atomic
  ci->localTypeCellCount[cell_new->cell_type]++;
  constructCell(cell_new, cell_base->cell_type, &ci->cellTypes[cell_base->cell_type]);
  constructCell(cell_base, cell_base->cell_type, &ci->cellTypes[cell_base->cell_type]);
  cell_new->size = new_size;
  cell_new->x += shift.x;
  cell_new->y += shift.y;
  cell_new->z += shift.z;

  cell_base->size = new_size;
  cell_base->x -= shift.x;
  cell_base->y -= shift.y;
  cell_base->z -= shift.z;
  cell_new->gid = State.MPIrank * s->id_range + State.localID;
  State.localID++;
}

/*!
 * This function updates cells' cycle phases.
 */
void updateCellCycles(struct cellsInfo *ci, const struct settings *s){
  uint64_t number_of_cells = ci->localCellCount.number_of_cells;
  struct cellData * cells = ci->cells;
  bool to_g0_phase;
  for(uint64_t i = 0; i < number_of_cells; i++){
    if (cells[i].phase == Necrotic_phase || cells[i].phase == const_cell)
      continue;

    if (outsideTheBox(&cells[i], s)) {
      switch (cells[i].phase){
        case G0_phase:
#pragma omp atomic
          ci->localCellCount.g0_phase_number_of_cells--;
          break;
        case G1_phase:
#pragma omp atomic
          ci->localCellCount.g1_phase_number_of_cells--;break;
        case S_phase:
#pragma omp atomic
          ci->localCellCount.s_phase_number_of_cells--;break;
        case G2_phase:
#pragma omp atomic
          ci->localCellCount.g2_phase_number_of_cells--;break;
        case M_phase:
#pragma omp atomic
          ci->localCellCount.m_phase_number_of_cells--;break;
        case Necrotic_phase:
#pragma omp atomic
          ci->localCellCount.necrotic_phase_number_of_cells--;break;
        case const_cell:break;
        case cell_to_destroy:break;
      }
      cells[i].phase = cell_to_destroy;
      continue;
    }


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
    if (cells[i].phase != G0_phase)
      cells[i].phasetime += s->global_fields_time_delta / 3600.0;
    switch (cells[i].phase){
      case G0_phase:
        if (cells[i].density <= ci->cellTypes[cells[i].cell_type].eps ){
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
        if (cells[i].density > ci->cellTypes[cells[i].cell_type].eps)
          to_g0_phase = true;
        for (size_t j = 0; j < s->numberOfEnvironments && (!to_g0_phase); j++){
          if (ci->cellFields[i][j] < s->environments[j].critical_level_1){
            to_g0_phase = true;
          }
        }
        if (to_g0_phase){
          cells[i].phase = G0_phase;
          cells[i].phasetime -= s->global_fields_time_delta / 3600.0;
          //cells[i].phasetime = 0;
#pragma omp atomic
          ci->localCellCount.g0_phase_number_of_cells++;
#pragma omp atomic
          ci->localCellCount.g1_phase_number_of_cells--;
          break;
        }{
          double max_size = ci->cellTypes[cells[i].cell_type].max_size;
          if (cells[i].size < max_size){
            cells[i].size +=  (max_size - pow(2.0, -(1.0 / 3.0)) * max_size) *
                    (s->global_fields_time_delta) / (3600.0 * cells[i].g1);
            //FIXME This should depend of phase length and biomass income
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
            if (ci->cellTypes[cells[i].cell_type].enable_random_death){
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
          if (ci->cellTypes[cells[i].cell_type].enable_random_death){
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
          mitosis(ci, i, s);
        }
        break;
      default: /*catch in first if in for loop */ ;
    }
  }
  MPI_Allreduce(&ci->localCellCount.number_of_cells, &ci->totalCellCount.number_of_cells, 1,
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
}

/*!
 * This function allocates tables responsible for carrying
 * informations about cells, their current state and evolution.
 */
void initiateCellsInfo(struct cellsInfo *ci, const struct settings *s){
  // s->maxCellsPerProc = 1.5 * s->maxCells; TODO Add in create setings
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
  memset(ci->cellFields[0], 0, cellFieldsSize);
  ci->velocity = (struct doubleVector3d *) aligned_alloc(64, sizeof(struct doubleVector3d) * s->maxCellsPerProc);
  if (!ci->cells){
    stopRun(106, "cellsInfo::velocity", __FILE__, __LINE__);
  }
  memset(ci->velocity, 0, sizeof(struct doubleVector3d) * s->maxCellsPerProc);
  ci->numberOfCellsInEachProcess = (uint64_t *) aligned_alloc(64, sizeof(int64_t) * s->MPI_group_size);
  ci->cellTypes = s->cellTypes;
  ci->cellTypeNumberDict = s->cellTypeNumberDict;
  // set counters to 0
  memset(&ci->localCellCount, 0, sizeof(struct cellCountInfo));
  memset(&ci->totalCellCount, 0, sizeof(struct cellCountInfo));
}


/*!
 * This function dealocates all tables allocated during initialization of cell data
 */
void freeCellsInfo(struct cellsInfo *ci){
  free(ci->localTypeCellCount);
  free(ci->totalTypeCellCount);
  free(ci->cells);
  free(ci->cellFields[0]);
  free(ci->cellFields);
  free(ci->velocity);
}

void changeCellType(struct cellsInfo *ci, struct cellData *cell, int cell_type){
  uint64_t * counter;
  counter = &ci->localTypeCellCount[cell->cell_type];
#pragma omp atomic
  (*counter)--;
  cell->cell_type = cell_type;
  counter = &ci->localTypeCellCount[cell->cell_type];
#pragma omp atomic
  (*counter)++;
}

void changeCellTypeByName(struct cellsInfo *ci, struct cellData *cell, char *name){
  changeCellType(ci, cell, get_value(ci->cellTypeNumberDict, name));
}

void middleCellInit(struct cellsInfo *ci, const struct settings *s, int type){
  struct cellData * cell;
  cell = &ci->cells[ci->localCellCount.number_of_cells];
  constructCell(cell, type, &s->cellTypes[type]);
  cell->x = s->size_x/2;
  cell->y = s->size_y/2;
  cell->z = s->size_z/2;
  cell->age = 0;
  cell->v = 0.0;
  cell-> density = 0.0;
  cell->h = s->neighbourhood;
  cell->halo = 0;

  ci->localCellCount.g1_phase_number_of_cells++;
  ci->localTypeCellCount[type]++;
}

