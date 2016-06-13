//
// Created by czaki on 10.06.16.
//

#include <stdlib.h>
#include "interpolation.h"

/*!
 * for each cell calculate process which handle envirorment information for this cell
 * @param ci
 * @param st
 * @param grid
 * @return array: cell index -> process number
 */
int * findCellsProcess(struct cellsInfo * ci, struct state *st, struct gridData * grid){
  int * res;
  size_t number_of_cells = ci->localCellCount.number_of_cells;
  struct cellData * cells = ci->cells;
  res = calloc(number_of_cells, sizeof(int));
  doubleVector3d position;
  struct intVector3d cords;
  for (size_t i = 0; i < number_of_cells; ++i){
    position.x = cells[i].x - grid->lower_corner.x;
    position.y = cells[i].y - grid->lower_corner.y;
    position.z = cells[i].z - grid->lower_corner.z;
    cords.x = (int) (position.x * st->MPIdim[0]/grid->global_size.x);
    cords.y = (int) (position.y * st->MPIdim[1]/grid->global_size.y);
    cords.z = (int) (position.z * st->MPIdim[2]/grid->global_size.z);
    int procNum = *getProcesNum(st, cords);
    if (position.x < grid->lower_indices[procNum].x){
      cords.x--;
      procNum = *getProcesNum(st,cords);
    }
    if (position.y < grid->lower_indices[procNum].y){
      cords.y--;
      procNum = *getProcesNum(st,cords);
    }
    if (position.z < grid->lower_indices[procNum].z){
      cords.z--;
      procNum = *getProcesNum(st,cords);
    }
    if (position.x > grid->upper_indices[procNum].x){
      cords.x++;
      procNum = *getProcesNum(st,cords);
    }
    if (position.y > grid->upper_indices[procNum].y){
      cords.y++;
      procNum = *getProcesNum(st,cords);
    }
    if (position.z > grid->upper_indices[procNum].z){
      cords.z++;
      procNum = *getProcesNum(st,cords);
    }

    res[i] = procNum;
  }
  return res;
}