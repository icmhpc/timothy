//
// Created by czaki on 10.06.16.
//

#include <stdlib.h>
#include "interpolation.h"


int * findCellsProcess(struct cellsInfo * ci, struct state *st, struct gridData * grid){
  int * res;
  size_t number_of_cells = ci->localCellCount.number_of_cells;
  struct cellData * cells = ci->cells;
  res = calloc(number_of_cells, sizeof(int));
  doubleVector3d position;
  for (size_t i = 0; i < number_of_cells; ++i){
    position.x = cells[i].x - grid->lower_corner.x;
    position.y = cells[i].y - grid->lower_corner.y;
    position.z = cells[i].z - grid->lower_corner.z;
  }
  return res;
}