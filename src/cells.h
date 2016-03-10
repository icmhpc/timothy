#ifndef TIMOTHY_CELLS_H
#define TIMOTHY_CELLS_H
void cellsAllocate();
void updateCellPositions(struct cellsInfo *ci, struct statisticsData *st);
void updateCellStates();
void cellsCleanup();
void cellsCycleInit();
int cellsRandomInit();
void initiateCellsInfo(struct cellsInfo *ci, const struct settings *s);
void freeCellsInfo(struct cellsInfo *ci);
void changeCellType(struct cellsInfo *ci, struct cellData *cell, int cell_type);
void changeCellTypeByName(struct cellsInfo *ci, struct cellData *cell, char *name);
void updateCellCycles2(struct cellsInfo *ci, const struct settings *s);
#endif //TIMOTHY_CELLS_H


