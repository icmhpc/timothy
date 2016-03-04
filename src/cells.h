#ifndef TIMOTHY_CELLS_H
#define TIMOTHY_CELLS_H
void cellsAllocate();
void updateCellPositions();
void updateCellStates();
void cellsCleanup();
void cellsCycleInit();
int cellsRandomInit();
void initiateCellsInfo(struct cellsInfo *ci, const struct settings *s);
void freeCellsInfo(struct cellsInfo *ci);
void changeCellType(const struct cellsInfo *ci, struct cellData * cell, int ctype);
void changeCellTypeByName(const struct cellsInfo *ci, struct cellData * cell, char * name);
#endif //TIMOTHY_CELLS_H


