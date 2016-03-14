#ifndef TIMOTHY_CELLS_H
#define TIMOTHY_CELLS_H

void updateCellPositions(struct cellsInfo *ci, struct statisticsData *st);
void initiateCellsInfo(struct cellsInfo *ci, const struct settings *s);
void freeCellsInfo(struct cellsInfo *ci);
void changeCellType(struct cellsInfo *ci, struct cellData *cell, int cell_type);
void changeCellTypeByName(struct cellsInfo *ci, struct cellData *cell, char *name);
void updateCellCycles(struct cellsInfo *ci, const struct settings *s);
void updateCellStates(struct cellsInfo *ci, const struct settings *s);
void middleCellInit(struct cellsInfo *ci, const struct settings *s, int type);
void updateCellCounters(struct cellsInfo *ci);
#endif //TIMOTHY_CELLS_H


