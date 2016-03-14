#ifndef TIMOTHY_COMM_H
#define TIMOTHY_COMM_H
void createExportList(struct cellsInfo *ci);
void commCleanup(uint64_t total_number_of_cells);
void cellsExchangeInit(struct cellsInfo *ci);
void cellsExchangeWait(struct cellsInfo *ci);
void densPotExchangeInit(struct cellsInfo *ci);
void densPotExchangeWait(struct cellsInfo *ci);
#endif /* ifndef TIMOTHY_COMM_H */
