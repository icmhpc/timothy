//
// Created by czaki on 03.06.16.
//

#ifndef TIMOTHY_FIELDS_H
#define TIMOTHY_FIELDS_H
#include "environment.h"

void allocateFieldGradient(const struct settings *set, const struct state *sta,
                           const struct gridData *grid);

void initFieldHaloExchange(const struct state *sta, const struct settings *set,
                           const struct gridData *grid);

void freeFieldGradient();

void fieldsInit(const struct settings *set, const struct state *sta,
                const struct gridData *grid)
#endif // TIMOTHY_FIELDS_H
