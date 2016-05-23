//
// Created by czaki on 14.04.16.
//

#ifndef TIMOTHY_FUNCTIONS_H
#define TIMOTHY_FUNCTIONS_H

#include "global.h"

typedef void (*cellFunction)(cellData * cell, float * env, void * data);

typedef double (*functionalEnv)(doubleVector3d, void * data);

cellFunction getCellFunction(std::string s);
functionalEnv getFunEnvFunction(std::string s);


#endif //TIMOTHY_FUNCTIONS_H
