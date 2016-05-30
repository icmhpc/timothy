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

/*! \file fields.h
 *  \brief contains variables and arrays for global fields
 */
#ifndef TIMOTHY_FIELDS_H
#define TIMOTHY_FIELDS_H
#include <mpi.h>
#include <stdint.h>
#define NFIELDS 6 //5
#define NIF 2 /* number of fields that need to be interpolated from discrete data */
#define NGLOB 3
#define NCHEM 3
#define SCALAR_FIELD 1
#define VECTOR_FIELD 2

int nfields;
int nif;
int nglob;
int nchem;

/* grid data */
int64_t gridI,gridJ,gridK;
struct int64Vector3d gridSize;
struct int64Vector3d *gridStartIdx,*gridEndIdx;
struct doubleVector3d *gridBuffer;
#define grid(i,j,k) (gridBuffer[gridSize.y*gridSize.z*i+gridSize.z*j+k])
double gridResolution;
struct doubleVector3d lowerGridCorner,upperGridCorner;
struct doubleVector3d globalGridSize;
double boxVolume;

/* field data */
char fieldName[NFIELDS][128];
double *fieldAddr[NFIELDS];
double *gradAddr[NCHEM];
int fieldType[NFIELDS];

double fieldDt[NFIELDS]; /* how many seconds per iteration */

int fieldNumberOfCriticalLevels[NFIELDS]; /* number of critical levels for each field */

//double *fieldCriticalLevel1; /* critical level 1 - for cells in G1 phase */
//double *fieldCriticalLevel2; /* critical level 2 - for cells in S, G2, M phases */

double fieldMin[NFIELDS];
double fieldMax[NFIELDS];
//double *fieldMin;
//double *fieldMax;

double *tissueField;
double *vesselField;
double *tempField;
double *chemField[NCHEM];
double *chemGradient[NCHEM];
//double **chemField;
//double **chemGradient;
double **fieldsPatchesCommBuff;
double **fieldsPatches;


MPI_Request reqFGSend[6];
MPI_Request reqFGRecv[6];

#define DENS 0
#define BVES 1
#define TEMP 2
#define OXYG 3
#define GLUC 4
#define HYDR 5



/* GLOBAL SETTINGS */
int gfields;
int oxygen;
int glucose;
int hydrogenIon;
int temperature;
double fieldDiffCoef[NFIELDS];
double fieldLambda[NFIELDS];
double fieldBC[NFIELDS]; /* boundary conditions - constant value for Dirichlet */
double fieldICMean[NFIELDS]; /* initial conditions - mean value */
double fieldICVar[NFIELDS]; /* initial conditions - variability */
double fieldConsumption[NFIELDS]; /* units - mol (cell s)^-1 */
double fieldProduction[NFIELDS];  /* units - mol (cell s)^-1 */
double fieldCriticalLevel1[NFIELDS]; /* critical level 1 - for cells in G1 phase */
double fieldCriticalLevel2[NFIELDS]; /* critical level 2 - for cells in S, G2, M phases */


void fieldsInit();
void fieldsSolve();
void fieldGradient();

#endif /* ifndef TIMOTHY_FIELDS_H */
