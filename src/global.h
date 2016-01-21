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
#ifndef TIMOTHY_GLOBAL_H
#define  TIMOTHY_GLOBAL_H

#include <zoltan.h>
#include<mpi.h>
#include<stdint.h>

/*! \file global.h
 *  \brief contains the most important global variables, arrays and defines
 */

#define  VERSION   "1.0"

#ifdef __MIC__
#define MIC_ATTR __attribute__((target(mic)))
#else
#define MIC_ATTR
#endif

#define POWER_OF_TWO(x) !(x&(x-1))

/* architecture */
#if defined(__ia64) || defined(__itanium__) || defined(_M_IA64)
#define CPUARCH "Itanium"
#endif
#if defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
#if defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__) || \
	defined(__64BIT__) || defined(_LP64) || defined(__LP64__)
#define CPUARCH "POWER64"
#else
#define CPUARCH "POWER32"
#endif
#endif
#if defined(__sparc)
#define CPUARCH "SPARC"
#endif
#if defined(__x86_64__) || defined(_M_X64)
#define CPUARCH "x86_64"
#elif defined(__i386) || defined(_M_IX86)
#define CPUARCH "x86_32"
#endif

/* cell data structure */
//#pragma pack(1)
struct cellData {
  int lifetime;         /* age of the cell */
  int phase;            /* actual phase of the cell (0=G0,1=G1,2=S,3=G2,4=M,5=Necrotic) */
  int age;              /* cell's age */
  int death;		/* 1 - dead cell, 0 - living cell */
  int halo;		/* cell on the border of parallel region */
  float phasetime;      /* actual phase time */
  float g1;		/* g1 phase lenght - randomly selected for each new cell */
  float s;		/* s phase lenght - randomly selected for each new cell */
  float g2;		/* g2 phase lenght - randomly selected for each new cell */
  float m;		/* m phase lenght - randomly selected for each new cell */
  float young;		/* Young modulus - randomly selected for each new cell */
  ZOLTAN_ID_TYPE gid;   /* global ID of the particle */
  double x;             /* x coordinate of the particle position */
  double y;             /* y coordinate of the particle position */
  double z;             /* z coordinate of the particle position */
  double size;          /* radius of the cell */
  double h;             /* neighbourhood of the cell */
  double v;             /* particle potential */
  double density;       /* particle density */
  double scalarField;   /* additional scalar field which can be used for analysis of results (printed to the output VTK files) */
  int ctype;		/* cell type 1=endothelial */
  int scstage;          /* stem cells stage (in the course of differentiation) */
  unsigned char tumor;	/* 1 - tumor cell, 0 - healthly cell */
};

/* !!!!!!!!!!!!!!!!!!!!!!! */
/* the most important data */
MIC_ATTR struct cellData *cells;		/* main array for keeping cell data */
double **cellFields;				/* fields value for each cell - interpolated from global fields in each step */
MIC_ATTR struct doubleVector3d *velocity;	/* velocity table - velocity of each cell modified in each step */
/* !!!!!!!!!!!!!!!!!!!!!!! */

int64_t maxCells;
#define numberOfCounts 10

MIC_ATTR int64_t localCellCount[numberOfCounts];
int64_t totalCellCount[numberOfCounts];

#define nc   totalCellCount[0]
#define g0nc totalCellCount[1]
#define g1nc totalCellCount[2]
#define snc  totalCellCount[3]
#define g2nc totalCellCount[4]
#define mnc  totalCellCount[5]
#define cnc  totalCellCount[6]
#define nnc  totalCellCount[7]
#define vc   totalCellCount[8]
#define bnc  totalCellCount[9]

#define lnc   localCellCount[0]
#define lg0nc localCellCount[1]
#define lg1nc localCellCount[2]
#define lsnc  localCellCount[3]
#define lg2nc localCellCount[4]
#define lmnc  localCellCount[5]
#define lcnc  localCellCount[6]
#define lnnc  localCellCount[7]
#define lvc   localCellCount[8]
#define lbnc  localCellCount[9]

int64_t *tlnc;

int scsim;
int bvsim;
int bnsim;

int nscstages; 		/* number of stem cells stages */
double *sctprob; 	/* stem cells stages transition probabilities */
int64_t *nscinst; 	/* local number of stem cells in different stages */
int64_t *gnscinst; 	/* global number of stem cells in different stages */
int64_t localbc;
int64_t globalbc;

int64_t middleCellIdx;
/* Parallelization */

//#pragma pack(1)
struct partData { /* this structure keeps cell data needed in potential computations */
  double x;
  double y;
  double z;
  double h;
  double size;
  double young;
  int ctype;
};

//#pragma pack(1)
struct densPotData { /* this structure keeps additional cell data (potential & density) */
  double v;
  double density;
};

//#define MAX_CELLS_PER_PROC 10485760
int maxCellsPerProc;

int MPIrank;                    /* MPI rank */
int MPIsize;                    /* MPI size */
int MPIdim[3];                  /* processor topology dimensions (MPI_Dims_create) */
int OMPthreads;                 /* number of OpenMP threads in use */

int MPINodeRank;
int MPINodeSize;
int memPerProc;

MPI_Comm MPI_CART_COMM;
int **MPIcoords;

struct Zoltan_Struct *ztn;

struct partData *sendData;
MIC_ATTR struct partData *recvData;
struct densPotData *sendDensPotData;
MIC_ATTR struct densPotData *recvDensPotData;

int numExp;
MIC_ATTR int numImp;

/* system */
int endian;		/* =0 - big endian, =1 - little endian */

/* model setup */
int MIC_ATTR sdim; 		/* dimensionality of the system */
int mitrand; 		/* mitosis random direction */
int MIC_ATTR nx; 		/* box x size */
int ny; 		/* box y size */
int nz; 		/* box z size */
char rstFileName[128]; 	/* restart file name */
char outdir[128];	/* output directory */
char logdir[128];       /* log directory */
char rng[3]; 		/* type of the Random Number Generator */
int nsteps; 		/* number of simulation steps */

/* simulation */
int simStart;  	        /* start simulation flag */
int step; 		/* step number */
float tstep; 		/* time step size */
float simTime;          /* time of the simulation */
float maxSpeed;         /* maximal displacement of cells in a single time step given by fraction of cell size */
float maxSpeedInUnits;  /* maximal displacement in cm/s */
char cOutType[3];
char fOutType[3];
int vtkout;
int povout;
int vnfout;

/* cell cycle */
float g1;               /* mean duration of G1 phase - healthy tissue */
float s;                /* mean duration of S phase - healthy tissue */
float g2;               /* mean duration of G2 phase - healthy tissue */
float m;                /* mean duration of M phase - healthy tissue */
float v;                /* variability of duration of cell cycles */
float rd;               /* random death probability */

float cg1;              /* mean duration of G1 phase - cancer cells */
float cs;               /* mean duration of S phase - cancer cells */
float cg2;              /* mean duration of G2 phase - cancer cells */
float cm;               /* mean duration of M phase - cancer cells */

double MIC_ATTR csize;           /* cell initial size, no units */
double csizeInUnits;    /* cell size in micrometers */
double cellVolume;      /* cell volume */
double MIC_ATTR h;               /* cell neighbourhood radius */
double MIC_ATTR h2;              /* 2nd power of h */
double MIC_ATTR h3;              /* 3rd power of h */
double MIC_ATTR h4;              /* 4th power of h */

int cancer;
int64_t rsum;

double densityCriticalLevel1;
double densityCriticalLevel2;

int rst;
int rstReset;

int statOutStep;
int rstOutStep;
int vtkOutStep;

int64_t nhs;            /* number of cells to activate random dying - homeostasis of cell colony */

int tgs;		/* - tumor growth simulation, 0 - no tumor growth */

struct doubleVector3d {
  double x;
  double y;
  double z;
};

struct int64Vector3d {
  int64_t x;
  int64_t y;
  int64_t z;
};

struct floatVector3d {
  float x;
  float y;
  float z;
};

struct uintVector3d {
  unsigned int x;
  unsigned int y;
  unsigned int z;
};

struct uintVector3d *locCode;

int64_t localID;

float dummy; /* dummy float parameter in the restart file (it can be used if necessary) */

double boxmin[3],boxmax[3];
double boxsize;

float secondsPerStep;

float gfDt;
float gfH;

int gfIter;
int gfIterPerStep;

/* statistics */
struct statisticsData {
  double minsize; /* Minimum size of cells */
  double maxsize; /* Maximum size of cells */
  double mindist; /* Minimum distance between neighbourhood cells */
  double maxvel;  /* Maximum velocity in the system */
  double minvel;  /* Minimum velocity in the system */
  double maxdens; /* Maximum density */
  double mindens; /* Minimum density */
  double densdev; /* Density deviation */
  double densavg; /* Average density */
};

struct statisticsData MIC_ATTR statistics;

double globalMinVel;
double globalMaxVel;

/* randomization */
#define SEED 985456376
int *stream;


#define N_LEVELS 30
#define ROOT_LEVEL N_LEVELS-1
#define MAXVALUE powf(2,ROOT_LEVEL)

typedef struct _octNode {
  unsigned int xcode;
  unsigned int ycode;
  unsigned int zcode;
  unsigned int xlimit;
  unsigned int ylimit;
  unsigned int zlimit;
  unsigned int level;
  int64_t father;
  int64_t child[8];
  int data;
} octNode;

octNode MIC_ATTR *octree;
int64_t octSize;

/* properties of the affine transformation */
struct doubleVector3d MIC_ATTR affShift;
double MIC_ATTR affScale;
int64_t root;

/*typedef struct _heap {
  int size;
  int count;
  struct bht_node **data;
} heap;
*/
typedef struct _octHeap {
  int size;
  int count;
  int *data;
} octHeap;

int MIC_ATTR tnc;

int ni;
#endif
