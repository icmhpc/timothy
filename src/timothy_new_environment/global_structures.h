//
// Created by Grzegorz Bokota on 10.06.2016.
//


#ifndef TIMOTHY_GLOBAL_STRUCTURES_H
#define TIMOTHY_GLOBAL_STRUCTURES_H

#include <mpi.h>
#include <stdint.h>
#ifdef __MIC__
#define MIC_ATTR __attribute__((target(mic)))
#else
#define MIC_ATTR
#endif

struct environment
{
  double diffusion_coefficient;
  double boundary_condition;
  double initial_condition_mean;
  double initial_condition_variance;
  double lambda_delay;
  double critical_level_1;
  double critical_level_2;
};
/* NEW */

struct doubleVector3d {
  double x;
  double y;
  double z;
};

#ifndef __cplusplus
typedef struct doubleVector3d  doubleVector3d;
#endif

struct int64Vector3d {
  int64_t x;
  int64_t y;
  int64_t z;
};

struct intVector3d {
  int x;
  int y;
  int z;
};

struct str_uint16_pair{
  char * first;
  uint16_t second;
};
typedef struct str_uint16_pair str_uint16_pair;

struct str_uint16_dict{
  str_uint16_pair * members;
  size_t size;
};

typedef struct str_uint16_dict str_uint16_dict;

struct cellCountInfo{
  uint64_t number_of_cells;
  uint64_t g0_phase_number_of_cells;
  uint64_t g1_phase_number_of_cells;
  uint64_t s_phase_number_of_cells;
  uint64_t g2_phase_number_of_cells;
  uint64_t m_phase_number_of_cells;
  uint64_t necrotic_phase_number_of_cells;
  //TODO check names, all needed?
};
enum cell_phase{
  G0_phase,
  G1_phase,
  S_phase,
  G2_phase,
  M_phase,
  Necrotic_phase,
  const_cell,
  cell_to_destroy
};
#define ZOLTAN_ID_TYPE int
struct cellData {
  int lifetime;         /* age of the cell */
  enum cell_phase phase;            /* actual phase of the cell (0=G0,1=G1,2=S,3=G2,4=M,5=Necrotic) */
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
  double h;             /* neighbourhood of the cell  */
  double v;             /* particle potential */
  double density;       /* particle density */
  double scalarField;   /* additional scalar field which can be used for analysis of results (printed to the output VTK files) */
  int cell_type;		/* cell type 1=endothelial */
  int scstage;          /* stem cells stage (in the course of differentiation) */
  unsigned char tumor;	/* 1 - tumor cell, 0 - healthly cell */
};

struct cellsInfo{
  struct cellCountInfo localCellCount;
  struct cellCountInfo totalCellCount;
  uint64_t * localTypeCellCount;
  uint64_t * totalTypeCellCount;
  uint64_t * numberOfCellsInEachProcess;
  struct cellData * cells;
  double ** cellFields;
  size_t numberOfCellTypes;
  struct doubleVector3d *velocity;
  struct cellTypeData * cellTypes;
  str_uint16_dict * cellTypeNumberDict;
  double minimal_density_to_keep_necrotic_cell;
};

struct settings{
  int64_t maxCells;	/* maximal number of cells (set in parameter file) */

  bool  vtkout; //TODO [czaki] setings?
  bool povout; //TODO [czaki] setings?
  bool vnfout; //TODO [czaki] setings?

  int MPI_group_size;
  size_t numberOfCellTypes;
  size_t numberOfEnvironments;
  uint64_t maxCellsPerProc;
  str_uint16_dict * envirormentNumberDict;
  str_uint16_dict * cellTypeNumberDict;
  struct cellTypeData * cellTypes;
  struct environment * environments;
  float global_fields_time_delta;
  //char dimension;

  bool enable_step_transformation;
  void (*step_transformation)(struct cellsInfo * ci, const struct settings * s);
  int numOfDimensions;
  int size_x;
  int size_y;
  int size_z;
  double neighbourhood;
  int MIC_ATTR treeNumberOfChildren;
  int number_of_steps;
  ZOLTAN_ID_TYPE id_range; //Calculate range based on number of mpi process

  bool mitosis_random_direction;
  float random_death;               /* random death probability */
  float minimal_number_of_cells_for_random_death;
  float secondsPerStep;
  int statOutStep;
  int rstOutStep;
  int vtkOutStep;
  float maxSpeed;
  int MIC_ATTR dimension;
  float gfDt;
  float gfH;
  int gfIterPerStep;


  float maxSpeedInUnits;  /* maximal displacement in cm/s */

};

struct state {
  int MPIrank;                    /* MPI rank */
  //unsigned int uMPIrank;                    /* MPI rank */
  int MPIsize;                    /* MPI size */
  unsigned int uMPIsize;           /* MPI size - unsigned version */
  int OMPthreads;
  int MPIdim[3];
  int MPINodeRank;
  int MPINodeSize;
  int memPerProc;
  MPI_Comm MPI_CART_COMM;
  int **MPIcoords;
  int *MPIreverseCords;
  int rstReset;
  ZOLTAN_ID_TYPE localID; //set to 1 on begin
  int simStart;  	        /* start simulation flag */
  int step; 		/* step number */
};

static inline int * getProcesNum(struct state *st, struct intVector3d p){
  return &st->MPIreverseCords[p.x*st->MPIdim[1]*st->MPIdim[2] + p.y * st->MPIdim[2] +p.z];
}

#undef ZOLTAN_ID_TYPE

#endif //TIMOTHY_GLOBAL_STRUCTURES_H
