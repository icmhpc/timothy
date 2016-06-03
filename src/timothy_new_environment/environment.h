#ifndef TIMOTHY_ENVIRONMENT_H
#define TIMOTHY_ENVIRONMENT_H
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <mpi.h>
#ifdef __MIC__
#define MIC_ATTR __attribute__((target(mic)))
#else
#define MIC_ATTR
#endif

/* NEW */
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

typedef struct _doubleVector3d {
  double x;
  double y;
  double z;
} doubleVector3d;

struct int64Vector3d {
  int64_t x;
  int64_t y;
  int64_t z;
};

/* NEW */
struct gridData{
  struct int64Vector3d global_size;
  doubleVector3d lower_corner,upper_corner;
  struct int64Vector3d local_size;
  double resolution;
  double voxel_volume;
  struct int64Vector3d *lower_indices,*upper_indices;
  doubleVector3d *buffer;
};
/* NEW */

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

#define ZOLTAN_ID_TYPE int

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
    int rstReset;
    ZOLTAN_ID_TYPE localID; //set to 1 on begin
    int simStart;  	        /* start simulation flag */
    int step; 		/* step number */
};
#endif //TIMOTHY_ENVIRONMENT_H
