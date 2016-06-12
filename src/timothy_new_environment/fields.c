//
// Created by czaki on 03.06.16.
//

#include "fields.h"
#include <mpi.h>
#include <stdlib.h>

struct exchangeBuffer {
  double *send;
  double *recv;
  int address;
};

struct buffers {
  struct exchangeBuffer x0;
  struct exchangeBuffer y0;
  struct exchangeBuffer z0;
  struct exchangeBuffer x1;
  struct exchangeBuffer y1;
  struct exchangeBuffer z1;
};

MPI_Request reqFGSend[6];
MPI_Request reqFGRecv[6];

struct buffers buff;
extern double **fieldAddr;
struct doubleVector3d **gradAddr;

//TODO remove while mergigng
void stopRun(int ierr __attribute__((unused)),
             const char *name __attribute__((unused)),
             const char *file __attribute__((unused)),
             int line __attribute__((unused))) {
  exit(1);
}

/*!
 * allocate memory for gradient.  
 */
void fieldsInit(const struct settings *set, const struct state *sta __attribute__((unused)),
                const struct gridData *grid) {
  gradAddr = (struct doubleVector3d **)calloc(set->numberOfEnvironments,
                                              sizeof(struct doubleVector3d *));
  for (size_t i = 0; i < set->numberOfEnvironments; ++i) {
    gradAddr[i] = calloc((size_t)grid->local_size.x * grid->local_size.y *
                             grid->local_size.z,
                         sizeof(struct doubleVector3d));
  }
}

/*!
 * allocate buffers needed to gradient calculation
 */
void allocateFieldGradient(const struct settings *set, const struct state *sta,
                           const struct gridData *grid) {
  if (set->numberOfEnvironments == 0)
    return;
  size_t buffer_size_x =
      grid->local_size.y * grid->local_size.z * set->numberOfEnvironments;
  size_t buffer_size_y =
      grid->local_size.x * grid->local_size.z * set->numberOfEnvironments;
  size_t buffer_size_z =
      grid->local_size.x * grid->local_size.y * set->numberOfEnvironments;
  MPI_Cart_shift(sta->MPI_CART_COMM, 0, 1, &buff.x0.address, &buff.x1.address);
  MPI_Cart_shift(sta->MPI_CART_COMM, 1, 1, &buff.y0.address, &buff.y1.address);
  MPI_Cart_shift(sta->MPI_CART_COMM, 2, 1, &buff.z0.address, &buff.z1.address);

  /* allocate send buffers */
  if (buff.x0.address != MPI_PROC_NULL)
    buff.x0.send = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.x1.address != MPI_PROC_NULL)
    buff.x1.send = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.y0.address != MPI_PROC_NULL)
    buff.y0.send = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.y1.address != MPI_PROC_NULL)
    buff.y1.send = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.z0.address != MPI_PROC_NULL)
    buff.z0.send = (double *)calloc(buffer_size_z, sizeof(double));
  if (buff.z1.address != MPI_PROC_NULL)
    buff.z1.send = (double *)calloc(buffer_size_z, sizeof(double));

  /* allocate receive buffers */
  if (buff.x0.address != MPI_PROC_NULL)
    buff.x0.recv = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.x1.address != MPI_PROC_NULL)
    buff.x1.recv = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.y0.address != MPI_PROC_NULL)
    buff.y0.recv = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.y1.address != MPI_PROC_NULL)
    buff.y1.recv = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.z0.address != MPI_PROC_NULL)
    buff.z0.recv = (double *)calloc(buffer_size_z, sizeof(double));
  if (buff.z1.address != MPI_PROC_NULL)
    buff.z1.recv = (double *)calloc(buffer_size_z, sizeof(double));

  return;
}

/*!
 * Cleaning all buffers used in gradient calculation
 */
void freeFieldGradient() {
  if (buff.x0.address != MPI_PROC_NULL) {
    free(buff.x0.send);
    free(buff.x0.recv);
  }
  if (buff.x1.address != MPI_PROC_NULL) {
    free(buff.x1.send);
    free(buff.x1.recv);
  }
  if (buff.y0.address != MPI_PROC_NULL) {
    free(buff.y0.send);
    free(buff.y0.recv);
  }
  if (buff.y1.address != MPI_PROC_NULL) {
    free(buff.y1.send);
    free(buff.y1.recv);
  }
  if (buff.z0.address != MPI_PROC_NULL) {
    free(buff.z0.send);
    free(buff.z0.recv);
  }
  if (buff.z1.address != MPI_PROC_NULL) {
    free(buff.z1.send);
    free(buff.z1.recv);
  }
}

/*!
 * This function start exchange data needed to calculate gradient of
 * environment elements
 */
void initFieldHaloExchange(const struct state *sta, const struct settings *set,
                           const struct gridData *grid) {
  /* function is longer than original to avoid reading fields inside box */
  if (set->numberOfEnvironments == 0)
    return;

  struct intVector3d gridSize = grid->local_size;
  int64_t buffer_x_step, buffer_y_step, buffer_en_step;
  size_t numberOfEnvironments;
  numberOfEnvironments = set->numberOfEnvironments;
  buffer_x_step = gridSize.y * gridSize.z;
  buffer_y_step = gridSize.z;
  /* copy to x0 buffer */
  if (buff.x0.address != MPI_PROC_NULL) {
    buffer_en_step = gridSize.y * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t y = 0; y < gridSize.y; y++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.x0.send[buffer_en_step * en + gridSize.z * y + z] =
              fieldAddr[en][buffer_y_step * y + z];
        }
      }
    }
    MPI_Isend(buff.x0.send, gridSize.y * gridSize.z, MPI_DOUBLE,
              buff.x0.address, sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM,
              &reqFGSend[0]);
    MPI_Irecv(buff.x0.recv, gridSize.y * gridSize.z, MPI_DOUBLE,
              buff.x0.address, sta->MPIsize + buff.x0.address,
              sta->MPI_CART_COMM, &reqFGRecv[0]);
  }
  /* copy to x1 buffer */
  if (buff.x1.address != MPI_PROC_NULL) {
    buffer_en_step = gridSize.y * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t y = 0; y < gridSize.y; y++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.x1.send[buffer_en_step * en + gridSize.z * y + z] =
              fieldAddr[en]
                       [buffer_x_step * (gridSize.x - 1) + gridSize.z * y + z];
        }
      }
    }
    MPI_Isend(buff.x1.send, gridSize.y * gridSize.z, MPI_DOUBLE,
              buff.x1.address, sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM,
              &reqFGSend[1]);
    MPI_Irecv(buff.x1.recv, gridSize.y * gridSize.z, MPI_DOUBLE,
              buff.x1.address, sta->MPIsize + buff.x1.address,
              sta->MPI_CART_COMM, &reqFGRecv[1]);
  }
  /* copy to y0 buffer */
  if (buff.y0.address != MPI_PROC_NULL) {
    buffer_en_step = gridSize.x * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t x = 0; x < gridSize.y; x++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.y0.send[buffer_en_step * en + gridSize.z * x + z] =
              fieldAddr[en][buffer_x_step * x + z];
        }
      }
    }
    MPI_Isend(buff.y0.send, gridSize.x * gridSize.z, MPI_DOUBLE,
              buff.y0.address, sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM,
              &reqFGSend[2]);
    MPI_Irecv(buff.y0.recv, gridSize.x * gridSize.z, MPI_DOUBLE,
              buff.y0.address, sta->MPIsize + buff.y0.address,
              sta->MPI_CART_COMM, &reqFGRecv[2]);
  }
  /* copy to y1 buffer */
  if (buff.y1.address != MPI_PROC_NULL) {
    buffer_en_step = gridSize.x * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t x = 0; x < gridSize.y; x++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.y1.send[buffer_en_step * en + gridSize.z * x + z] =
              fieldAddr[en][buffer_x_step * x +
                            buffer_y_step * (gridSize.y - 1) + z];
        }
      }
    }
    MPI_Isend(buff.y1.send, gridSize.x * gridSize.z, MPI_DOUBLE,
              buff.y1.address, sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM,
              &reqFGSend[3]);
    MPI_Irecv(buff.y1.recv, gridSize.x * gridSize.z, MPI_DOUBLE,
              buff.y1.address, sta->MPIsize + buff.y1.address,
              sta->MPI_CART_COMM, &reqFGRecv[3]);
  }
  /* copy to z0 buffer */
  if (buff.z0.address != MPI_PROC_NULL) {
    buffer_en_step = gridSize.x * gridSize.y;
    buffer_x_step = gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t x = 0; x < gridSize.y; x++) {
        for (int64_t y = 0; y < gridSize.z; y++) {
          buff.y0.send[buffer_en_step * en + gridSize.y * x + y] =
              fieldAddr[en][buffer_x_step * x + buffer_y_step * y];
        }
      }
    }
    MPI_Isend(buff.z0.send, gridSize.x * gridSize.y, MPI_DOUBLE,
              buff.z0.address, sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM,
              &reqFGSend[4]);
    MPI_Irecv(buff.z0.recv, gridSize.x * gridSize.y, MPI_DOUBLE,
              buff.z0.address, sta->MPIsize + buff.z0.address,
              sta->MPI_CART_COMM, &reqFGRecv[4]);
  }
  /* copy to z1 buffer */
  if (buff.z1.address != MPI_PROC_NULL) {
    buffer_en_step = gridSize.x * gridSize.y;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t x = 0; x < gridSize.y; x++) {
        for (int64_t y = 0; y < gridSize.z; y++) {
          buff.y1.send[buffer_en_step * en + gridSize.y * x + y] =
              fieldAddr[en][buffer_x_step * x + buffer_y_step * y +
                            (gridSize.z - 1)];
        }
      }
    }
    MPI_Isend(buff.z1.send, gridSize.x * gridSize.y, MPI_DOUBLE,
              buff.z1.address, sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM,
              &reqFGSend[5]);
    MPI_Irecv(buff.z1.recv, gridSize.x * gridSize.y, MPI_DOUBLE,
              buff.z1.address, sta->MPIsize + buff.z1.address,
              sta->MPI_CART_COMM, &reqFGRecv[5]);
  }
}

void waitFieldHaloRecv() {
  MPI_Status status;
  if (buff.x0.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[0], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.x1.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[1], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.y0.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[2], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.y1.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[3], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.z0.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[4], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.z1.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[5], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  return;
}

void waitFieldHaloSend() {
  MPI_Status status;
  if (buff.x0.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[0], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.x1.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[1], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.y0.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[2], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.y1.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[3], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.z0.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[4], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if (buff.z1.address != MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[5], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  return;
}

/*!
 * This function compute gradients of envirorments elements
 */
void computeFieldGradient(const struct state *sta, const struct settings *set,
                          const struct gridData *grid) {
  initFieldHaloExchange(sta, set, grid);
  struct intVector3d gridSize = grid->local_size;
  int x_step = gridSize.y * gridSize.z;
  int y_step = gridSize.z;
  for (size_t en = 0; en < set->numberOfEnvironments; en++) {
    for (int x = 0; x < gridSize.x; ++x) {
      for (int y = 0; y < gridSize.y; ++y) {
        for (int z = 0; z < gridSize.z; ++z) {
          /* x coord gradient */
          if (x != 0 && x != gridSize.x - 1)
            gradAddr[en][x_step * x + y_step * y + z].x =
                (fieldAddr[en][x_step * (x + 1) + y_step * y + z] -
                 fieldAddr[en][x_step * (x - 1) + y_step * y + z]) /
                grid->resolution;
          /* y coord gradient */
          if (y != 0 && y != gridSize.y - 1)
            gradAddr[en][x_step * x + y_step * y + z].y =
                (fieldAddr[en][x_step * x + y_step * (y + 1) + z] -
                 fieldAddr[en][x_step * x + y_step * (y - 1) + z]) /
                grid->resolution;
          /* z coord gradient */
          if (z != 0 && z != gridSize.z - 1)
            gradAddr[en][x_step * x + y_step * y + z].z =
                (fieldAddr[en][x_step * x + y_step * y + z + 1] -
                 fieldAddr[en][x_step * x + y_step * y + z - 1]) /
                grid->resolution;
        }
      }
    }
  }
  waitFieldHaloRecv();
  for (size_t en = 0; en < set->numberOfEnvironments; en++) {
    int buffer_en_step = gridSize.y * gridSize.z;
    for (int y = 0; y < gridSize.y; ++y) {
      for (int z = 0; z < gridSize.z; ++z) {
        double x0, x1;
        if (buff.x0.address != MPI_PROC_NULL)
          x0 = buff.x0.recv[buffer_en_step * en + gridSize.z * y + z];
        else
          x0 = fieldAddr[en][x_step * 0 + y_step * y + z];
        if (buff.x1.address != MPI_PROC_NULL)
          x1 = buff.x1.recv[buffer_en_step * en + gridSize.z * y + z];
        else
          x1 = fieldAddr[en][x_step * (gridSize.x - 1) + y_step * y + z];
        gradAddr[en][x_step * 0 + y_step * y + z].x =
            (fieldAddr[en][x_step * 1 + y_step * y + z] - x0) /
            grid->resolution;
        gradAddr[en][x_step * (gridSize.x - 1) + y_step * y + z].x =
            (x1 - fieldAddr[en][x_step * (gridSize.x - 2) + y_step * y + z]) /
            grid->resolution;
      }
    }
  }
  for (size_t en = 0; en < set->numberOfEnvironments; en++) {
    int buffer_en_step = gridSize.x * gridSize.z;
    for (int x = 0; x < gridSize.y; ++x) {
      for (int z = 0; z < gridSize.z; ++z) {
        double y0, y1;
        if (buff.y0.address != MPI_PROC_NULL)
          y0 = buff.y0.recv[buffer_en_step * en + gridSize.z * x + z];
        else
          y0 = fieldAddr[en][x_step * x + y_step * 0 + z];
        if (buff.y1.address != MPI_PROC_NULL)
          y1 = buff.y1.recv[buffer_en_step * en + gridSize.z * x + z];
        else
          y1 = fieldAddr[en][x_step * x + y_step * (gridSize.y - 1) + z];
        gradAddr[en][x_step * x + y_step * 0 + z].y =
            (fieldAddr[en][x_step * x + y_step * 1 + z] - y0) /
            grid->resolution;
        gradAddr[en][x_step * x + y_step * (gridSize.y - 1) + z].z =
            (y1 - fieldAddr[en][x_step * x + y_step * (gridSize.y - 2) + z]) /
            grid->resolution;
      }
    }
  }
  for (size_t en = 0; en < set->numberOfEnvironments; en++) {
    int buffer_en_step = gridSize.x * gridSize.z;
    for (int x = 0; x < gridSize.y; ++x) {
      for (int y = 0; y < gridSize.z; ++y) {
        double z0, z1;
        if (buff.z0.address != MPI_PROC_NULL)
          z0 = buff.z0.recv[buffer_en_step * en + gridSize.y * x + y];
        else
          z0 = fieldAddr[en][x_step * x + y_step * y + 0];
        if (buff.z1.address != MPI_PROC_NULL)
          z1 = buff.z1.recv[buffer_en_step * en + gridSize.y * x + y];
        else
          z1 = fieldAddr[en][x_step * x + y_step * y + (gridSize.z - 1)];
        gradAddr[en][x_step * x + y_step * y].z =
            (fieldAddr[en][x_step * x + y_step * y + 1] - z0) /
            grid->resolution;
        gradAddr[en][x_step * x + y_step * y + (gridSize.z - 1)].z =
            (z1 - fieldAddr[en][x_step * x + y_step * y + (gridSize.z - 2)]) /
            grid->resolution;
      }
    }
  }

  waitFieldHaloSend();
}