//
// Created by czaki on 03.06.16.
//

#include "fields.h"
#include <mpi.h>
#include <stdlib.h>

struct exchangeBuffer {
  double *send;
  double *recv;
};

struct buffers {
  struct exchangeBuffer x0;
  struct exchangeBuffer y0;
  struct exchangeBuffer z0;
  struct exchangeBuffer x1;
  struct exchangeBuffer y1;
  struct exchangeBuffer z1;
  int rX0, rX1, rY0, rY1, rZ0, rZ1;
};

MPI_Request reqFGSend[6];
MPI_Request reqFGRecv[6];

struct buffers buff;
extern double **fieldAddr;

void stopRun(int ierr __attribute__((unused)), const char *name __attribute__((unused)),
             const char *file __attribute__((unused)), int line  __attribute__((unused))) {
  exit(1);
}

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
  MPI_Cart_shift(sta->MPI_CART_COMM, 0, 1, &buff.rX0, &buff.rX1);
  MPI_Cart_shift(sta->MPI_CART_COMM, 1, 1, &buff.rY0, &buff.rY1);
  MPI_Cart_shift(sta->MPI_CART_COMM, 2, 1, &buff.rZ0, &buff.rZ1);

  /* allocate send buffers */
  if (buff.rX0 != MPI_PROC_NULL)
    buff.x0.send = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.rX1 != MPI_PROC_NULL)
    buff.x1.send = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.rY0 != MPI_PROC_NULL)
    buff.y0.send = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.rY1 != MPI_PROC_NULL)
    buff.y1.send = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.rZ0 != MPI_PROC_NULL)
    buff.z0.send = (double *)calloc(buffer_size_z, sizeof(double));
  if (buff.rZ1 != MPI_PROC_NULL)
    buff.z1.send = (double *)calloc(buffer_size_z, sizeof(double));

  /* allocate receive buffers */
  if (buff.rX0 != MPI_PROC_NULL)
    buff.x0.recv = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.rX1 != MPI_PROC_NULL)
    buff.x1.recv = (double *)calloc(buffer_size_x, sizeof(double));
  if (buff.rY0 != MPI_PROC_NULL)
    buff.y0.recv = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.rY1 != MPI_PROC_NULL)
    buff.y1.recv = (double *)calloc(buffer_size_y, sizeof(double));
  if (buff.rZ0 != MPI_PROC_NULL)
    buff.z0.recv = (double *)calloc(buffer_size_z, sizeof(double));
  if (buff.rZ1 != MPI_PROC_NULL)
    buff.z1.recv = (double *)calloc(buffer_size_z, sizeof(double));

  return;
}

void freeFieldGradient() {
  free(buff.x0.send);
  free(buff.x0.recv);
  free(buff.x1.send);
  free(buff.x1.recv);
  free(buff.y0.send);
  free(buff.y0.recv);
  free(buff.y1.send);
  free(buff.y1.recv);
  free(buff.z0.send);
  free(buff.z0.recv);
  free(buff.z1.send);
  free(buff.z1.recv);
}

void initFieldHaloExchange(const struct state *sta, const struct settings *set,
                           const struct gridData *grid) {
  /* function is longer than original to avoid reading fields inside box */
  if (set->numberOfEnvironments == 0)
    return;

  struct int64Vector3d gridSize = grid->local_size;
  int64_t buffer_x_step, buffer_y_step, buffer_en_step;
  size_t numberOfEnvironments;
  numberOfEnvironments = set->numberOfEnvironments;
  buffer_x_step = gridSize.y * gridSize.z;
  buffer_y_step = gridSize.z;
  /* copy to x0 buffer */
  if (buff.rX0 != MPI_PROC_NULL) {
    buffer_en_step = gridSize.y * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t y = 0; y < gridSize.y; y++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.x0.send[buffer_en_step * en + gridSize.z * y + z] =
              fieldAddr[en][buffer_y_step * y + z];
        }
      }
    }
    MPI_Isend(buff.x0.send, gridSize.y * gridSize.z, MPI_DOUBLE, buff.rX0,
              sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM, &reqFGSend[0]);
    MPI_Irecv(buff.x0.recv, gridSize.y * gridSize.z, MPI_DOUBLE, buff.rX0,
              sta->MPIsize + buff.rX0, sta->MPI_CART_COMM, &reqFGRecv[0]);
  }
  /* copy to x1 buffer */
  if (buff.rX1 != MPI_PROC_NULL) {
    buffer_en_step = gridSize.y * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t y = 0; y < gridSize.y; y++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.x1.send[buffer_en_step * en + gridSize.z * y + z] =
              fieldAddr[en][buffer_x_step * (gridSize.x - 1) + gridSize.z * y +
                            z];
        }
      }
    }
    MPI_Isend(buff.x1.send, gridSize.y * gridSize.z, MPI_DOUBLE, buff.rX1,
              sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM, &reqFGSend[1]);
    MPI_Irecv(buff.x1.recv, gridSize.y * gridSize.z, MPI_DOUBLE, buff.rX1,
              sta->MPIsize + buff.rX1, sta->MPI_CART_COMM, &reqFGRecv[1]);
  }
  /* copy to y0 buffer */
  if (buff.rY0 != MPI_PROC_NULL) {
    buffer_en_step = gridSize.x * gridSize.z;
    for (size_t en = 0; en < numberOfEnvironments; en++) {
      for (int64_t x = 0; x < gridSize.y; x++) {
        for (int64_t z = 0; z < gridSize.z; z++) {
          buff.y0.send[buffer_en_step * en + gridSize.z * x + z] =
              fieldAddr[en][buffer_x_step * x + z];
        }
      }
    }
    MPI_Isend(buff.y0.send, gridSize.x * gridSize.z, MPI_DOUBLE, buff.rY0,
              sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM, &reqFGSend[2]);
    MPI_Irecv(buff.y0.recv, gridSize.x * gridSize.z, MPI_DOUBLE, buff.rY0,
              sta->MPIsize + buff.rY0, sta->MPI_CART_COMM, &reqFGRecv[2]);
  }
  /* copy to y1 buffer */
  if (buff.rY1 != MPI_PROC_NULL) {
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
    MPI_Isend(buff.y1.send, gridSize.x * gridSize.z, MPI_DOUBLE, buff.rY1,
              sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM, &reqFGSend[3]);
    MPI_Irecv(buff.y1.recv, gridSize.x * gridSize.z, MPI_DOUBLE, buff.rY1,
              sta->MPIsize + buff.rY1, sta->MPI_CART_COMM, &reqFGRecv[3]);
  }
  /* copy to z0 buffer */
  if (buff.rZ0 != MPI_PROC_NULL) {
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
    MPI_Isend(buff.z0.send, gridSize.x * gridSize.y, MPI_DOUBLE, buff.rZ0,
              sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM, &reqFGSend[4]);
    MPI_Irecv(buff.z0.recv, gridSize.x * gridSize.y, MPI_DOUBLE, buff.rZ0,
              sta->MPIsize + buff.rZ0, sta->MPI_CART_COMM, &reqFGRecv[4]);
  }
  /* copy to z1 buffer */
  if (buff.rZ1 != MPI_PROC_NULL) {
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
    MPI_Isend(buff.z1.send, gridSize.x * gridSize.y, MPI_DOUBLE, buff.rZ1,
              sta->MPIsize + sta->MPIrank, sta->MPI_CART_COMM, &reqFGSend[5]);
    MPI_Irecv(buff.z1.recv, gridSize.x * gridSize.y, MPI_DOUBLE, buff.rZ1,
              sta->MPIsize + buff.rZ1, sta->MPI_CART_COMM, &reqFGRecv[5]);
  }
}

void waitFieldHaloRecv()
{
  MPI_Status status;
  if(buff.rX0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[0], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rX1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[1], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rY0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[2], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rY1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[3], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rZ0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[4], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rZ1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGRecv[5], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  return;
}

void waitFieldHaloSend()
{
  MPI_Status status;
  if(buff.rX0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[0], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rX1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[1], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rY0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[2], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rY1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[3], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rZ0!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[4], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  if(buff.rZ1!=MPI_PROC_NULL) {
    if (MPI_Wait(&reqFGSend[5], &status) != MPI_SUCCESS)
      stopRun(103, "reqFGSend", __FILE__, __LINE__);
  }
  return;
}


void computeFieldGradient(const struct state *sta, const struct settings *set,
                          const struct gridData *grid) {
  waitFieldHaloRecv();

  waitFieldHaloSend();

}