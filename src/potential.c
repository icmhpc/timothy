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

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "global.h"
#include "octree.h"
#include "fields.h"
#include "inline.h"

/*! \file potential.c
 *  \brief contains functions that compute the potential
 */

/*!
 * This function computes potential for two neighbour cells.
 * mode=0 - computations for two local cells
 * mode=1 - computations for one local and one remote cell
 */
MIC_ATTR double ccPot(int p1, int p2, int mode,double *mindist)
{
  double pot;
  double dist=0.0;
  double size;
  double x, y, z;
  double xc;
  double D;
  double poisson = 0.33;
  double young;
  int ctype;


  if (mode == 0 && p1 == p2)
    return 0.0;

  if (mode == 0) {
    x = cells[p1].x;
    y = cells[p1].y;
    z = cells[p1].z;
    size = cells[p1].size;
    young = cells[p1].young;
    ctype = cells[p1].ctype;
  }
  if (mode == 1) {
    x = recvData[p1].x;
    y = recvData[p1].y;
    z = recvData[p1].z;
    size = recvData[p1].size;
    young = recvData[p1].young;
    ctype = recvData[p1].ctype;
  }

  /* compute the distance between two cells */
  if (sdim == 2)
    dist =
      sqrt((x - cells[p2].x) * (x - cells[p2].x) +
           (y - cells[p2].y) * (y - cells[p2].y));
  if (sdim == 3)
    dist =
      sqrt((x - cells[p2].x) * (x - cells[p2].x) +
           (y - cells[p2].y) * (y - cells[p2].y) +
           (z - cells[p2].z) * (z - cells[p2].z));

  if (mindist[0] > dist && ctype!=1 && cells[p2].ctype!=1) {
    mindist[0] = dist;
  }

  if (dist <= h) {

    double r01, r02;
    double area;
    double sc=1.0;

    if (sdim == 2)
      sc = h2;
    if (sdim == 3)
      sc = h3;

    if (mode == 0) {
      cells[p1].density +=
        sc * (cells[p2].size / csize) * sph_kernel(dist);
    }
    if (mode == 1) {
      #pragma omp atomic
      cells[p2].density += sc * (size / csize) * sph_kernel(dist);
    }

    xc = size + cells[p2].size - dist;

    if (xc <= 0.0)
      return 0.0;

    D = 0.75 * ((1.0 - poisson * poisson) / young +
                (1.0 - poisson * poisson / cells[p2].young));

    /* adhesion */
    r01 =
      (size * size - cells[p2].size * cells[p2].size +
       dist * dist) / (2 * dist);
    r02 = dist - r01;

    area =
      M_PI *
      ((size * size * (size - r01) -
        (size * size * size - r01 * r01 * r01) / 3) +
       (cells[p2].size * cells[p2].size * (cells[p2].size - r02) -
        (cells[p2].size * cells[p2].size * cells[p2].size -
         r02 * r02 * r02) / 3));

    /* compute potential */
    pot =
      (2.0 * pow(xc, 5 / 2) / (5.0 * D)) * sqrt((size * cells[p2].size) /
          (size +
           cells[p2].size)) +
      area * 0.1;

    if(ctype==1 && cells[p2].ctype==1) pot = 0.0;

    return pot;

  } else
    return 0.0;

}

/*!
 * This function implements tree walk algorithm for each local cell.
 * Function ccPot(..) is called for each pair of neighbours.
 */
MIC_ATTR void compPot()
{
  if(lnc<=1) return;
  #pragma omp parallel
  {
    int p;
    double mindist;
    mindist=nx;

    #pragma omp for schedule(dynamic,64)
    for (p = 0; p < lnc; p++) {
      int64_t cellIdx;
      int newIdx;
      int s;
      int octIndex;
      struct uintVector3d minLocCode,maxLocCode;
      octHeap octh;
      octHeapInit(&octh);

      cells[p].density = 0.0;
      cells[p].v = 0.0;

      octComputeBox(p,&minLocCode,&maxLocCode);
      octIndex=octLocateRegion(minLocCode,maxLocCode);

      for(s=0; s<tnc; s++) {
        int idx;
        idx=octree[octIndex].child[s];
        if(idx!=-1 && octNodeIntersection(idx,minLocCode,maxLocCode))
          octHeapPush(&octh,idx);
      }

      while(octh.count>0) {
        int idx;
        idx=octHeapPop(&octh);
        cellIdx=octree[idx].data;
        if(cellIdx>=0) {
          cells[p].v += ccPot(p,cellIdx,0,&mindist);
        } else {
          for(s=0; s<tnc; s++) {
            newIdx=octree[idx].child[s];
            if(newIdx!=-1 && octNodeIntersection(newIdx,minLocCode,maxLocCode))
              octHeapPush(&octh,newIdx);
          }
        }
      }
      octHeapFree(&octh);
    }
    #pragma omp critical
    statistics.mindist=(statistics.mindist>mindist?mindist:statistics.mindist);
  }
}

/*!
 * This function implements tree walk algorithm for each remote cell.
 * Function ccPot(..) is called for each pair of neighbours.
 */
MIC_ATTR void compRPot()
{
  int rp;
  if(numImp<=0) return;
  #pragma omp parallel
  {
    double mindist;
    mindist=nx;
    #pragma omp for schedule(dynamic,64)
    for (rp = 0; rp < numImp; rp++) {
      int64_t cellIdx;
      int newIdx;
      int s;
      double v;
      struct uintVector3d minLocCode,maxLocCode;
      octHeap octh;
      octHeapInit(&octh);
      octComputeBoxR(rp,&minLocCode,&maxLocCode);
      octHeapPush(&octh,0);
      while(octh.count>0) {
        int idx;
        idx=octHeapPop(&octh);
        cellIdx=octree[idx].data;
        if(cellIdx>=0) {
          v=ccPot(rp,cellIdx,1,&mindist);
          #pragma omp atomic
          cells[cellIdx].v += v;
        } else {
          for(s=0; s<tnc; s++) {
            newIdx=octree[idx].child[s];
            if(newIdx!=-1 && octNodeIntersection(newIdx,minLocCode,maxLocCode))
              octHeapPush(&octh,newIdx);
          }
        }
      }
      octHeapFree(&octh);
    }
    #pragma omp critical
    statistics.mindist=(statistics.mindist>mindist?mindist:statistics.mindist);
  }
}

/*!
 *  This function computes potential gradient for two neighbour cells.
 * mode=0 - computations for local cells
 * mode=1 - computations for remote cells
 */
MIC_ATTR void ccPotGrad(int p1, int p2, int mode)
{
  double grad[3];
  /* we assume that cells' mass is always 1.0 */
  double m = 1.0;
  double sc;
  double x1, x2, y1, y2, z1, z2;
  double v, density, size;
  double dist;
  int ctype; // TODO Check

  if (p1 == p2 && mode == 0)
    return;

  if (mode == 0) {
    x1 = cells[p1].x;
    x2 = cells[p2].x;
    y1 = cells[p1].y;
    y2 = cells[p2].y;
    z1 = cells[p1].z;
    z2 = cells[p2].z;
    v = cells[p2].v;
    density = cells[p2].density;
    size = cells[p2].size;
    ctype = cells[p1].ctype;
  }
  if (mode == 1) {
    x1 = recvData[p1].x;
    x2 = cells[p2].x;
    y1 = recvData[p1].y;
    y2 = cells[p2].y;
    z1 = recvData[p1].z;
    z2 = cells[p2].z;
    v = recvDensPotData[p1].v;
    density = recvDensPotData[p1].density;
    size = recvData[p1].size;
    ctype = recvData[p1].ctype;
  } // Is possible to have more modes, why not else?
  
  dist = 1; // TODO check
  if (sdim == 2)
    dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  if (sdim == 3)
    dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) +
                (z1 - z2) * (z1 - z2));

  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = 0.0;

  /* compute the gradient of SPH kernel function */
  sph_kernel_gradient(p1, p2, grad, mode, dist);

  if (density == 0.0)
    sc = 0.0;
  else {
    m = size / csize;
    sc = (m / density) * v;
  }

  /* update forces */
  if (mode == 0) {
    velocity[p1].x += sc * grad[0];
    velocity[p1].y += sc * grad[1];
    velocity[p1].z += sc * grad[2];
  } else {
    #pragma omp atomic
    velocity[p2].x -= sc * grad[0];
    #pragma omp atomic
    velocity[p2].y -= sc * grad[1];
    #pragma omp atomic
    velocity[p2].z -= sc * grad[2];
  }
}

/*!
 * This function implements tree walk algorithm for each local cell to compute potential gradient.
 * Function ccPotGrad(..) is called for each pair of neighbours.
 */
MIC_ATTR void compPotGrad()
{
  int p;
  if(lnc<=1) return;
  #pragma omp parallel for schedule(dynamic,64)
  for(p=0; p<lnc; p++) {
    int64_t cellIdx;
    int newIdx;
    int s;
    int octIndex;
    struct uintVector3d minLocCode,maxLocCode;
    octHeap octh;

    octHeapInit(&octh);
    velocity[p].x=0.0;
    velocity[p].y=0.0;
    velocity[p].z=0.0;

    octComputeBox(p,&minLocCode,&maxLocCode);
    octIndex=octLocateRegion(minLocCode,maxLocCode);

    for(s=0; s<tnc; s++) {
      int idx;
      idx=octree[octIndex].child[s];
      if(idx!=-1 && octNodeIntersection(idx,minLocCode,maxLocCode))
        octHeapPush(&octh,idx);
    }

    while(octh.count>0) {
      int idx;
      idx=octHeapPop(&octh);
      cellIdx=octree[idx].data;
      if(cellIdx>=0) {
        ccPotGrad(p,cellIdx,0);
      } else {
        for(s=0; s<tnc; s++) {
          newIdx=octree[idx].child[s];
          if(newIdx!=-1 && octNodeIntersection(newIdx,minLocCode,maxLocCode))
            octHeapPush(&octh,newIdx);
        }
      }
    }
    octHeapFree(&octh);
  }
}

/*!
 * This function implements tree walk algorithm for each remote cell to compute potential gradient.
 * Function ccPotGrad(..) is called for each pair of neighbours.
 */
MIC_ATTR void compRPotGrad()
{
  int rp;
  if(numImp<=0) return;
  #pragma omp parallel for schedule(dynamic,64)
  for (rp = 0; rp < numImp; rp++) {
    int64_t cellIdx;
    int newIdx;
    int s;
    struct uintVector3d minLocCode,maxLocCode;
    octHeap octh;
    octHeapInit(&octh);
    octComputeBoxR(rp,&minLocCode,&maxLocCode);
    octHeapPush(&octh,0);
    while(octh.count>0) {
      int idx;
      idx=octHeapPop(&octh);
      cellIdx=octree[idx].data;
      if(cellIdx>=0) {
        ccPotGrad(rp,cellIdx,1);
      } else {
        for(s=0; s<tnc; s++) {
          newIdx=octree[idx].child[s];
          if(newIdx!=-1 && octNodeIntersection(newIdx,minLocCode,maxLocCode))
            octHeapPush(&octh,newIdx);
        }
      }
    }
    octHeapFree(&octh);
  }
}
