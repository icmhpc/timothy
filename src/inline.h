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

/*! \file inline.h
 *  \brief contains various Timothy inline functions
 */
#ifndef TIMOTHY_INLINE_H
#define TIMOTHY_INLINE_H
#include <math.h>
#include "global.h"
MIC_ATTR static inline int octNodeIntersection(int idx,struct uintVector3d minLocCode,struct uintVector3d maxLocCode)
{
  if( maxLocCode.x<=octree[idx].xcode || minLocCode.x>=octree[idx].xlimit ) return 0;
  if( maxLocCode.y<=octree[idx].ycode || minLocCode.y>=octree[idx].ylimit ) return 0;
  if( maxLocCode.z<=octree[idx].zcode || minLocCode.z>=octree[idx].zlimit ) return 0;
  return 1;
}

/*!
 * This function returns the value of the SPH kernel function for given two cells p1 and p2.
 * Kernel function as in:
 * - J.J.Monaghan,"Smoothed Particle Hydrodynamics",Annu.Rev.Astron.Astrophys.,1992,30:543-74
 * - V.Springel,"Smoothed Particle Hydrodynamics in Astrophysics",arXiv:1109.2219v1
 */
MIC_ATTR static inline float sph_kernel(double r)
{
  double u,c=1.0;

  if(r<0.0) return 0.0;

  u=r/h;
  if(sdim==2) c=40/(7*M_PI*h2);
  if(sdim==3) c=8/(M_PI*h3);

  if(u>=0.0 && u<=0.5) return c*(1-6*u*u+6*u*u*u);
  if(u>0.5 && u<=1.0) return c*(2*(1-u)*(1-u)*(1-u));
  if(u>1.0) return 0.0;
  return 0.0;
}

/*!
 * This function return the value of the SPH kernel function gradient.
 * mode=0 - computations for local particles
 * mode=1 - computations for remote particles
 * Kernel function as in
 * - J.J.Monaghan,"Smoothed Particle Hydrodynamics",Annu.Rev.Astron.Astrophys.,1992,30:543-74
 * - V.Springel,"Smoothed Particle Hydrodynamics in Astrophysics",arXiv:1109.2219v1
 */
MIC_ATTR static inline int sph_kernel_gradient(int p1, int p2, double grad[3],int mode,double r)
{
  double u=1.0,c=1.0,w=1.0;
  double x1,x2,y1,y2,z1,z2;

  if(mode==0) {
    x1=cells[p1].x;
    x2=cells[p2].x;
    y1=cells[p1].y;
    y2=cells[p2].y;
    z1=cells[p1].z;
    z2=cells[p2].z;
  } else {
    x1=recvData[p1].x;
    x2=cells[p2].x;
    y1=recvData[p1].y;
    y2=cells[p2].y;
    z1=recvData[p1].z;
    z2=cells[p2].z;
  }

  if(r>=0.0 && r<=h) {

    u=r/h;
    if(sdim==2) c=(40*8)/(7*M_PI*h3);
    if(sdim==3) c=48/(M_PI*h4);

    if(u>=0.0 && u<=0.5) w=-2.0*u+3.0*u*u;
    if(u>0.5 && u<=1.0) w=-(1.0-u)*(1.0-u);
    if(u>1.0) w=0.0;

    grad[0]=w*c*(x2-x1)/r;
    grad[1]=w*c*(y2-y1)/r;
    if(sdim==3) grad[2]=w*c*(z2-z1)/r;
    if(sdim==2) grad[2]=0.0;

  } else {

    grad[0]=0.0;
    grad[1]=0.0;
    grad[2]=0.0;

  }
  return 0;
}

#endif /* ifndef TIMOTHY_INLINE_H  */
