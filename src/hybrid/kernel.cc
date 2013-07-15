// Copyright (c) 2013, Timothy Roberts and Amanda Ng
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Timothy Roberts and Amanda Ng ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author: timothy.daniel.roberts@gmail.com, amanda.ng@gmail.com

/*! \file kernel.cc
    \brief Kernel class file.

    Implementation of the Kernel class.
*/

#include "hybrid/kernel.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

#include "hybrid/dataspec.h"
#include "hybrid/modelmap.h"

Kernel::Kernel() {
  kernel = NULL;
  nnz = 0;
  skernel = NULL;
  ckernel = NULL;
  ctr = NULL;
  sin2beta = NULL;
  gx = NULL;
  gy = NULL;
  gz = NULL;
}
Kernel::~Kernel() {}
void Kernel::CreateSphericalKernel(const DataSpec &dspec) {
  usedtype SPHa = dspec.B0 / (4 * M_PI);
  usedtype minkernelvalue, maxkernelvalue = 0;

  for (int x = -halfsize; x <= halfsize; x++) {
    for (int y = -halfsize; y <= halfsize; y++) {
      for (int z = -halfsize; z <= halfsize; z++) {
        int p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
        if (x == 0 && y == 0 && z == 0) {
          skernel[p] = 0;
          nnz--;
        } else {
          usedtype rmag = sqrt(x*x + y*y + z*z);
          usedtype rx = x/rmag;
          usedtype ry = y/rmag;
          usedtype rz = z/rmag;
          usedtype br = dspec.bhat[0]*rx + dspec.bhat[1]*ry + dspec.bhat[2]*rz;
          skernel[p] = SPHa * (3*br*br - 1) / (rmag*rmag*rmag);

          if (fabs(skernel[p]) > maxkernelvalue)
            maxkernelvalue = fabs(skernel[p]);
        }
      }
    }
  }
  // if (rank==0) printroot("   max value : %0.20e\n", maxkernelvalue);
  minkernelvalue = threshold*maxkernelvalue;
  // if (rank==0) printroot("   threshold value : %0.20e\n", minkernelvalue);
  for (int i = 0; i < N; i++) {
    if (fabs(skernel[i]) < minkernelvalue) {
      skernel[i] = 0;
      nnz--;
    }
  }
}

// returns the number of non-zeros
void Kernel::CreateCylindricalKernel(const DataSpec &dspec) {
  usedtype a = 0.475;
  usedtype CYL2alpha = 2*a;
  usedtype CYL3alpha = 3*a;
  usedtype CYL4alpha3 = 1/(4*a*a*a);
  usedtype CYLa = dspec.B0 / 2 / M_PI;
  usedtype cmag = sqrt(dspec.caxis[0]*dspec.caxis[0] +
                       dspec.caxis[1]*dspec.caxis[1] +
                       dspec.caxis[2]*dspec.caxis[2]);
  usedtype cx = dspec.caxis[0]/cmag;
  usedtype cy = dspec.caxis[1]/cmag;
  usedtype cz = dspec.caxis[2]/cmag;
  usedtype cosbeta = cx*dspec.bhat[0] +
      cy*dspec.bhat[1] +
      cz*dspec.bhat[2];
  usedtype CYLsin2beta = 1 - cosbeta*cosbeta;
  usedtype gx = dspec.bhat[0] - cosbeta * cx;
  usedtype gy = dspec.bhat[1] - cosbeta * cy;
  usedtype gz = dspec.bhat[2] - cosbeta * cz;
  usedtype minkernelvalue, maxkernelvalue = 0;

  for (int x = -halfsize; x <= halfsize; x++) {
    for (int y = -halfsize; y <= halfsize; y++) {
      for (int z = -halfsize; z <= halfsize; z++) {
        int p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
        if (x == 0 && y == 0 && z == 0) {
          ckernel[p] = dspec.B0/6 * (3 * cosbeta*cosbeta - 1);
        } else {
          usedtype rc = x*cx + y*cy + z*cz;
          if (fabs(rc) <= CYL2alpha) {
            usedtype r2 = 1/(x*x + y*y + z*z - rc*rc);
            usedtype rg = x*gx + y*gy + z*gz;

            usedtype h = CYL2alpha-fabs(rc);
            usedtype q = h*h * (CYL3alpha-h) * CYL4alpha3;

            ckernel[p] = CYLa * (2*rg*rg*r2 - CYLsin2beta) * r2 * q;
          } else {
            ckernel[p] = 0.0;
          }

          if (fabs(ckernel[p]) > maxkernelvalue)
            maxkernelvalue = fabs(ckernel[p]);
        }
      }
    }
  }

  minkernelvalue = threshold*dspec.B0;
  for (int i = 0; i < N; i++) {
    if (fabs(ckernel[i]) < minkernelvalue) {
      ckernel[i] = 0;
      nnz--;
    }
  }
}
void Kernel::InitMixedModel(const DataSpec &dspec) {
  usedtype        a = 0.475;
  CYL2alpha = 2*a;
  CYL3alpha = 3*a;
  CYL4alpha3 = 1/(4*a*a*a);
  CYLa = dspec.B0 / 2 / M_PI;

  for (int i = 0; i < modelmap.ncyls; i++) {
    usedtype cmag = sqrt(modelmap.x[i]*modelmap.x[i] +
                         modelmap.y[i]*modelmap.y[i] +
                         modelmap.z[i]*modelmap.z[i]);
    usedtype cx = modelmap.x[i]/cmag;
    usedtype cy = modelmap.y[i]/cmag;
    usedtype cz = modelmap.z[i]/cmag;

    modelmap.x[i] = cx;
    modelmap.y[i] = cy;
    modelmap.z[i] = cz;

    usedtype cosbeta = cx*dspec.bhat[0] + cy*dspec.bhat[1] + cz*dspec.bhat[2];

    sin2beta[i] = 1 - cosbeta*cosbeta;

    ctr[i] = dspec.B0/6 * (3 * cosbeta*cosbeta - 1);
    gx[i] = dspec.bhat[0] - cosbeta * cx;
    gy[i] = dspec.bhat[1] - cosbeta * cy;
    gz[i] = dspec.bhat[2] - cosbeta * cz;
  }
}

bool Kernel::Create(const models &model,
               const DataSpec &dspec,
               const usedtype &threshold) {
  // if (rank==0) printroot("Creating kernel\n");

  this->model = model;
  this->threshold = threshold;
  B0 = dspec.B0;

  if (model == MODEL_SPHERICAL) {
    size = ceil(pow(1.0/threshold, 1.0/3.0));
  } else {
  //    else if (model == MODEL_CYLINDRICAL) {
  //            size = ceil(pow(1.0*threshold, -1.0/2.0));
  //    }
    // model == MODEL_MIXED
    // Create a spherical kernel. Cylinders will be calculated on the fly.
    size = ceil(pow(1.0/threshold, 1.0/3.0));
  }

  // if (iseven(size)) size++;
  if (size/2 == size/2.0) size++;

  long maxdatasize = dspec.size[0] > dspec.size[1] &
      dspec.size[0] > dspec.size[2] ? dspec.size[0] :
      dspec.size[1] > dspec.size[2] ? dspec.size[1] : dspec.size[2];
  if (size > maxdatasize * 2 + 1) size = maxdatasize * 2 + 1;

  //  if (rank==0) printroot("   threshold = %0.3e\n", threshold);
  //  if (rank==0) printroot("   kernel size = %d\n", size);
  halfsize = size/2;
  N = size*size*size;
  nnz = N;
  if (model == MODEL_SPHERICAL) {
    skernel = reinterpret_cast<usedtype*>(calloc(nnz, sizeof(usedtype)));
  } else {
  //    else if (model == MODEL_CYLINDRICAL) {
  //            ckernel = (usedtype*) calloc(nnz, sizeof(usedtype));
  //    }
    skernel = reinterpret_cast<usedtype*>(calloc(nnz, sizeof(usedtype)));
    ctr = reinterpret_cast<usedtype*>(calloc(modelmap.ncyls, sizeof(usedtype)));
    sin2beta = reinterpret_cast<usedtype*>(
        calloc(modelmap.ncyls, sizeof(usedtype)));
    gx = reinterpret_cast<usedtype*>(calloc(modelmap.ncyls, sizeof(usedtype)));
    gy = reinterpret_cast<usedtype*>(calloc(modelmap.ncyls, sizeof(usedtype)));
    gz = reinterpret_cast<usedtype*>(calloc(modelmap.ncyls, sizeof(usedtype)));
  }

  yoffset = size;
  zoffset = size*size;

  // Create kernel
  if (model == MODEL_SPHERICAL) {
    CreateSphericalKernel(dspec);

  //    else if (model == MODEL_CYLINDRICAL)
  //            CreateCylindricalKernel(dspec);

  } else {
    CreateSphericalKernel(dspec);
    InitMixedModel(dspec);
  }

  // if (rank==0) printroot("   non-zeros = %d\n", nnz);

  return true;
}
usedtype Kernel::Get(int x, int y, int z, int o) {
  int ix = x + y*yoffset + z*zoffset;
  int mix;
  if (model == MODEL_MIXED &&  (mix = modelmap.mask[o]) != -1) {
    x -= halfsize;
    y -= halfsize;
    z -= halfsize;
    if (x == 0 && y == 0 && z == 0) {
      return ctr[mix];
    } else {
      usedtype rc = x*modelmap.x[mix] + y*modelmap.y[mix] + z*modelmap.z[mix];
      if (fabs(rc) <= CYL2alpha) {
        usedtype r2 = 1/(x*x + y*y + z*z - rc*rc);
        usedtype rg = x*gx[mix] + y*gy[mix] + z*gz[mix];

        usedtype h = CYL2alpha-fabs(rc);
        usedtype q = h*h * (CYL3alpha-h) * CYL4alpha3;
        usedtype retval = CYLa * (2*rg*rg*r2 - sin2beta[mix]) * r2 * q;
        if (fabs(retval) < threshold*B0)
          return 0.0;
        else
          return retval;
      } else {
        return 0.0;
      }
    }
  } else {
  //    else if (model == MODEL_CYLINDRICAL) {
  //            return ckernel[ix];
  //    }
    return skernel[ix];
  }
}

usedtype Kernel::GetCyl(int mix, int x, int y, int z) {
  usedtype rc = x*modelmap.x[mix] + y*modelmap.y[mix] + z*modelmap.z[mix];
  if (fabs(rc) <= CYL2alpha) {
    usedtype r2 = 1/(x*x + y*y + z*z - rc*rc);
    usedtype rg = x*gx[mix] + y*gy[mix] + z*gz[mix];

    usedtype h = CYL2alpha-fabs(rc);
    usedtype q = h*h * (CYL3alpha-h) * CYL4alpha3;

    usedtype retval = CYLa * (2*rg*rg*r2 - sin2beta[mix]) * r2 * q;
    if (fabs(retval) < threshold*B0)
      return 0.0;
    else
      return retval;
  } else {
    return 0.0;
  }
}

void Kernel::close(void) {
  modelmap.close();
  if (kernel != NULL) free(kernel);
  if (skernel != NULL) free(skernel);
  if (ckernel != NULL) free(ckernel);
  if (ctr != NULL) free(ctr);
  if (sin2beta != NULL) free(sin2beta);
  if (gx != NULL) free(gx);
  if (gy != NULL) free(gy);
  if (gz != NULL) free(gz);
}
