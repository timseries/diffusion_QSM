// Copyright (c) 2013, Amanda Ng and Timothy Roberts
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
// THIS SOFTWARE IS PROVIDED BY Amanda Ng and Timothy Roberts ''AS IS'' AND ANY
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
// Author: amanda.ng@gmail.com, timothy.daniel.roberts@gmail.com

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

#ifdef USE_FOURIER_SPHERES
#include <fftw3.h>
#endif

Kernel::Kernel() {
  kernel = NULL;
  nnz = 0;
  _nnz = 0;
  skernel = NULL;
#ifdef USE_FOURIER_SPHERES
  skernel_fft=NULL;
#endif
  ckernel = NULL;
  ctr = NULL;
  sin2beta = NULL;
  gx = NULL;
  gy = NULL;
  gz = NULL;
}
Kernel::~Kernel() {}
void Kernel::CreateSphericalKernel(const DataSpec &dspec) {
  Real SPHa = dspec.B0 / (4 * M_PI);
  Real minkernelvalue, maxkernelvalue = 0;
  for (int x = -halfsize; x <= halfsize; x++) {
    for (int y = -halfsize; y <= halfsize; y++) {
      for (int z = -halfsize; z <= halfsize; z++) {
        int p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
        if (x == 0 && y == 0 && z == 0) {
          skernel[p] = 0;
          nnz--;
        } else {
          Real rmag = sqrt(x*x + y*y + z*z);
          Real rx = x/rmag;
          Real ry = y/rmag;
          Real rz = z/rmag;
          Real br = dspec.bhat[0]*rx + dspec.bhat[1]*ry + dspec.bhat[2]*rz;
          skernel[p] = SPHa * (3*br*br - 1) / (rmag*rmag*rmag);

          if (fabs(skernel[p]) > maxkernelvalue)
            maxkernelvalue = fabs(skernel[p]);
        }
      }
    }
  }
  minkernelvalue = threshold*maxkernelvalue;
  for (int i = 0; i < N; i++) {
    if (fabs(skernel[i]) < minkernelvalue) {
      skernel[i] = 0;
      nnz--;
    }
  }
  //take the FFT of the spherical kernel here. This may still be incorrect.
#ifdef USE_FOURIER_SPHERES
  int p = 0;
  int p_skernel = 0;
  //copy the skernel into a full N-sized array to take the appropriately-sized FFT. 
  //this copying is done in a circularly-shifted manner so as to avoid edge wrapping effects
  for (int x = -halfsize; x <= halfsize; x++) {
    for (int y = -halfsize; y <= halfsize; y++) {
      for (int z = -halfsize; z <= halfsize; z++) {
        p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
        //compute the wrapped index for skernel_size_N
        if ((x < 0) || (y < 0) || (z < 0)) {
          p_skernel = dspec.N;
          p_skernel -= abs(x + y*yoffset + z*zoffset);
        } else {
          p_skernel = 0;
          p_skernel += x + y*yoffset + z*zoffset;
        }
          // printroot("p_skernel: %d\n",p_skernel);
        skernel_size_N[p_skernel]=skernel[p];
    }
  }
}
fftwf_execute(skernel_plan_forward);
free(skernel_size_N);
#endif
}

// returns the number of non-zeros
void Kernel::CreateCylindricalKernel(const DataSpec &dspec) {
  Real a = 0.475;
  Real CYL2alpha = 2*a;
  Real CYL3alpha = 3*a;
  Real CYL4alpha3 = 1/(4*a*a*a);
  Real CYLa = dspec.B0 / 2 / M_PI;
  Real cmag = sqrt(dspec.caxis[0]*dspec.caxis[0] +
                       dspec.caxis[1]*dspec.caxis[1] +
                       dspec.caxis[2]*dspec.caxis[2]);
  Real cx = dspec.caxis[0]/cmag;
  Real cy = dspec.caxis[1]/cmag;
  Real cz = dspec.caxis[2]/cmag;
  Real cosbeta = cx*dspec.bhat[0] +
      cy*dspec.bhat[1] +
      cz*dspec.bhat[2];
  Real CYLsin2beta = 1 - cosbeta*cosbeta;
  Real gx = dspec.bhat[0] - cosbeta * cx;
  Real gy = dspec.bhat[1] - cosbeta * cy;
  Real gz = dspec.bhat[2] - cosbeta * cz;
  Real minkernelvalue, maxkernelvalue = 0;

  for (int x = -halfsize; x <= halfsize; x++) {
    for (int y = -halfsize; y <= halfsize; y++) {
      for (int z = -halfsize; z <= halfsize; z++) {
        int p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
        if (x == 0 && y == 0 && z == 0) {
          ckernel[p] = dspec.B0/6 * (3 * cosbeta*cosbeta - 1);
        } else {
          Real rc = x*cx + y*cy + z*cz;
          if (fabs(rc) <= CYL2alpha) {
            Real r2 = 1/(x*x + y*y + z*z - rc*rc);
            Real rg = x*gx + y*gy + z*gz;

            Real h = CYL2alpha-fabs(rc);
            Real q = h*h * (CYL3alpha-h) * CYL4alpha3;

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
  Real        a = 0.475;
  CYL2alpha = 2*a;
  CYL3alpha = 3*a;
  CYL4alpha3 = 1/(4*a*a*a);
  CYLa = dspec.B0 / 2 / M_PI;

  for (int i = 0; i < modelmap.ncyls; i++) {
    Real cmag = sqrt(modelmap.x[i]*modelmap.x[i] +
                         modelmap.y[i]*modelmap.y[i] +
                         modelmap.z[i]*modelmap.z[i]);
    Real cx = modelmap.x[i]/cmag;
    Real cy = modelmap.y[i]/cmag;
    Real cz = modelmap.z[i]/cmag;

    modelmap.x[i] = cx;
    modelmap.y[i] = cy;
    modelmap.z[i] = cz;

    Real cosbeta = cx*dspec.bhat[0] + cy*dspec.bhat[1] + cz*dspec.bhat[2];

    sin2beta[i] = 1 - cosbeta*cosbeta;

    ctr[i] = dspec.B0/6 * (3 * cosbeta*cosbeta - 1);
    gx[i] = dspec.bhat[0] - cosbeta * cx;
    gy[i] = dspec.bhat[1] - cosbeta * cy;
    gz[i] = dspec.bhat[2] - cosbeta * cz;
  }
}

bool Kernel::Create(const models &model,
               const DataSpec &dspec,
               const Real &threshold) {
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

  //ensure kernel width is odd to ensure correct centering
  if (size/2 == size/2.0) size++;

  //this is wrong, needs to be min datasize, consult with Amanda
  long maxdatasize = dspec.size[0] > dspec.size[1] &
      dspec.size[0] > dspec.size[2] ? dspec.size[0] :
      dspec.size[1] > dspec.size[2] ? dspec.size[1] : dspec.size[2];
  long mindatasize = dspec.size[0] < dspec.size[1] &
      dspec.size[0] < dspec.size[2] ? dspec.size[0] :
      dspec.size[1] < dspec.size[2] ? dspec.size[1] : dspec.size[2];
  if (size > maxdatasize * 2 + 1) size = maxdatasize * 2 + 1;
  // if (size > mindatasize * 2 + 1) size = mindatasize;
  // if (size/2 == size/2.0) size--;

  //  if (rank==0) printroot("   threshold = %0.3e\n", threshold);
  //  if (rank==0) printroot("   kernel size = %d\n", size);
  halfsize = size/2;
  N = size*size*size;
  nnz = _nnz = N;
  skernel = reinterpret_cast<Real*>(calloc(nnz, sizeof(Real)));
#ifdef USE_FOURIER_SPHERES
  //allocate k-space kernel skernel_fft and temporary input array skernel_size_N
  skernel_fft = fftwf_alloc_complex(dspec.N_fft);
  skernel_size_N = (Real*) calloc(dspec.N, sizeof(Real));
  //create fftw plan
  skernel_plan_forward = fftwf_plan_dft_r2c_3d(dspec.size[2], dspec.size[1], dspec.size[0],
                                              skernel_size_N, skernel_fft, FFTW_MEASURE);
#endif 
  if (model == MODEL_MIXED) {
    ctr = reinterpret_cast<Real*>(calloc(modelmap.ncyls, sizeof(Real)));
    sin2beta = reinterpret_cast<Real*>(
        calloc(modelmap.ncyls, sizeof(Real)));
    gx = reinterpret_cast<Real*>(calloc(modelmap.ncyls, sizeof(Real)));
    gy = reinterpret_cast<Real*>(calloc(modelmap.ncyls, sizeof(Real)));
    gz = reinterpret_cast<Real*>(calloc(modelmap.ncyls, sizeof(Real)));
  }
  yoffset = size;
  zoffset = size*size;

  // Create kernel
  if (model == MODEL_SPHERICAL) {
    CreateSphericalKernel(dspec);
  } else {
    CreateSphericalKernel(dspec);
    InitMixedModel(dspec);
  }
  return true;
}
Real Kernel::Get(int x, int y, int z, int o) {
  int ix = x + y*yoffset + z*zoffset;
  int mix;
  if (model == MODEL_MIXED && (mix = modelmap.mask[o]) != -1) {
    x -= halfsize;
    y -= halfsize;
    z -= halfsize;
    if (x == 0 && y == 0 && z == 0) {
      return ctr[mix];
    } else {
      Real rc = x*modelmap.x[mix] + y*modelmap.y[mix] + z*modelmap.z[mix];
      if (fabs(rc) <= CYL2alpha) {
        Real r2 = 1/(x*x + y*y + z*z - rc*rc);
        Real rg = x*gx[mix] + y*gy[mix] + z*gz[mix];

        Real h = CYL2alpha-fabs(rc);
        Real q = h*h * (CYL3alpha-h) * CYL4alpha3;
        Real retval = CYLa * (2*rg*rg*r2 - sin2beta[mix]) * r2 * q;
        if (fabs(retval) < threshold*B0)
          return 0.0;
        else
          return retval;
      } else {
        return 0.0;
      }
    }
  } else {
    return skernel[ix];
  }
}

Real Kernel::GetCyl(int mix, int x, int y, int z) {
  Real rc = x*modelmap.x[mix] + y*modelmap.y[mix] + z*modelmap.z[mix];
  if (fabs(rc) <= CYL2alpha) {
    Real r2 = 1/(x*x + y*y + z*z - rc*rc);
    Real rg = x*gx[mix] + y*gy[mix] + z*gz[mix];

    Real h = CYL2alpha-fabs(rc);
    Real q = h*h * (CYL3alpha-h) * CYL4alpha3;

    Real retval = CYLa * (2*rg*rg*r2 - sin2beta[mix]) * r2 * q;
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
