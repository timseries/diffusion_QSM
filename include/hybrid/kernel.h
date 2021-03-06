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
// THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
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

/*! \file kernel.h
    \brief Kernel class definitions file.

*/

#ifndef INCLUDE_KERNEL_H_
#define INCLUDE_KERNEL_H_

#include "hybrid/common.h"
#include "hybrid/dataspec.h"
#include "hybrid/modelmap.h"

#ifdef USE_FOURIER_SPHERES
#include <fftw3.h>
#endif

/**
* Kernel class. Used to create, store, and compute the spherical and cylinderical kernel elements of the matrix \f$ \mathbf{A} \f$.
*/
class Kernel {
 public:
/**
* Kernel constructor.
*/
  Kernel();
/**
* Kernel destructor.
*/
  virtual ~Kernel();
/**
* Creates a 3D spherical kernel with a threshold-restricted support in the spatial domain (and Fourier domain if enabled) with maximum size the largest dimension of \f$ \mathbf{\Delta B}\f$ according to: 
*/
  void CreateSphericalKernel(const DataSpec &dspec);
/**
* Creates a 3D cylindrical kernel with a threshold-restricted support in the spatial domain.
*/
  void CreateCylindricalKernel(const DataSpec &dspec);
/**
* Creates a 3D cylindrical kernel with a threshold-restricted support in the spatial domain.
*/
  void InitMixedModel(const DataSpec &dspec);
/**
* Initialization method for data members of class. Computes kernel sizes (and halfwidths) and allocates/initializes their memory space.
* @param &model pointer to an object of enum models.
* @param &dspec pointer to an object of Dataspec.
* @param threshold pointer to threshold for determining kernel spatial extent (when the values are close enough to zero).
* @return True if successful, false otherwise.
*/
  bool Create(const models &model,
              const DataSpec &dspec,
              const Real &threshold);
/**
* Method to compute/get the cylindrical or spherical x,y,z kernel element based on voxel index o. 
* @param x The x component of the kernel.
* @param y The y component of the kernel.
* @param z The z component of the kernel.
* @param o The index of the foreground voxel.
* @return The kernel element.
*/
  Real Get(int x, int y, int z, int o);
/**
* Method to compute/get the cylindrical x,y,z kernel element based on the foreground index mix.
* @param x The x component of the kernel.
* @param y The y component of the kernel.
* @param z The z component of the kernel.
* @param mix The index of the foreground voxel.
* @return The kernel element.
*/
  Real GetCyl(int mix, int x, int y, int z);
/**
* Method to deallocate memory of class data members.
*/
  void close(void);

  models model;
  Real *skernel;
  Real *ckernel;
#ifdef USE_FOURIER_SPHERES
  Real *skernel_size_N; ///< Temporary array used to create the FFT of skernel.
  fftwf_complex *skernel_fft; ///< FFT of spherical kernel.
  fftwf_plan skernel_plan_forward; ///< FFTW plan to do forward transform of spherical kernel.
#endif
  Real *kernel;
  Real threshold;
  Real B0;
  Real CYL2alpha;

  Real CYL3alpha;
  Real CYL4alpha3;
  Real CYLa;
  Real *ctr;
  Real *sin2beta;
  Real *gx;
  Real *gy;
  Real *gz;
  int yoffset;
  int zoffset;
  ModelMap modelmap;
  int nnz;
  int _nnz;
  int size;
  int halfsize;
  int N;
};
#endif  // INCLUDE_KERNEL_H_
