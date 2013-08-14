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

//#include "hybrid/basictypes.h"

#include "hybrid/common.h"
#include "hybrid/dataspec.h"
#include "hybrid/modelmap.h"

class Kernel {
 public:
  Kernel();
  virtual ~Kernel();
  void CreateSphericalKernel(const DataSpec &dspec);
  void CreateCylindricalKernel(const DataSpec &dspec);
  void InitMixedModel(const DataSpec &dspec);
  bool Create(const models &model,
              const DataSpec &dspec,
              const Real &threshold);
  Real Get(int x, int y, int z, int o);
  Real GetCyl(int mix, int x, int y, int z);
  void close(void);

  models model;
  Real *skernel, *ckernel;
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
