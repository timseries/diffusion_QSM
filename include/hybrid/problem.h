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

/*! \file problem.h
    \brief Problem class definitions file.

*/

#ifndef INCLUDE_PROBLEM_H_
#define INCLUDE_PROBLEM_H_

//#include "hybrid/basictypes.h"

#include "hybrid/basictypes.h"
#include "hybrid/dataspec.h"
#include "hybrid/kernel.h"

class Problem {
 public:
  Problem();
  Problem(Kernel &kernel, DataSpec &dspec, Real tau, Real alpha, Real beta);
  virtual ~Problem();
  DataSpec &dspec;
  int dStart, dEnd;
  int dN;  

  Real *Ax_b;  
  Real *AtAx_b;
  Real *Dx;
  Real *DtDx;
  Real *x;

  bool PreCalcCylinders;
  Real *cylColumns;
  int *FGindices; // indices corresponding to foreground elements

#ifdef USE_OPENCL
  OpenCL* cl;
  OpenCLProfile profile1, profile2, profile3;

  size_t cl_size_fg;
  size_t cl_size_n;
  size_t cl_size_cyl;

  //Kernels
  cl_kernel kernel_iterate1; 
  cl_kernel kernel_iterate2;
  cl_kernel kernel_zero;
  cl_kernel kernel_delta_b;
  cl_kernel kernel_rms_new_x;
  cl_kernel kernel_collect;

  //Buffers
  cl_mem cl_Ax_b, cl_x, cl_AtAx_b, cl_Dx, cl_DtDx;
  cl_mem cl_FGindices;
  cl_mem cl_skernel, cl_mask, cl_deltab;
  cl_mem cl_gx, cl_gy, cl_gz, cl_ctr, cl_sin2beta;
  cl_mem cl_mapx, cl_mapy, cl_mapz;
  cl_mem cl_out;

  int threads;
  float divide;
  char* clkern;

  OutputCL out; //output var buffer
#endif //USE_OPENCL
};
#endif  // INCLUDE_PROBLEM_H_
