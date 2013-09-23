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
#ifdef USE_OPENCL
#include "hybrid/opencl_base.h"
#endif
#include "hybrid/common.h"
#include "hybrid/dataspec.h"
#include "hybrid/kernel.h"

#ifdef USE_FOURIER_SPHERES
#include <fftw3.h>
#endif

/**
* Problem class. Used to store the temporary data members specific to solving the landweber iterations of the process class. 
*/
class Problem {
 public:
/**
* Problem constructor.
*/
  Problem();
/**
* Problem constructor.
* @param &kernel pointer to an object of Kernel.
* @param &dspec pointer to an object of DataSpec.
* @param &arghandler pointer to an object of ArgHandler.
* @param tau regularization parameter \f$ \tau \f$.
* @param alpha regualrization parameter \f$ \alpha \f$.
* @param beta regualrization parameter \f$ \beta \f$.
* @param rank Rank of this process.
*/
  Problem(Kernel &kernel, DataSpec &dspec, ArgHandler &arghandler, Real tau, Real alpha, Real beta, int rank);
/**
* Method to reallocate the elements of this problem class based on new ranges of \f$ \mathbf{\Delta B} \f$ for each process.
* @param &kernel pointer to an object of Kernel.
* @param &dspec pointer to an object of DataSpec.
*/
  void Reallocate(Kernel &kernel, DataSpec &dspec);
/**
* Method to allocate the indices of foreground spheres and cylinders into an array FGindicesUniform such that the ratio of cylinders/spheres in any contiguous portion of this array is relatively constant. This is used for load balancing the thread parallelization of the Mulatadd method in Process.
* @param mask pointer to the foreground index binary mask.
* @param rank rank of this process
* @param &kernel pointer to an object of Kernel.
* @param &dspec pointer to an object of DataSpec.
* @see Process::MultAdd
*/
  void UniformFGIndices(bool* mask, int rank, Kernel &kernel, DataSpec &dspec);
/**
* Problem destructor.
*/
  virtual ~Problem();
  DataSpec &dspec;

  Real *Ax_b; ///< The result of computing Ax-b.
  Real *AtAx_b; ///< The result of computing A'(Ax-b).
#ifdef USE_FOURIER_SPHERES
  Real *Ax_spheres; ///< Store N-sized version of x for computing Ax in Fourier domain for spheres.
  fftwf_complex *x_full_fft_out; ///< FFT of \f$ \mathbf{x} \f$.
  fftwf_plan x_fft_plan_forward; ///< FFTW plan for FFT of \f$ \mathbf{x} \f$.
  fftwf_plan x_fft_plan_inverse; ///< FFTW plan for IFFT of \f$ \mathbf{x} \f$.
  Real *AtAx_spheres; ///< Store N-sized version of Ax_b. Copy necessary to operate in the Fourier domain and the spatial domain.
  fftwf_complex *Ax_b_fft; ///< FFT of \f$ \mathbf{Ax-b} \f$.
  fftwf_plan Ax_b_fft_plan_forward; ///< FFTW plan for FFT of \f$\mathbf{\Detla B}\f$.
  fftwf_plan Ax_b_fft_plan_inverse; ///< FFTW plan for IFFT of \f$\mathbf{\Detla B}\f$.
#endif
  Real *Dx; ///< Matrix to store the discrete Laplacian kernel.
  Real *DtDx; ///< Matrix to store the inverse discrete Laplacian kernel.
  Real *x; ///< The current solution \f$ \mathbf{x} \f$ of \f$ \mathbf{y=Ax+b} \f$. 
  int *workmatrix; ///<cumulative work matrix over the volume.
  
  bool PreCalcCylinders; ///< Binary flag, true if memory allows for cylinder elements of kernel can be pre-computed. False otherwise.
  Real *cylColumns; ///< Array storing the columns of cylinder kernels in \f$ \mathbf{A} \f$.
  int *FGindices; ///< indices corresponding to foreground elements
  int *FGindicesUniform; ///< indices corresponding to foreground elements, reordered so that the distribution of indices of cylinders and spheres is as uniform as possible along the array
  int *FGindicesCyl; ///< indices corresponding to foreground elements, reordered so that the distribution of indices of cylinders and spheres is as uniform as possible along the array
  int *FGindicesSphere; ///< indices corresponding to foreground elements, reordered so that the distribution of indices of cylinders and spheres is as uniform as possible along the array

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
