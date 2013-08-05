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

/*! \file problem.cc
  \brief Problem class file.

  Implementation of the Problem class.
*/

#include "hybrid/problem.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

Problem::Problem(Kernel &kernel, DataSpec &dspec, Real tau, Real alpha, Real beta) : dspec(dspec) {
  // Get start and end of local portion of Deltab array
    
  //==================================================================================================================
  // Create Ax_b array
  printroot("   Creating Ax_b array ...\n");
  Ax_b = (Real*) calloc(dspec.range, sizeof(Real));
    
  //==================================================================================================================
  // Create AtAx_b array
  printroot("   Creating AtAx_b array ...\n");
  AtAx_b = (Real*) calloc(dspec.nFG, sizeof(Real));
    
  //==================================================================================================================
  // Create Dx array
  printroot("   Creating Dx array ...\n");
  Dx = (Real*) calloc(dspec.range, sizeof(Real));
    
  //==================================================================================================================
  // Create DtDx array
  printroot("   Creating DtDx array ...\n");
  DtDx = (Real*) calloc(dspec.nFG, sizeof(Real));
    
  //==================================================================================================================
  // Create x array
  printroot("   Creating x array ...\n");
  x = (Real*) calloc(dspec.nFG, sizeof(Real));
    
  //==================================================================================================================
  // Create foreground indices array
  printroot("   Creating foreground indices array ...\n");
  FGindices = (int*) calloc(dspec.N, sizeof(int));

  //==================================================================================================================
  // Create cylColumns array
  double M = (double)kernel.modelmap.ncyls * dspec.range * sizeof(Real);
  cylColumns = NULL;
  //if (M < 8000000000 / size)  //Hack, hard coded to 8GB total mem limit for now, should check GPU mem avail
  if (1)  //Force cylinder calc on the fly for testing
  {
    printf("Require %.0f bytes for cylColumns array, attempting allocation...", M);
    cylColumns = (Real*) malloc(M);
  }
  if (cylColumns == NULL) {
    printroot("   Not enough memory to pre-calculate cylinder kernels.\n");
    printroot("   Cylinder kernels will be calculated on-the-fly, which will take longer to process.\n");
    PreCalcCylinders = false;
  }
  else {
    PreCalcCylinders = true;
  }
  
  /* Init OpenCL */
#ifdef USE_OPENCL
  //Create OpenCL context
  threads = 128;
  arghandler.GetArg("-threads", threads);
  cl = cl_new(dspec.N, 1, threads, rank); //Initial size set as full voxel count
  assert(cl);

  divide = 1.0;
  arghandler.GetArg("-divide", divide);

  char* defkern = "landweber.cl";
  clkern = defkern;
  arghandler.GetArg("-cl_kernel", clkern);

  //Pass in these values as constants defined in kernel source.
  //Results in much faster processing on pre-fermi nvidia cards
  // as opposed to passing in global variables to kernel, seems newer
  // compilers are better at automatically caching reads of constant global values?
  cl_define_integer_constant(cl, "nFG", dspec.nFG);
  cl_define_integer_constant(cl, "dsize0", dspec.size[0]);
  cl_define_integer_constant(cl, "dsize1", dspec.size[1]);
  cl_define_integer_constant(cl, "zoffset", dspec.zoffset);
  cl_define_integer_constant(cl, "yoffset", dspec.yoffset);
  cl_define_integer_constant(cl, "kyoffset", kernel.yoffset);
  cl_define_integer_constant(cl, "kzoffset", kernel.zoffset);
  cl_define_integer_constant(cl, "khalfsize", kernel.halfsize);
  cl_define_real_constant(cl, "C_10_96", 10.0/96.0);
  cl_define_real_constant(cl, "C_3_96", 3.0/96.0);
  cl_define_real_constant(cl, "CYL2alpha", kernel.CYL2alpha);
  cl_define_real_constant(cl, "CYL3alpha", kernel.CYL3alpha);
  cl_define_real_constant(cl, "CYL4alpha3", kernel.CYL4alpha3);
  cl_define_real_constant(cl, "CYLa", kernel.CYLa);
  cl_define_real_constant(cl, "thresholdB0", kernel.threshold*kernel.B0);
  cl_define_real_constant(cl, "alpha", alpha);
  cl_define_real_constant(cl, "beta", beta);
  cl_define_real_constant(cl, "tau", tau);
  //printf(cl->constants);

  //Read & build the sources
  const char* kernel_src = cl_readfile(clkern);
  const char* common_src = cl_readfile("common.h");
  //Constants must be included after data type definitions from common.h
  const char *srcptrs[]={common_src, cl->constants, kernel_src};
  //Build source from 5 files, no additional compiler options
  cl_build(cl, srcptrs, "", 3);
  //Print log on root proc
  if (rank == 0) fprintf(stderr, "%s\n", cl->log);

  //Ensure enough local memory available for chosen kernel
  //cl_checklocalmem(cl, cl->nthreads * sizeof(Real) * (nlaw ? 7 : 4));

  // Allocate GPU global memory buffers
  cl_size_n = sizeof(Real) * dspec.range;
  cl_size_fg = sizeof(Real) * dspec.nFG;
  cl_size_cyl = sizeof(Real) * kernel.modelmap.ncyls;
  cl_Ax_b = cl_new_buffer(cl, CL_MEM_READ_WRITE, cl_size_n);
  cl_Dx = cl_new_buffer(cl, CL_MEM_READ_WRITE, cl_size_n);
  cl_AtAx_b = cl_new_buffer(cl, CL_MEM_READ_WRITE, cl_size_fg);
  cl_DtDx = cl_new_buffer(cl, CL_MEM_READ_WRITE, cl_size_fg);
  cl_x = cl_new_buffer(cl, CL_MEM_READ_WRITE, cl_size_fg);
  cl_FGindices = cl_new_buffer(cl, CL_MEM_READ_ONLY, sizeof(int) * dspec.N);
  cl_gx = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_gy = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_gz = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_ctr = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_sin2beta = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_mapx = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_mapy = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_mapz = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_cyl);
  cl_skernel = cl_new_buffer(cl, CL_MEM_READ_ONLY, kernel._nnz * sizeof(Real));
  cl_mask = cl_new_buffer(cl, CL_MEM_READ_ONLY, sizeof(int) * dspec.N);
  cl_deltab = cl_new_buffer(cl, CL_MEM_READ_ONLY, cl_size_n);
  //Resize to get new group calc...
  cl_size(cl, dspec.nFG, 0, threads, rank);
  cl_out = cl_new_buffer(cl, CL_MEM_WRITE_ONLY, sizeof(OutputCL) * cl->ngroups);

  // Get handles to the kernels and setup input/output buffer args
  kernel_iterate1 = cl_get_kernel(cl, "iterate1", 14, 
                                  cl_Ax_b, cl_x, cl_Dx, cl_FGindices, cl_gx, cl_gy, 
                                  cl_gz, cl_ctr, cl_sin2beta, cl_mapx, cl_mapy, cl_mapz,
                                  cl_skernel, cl_mask);
  kernel_iterate2 = cl_get_kernel(cl, "iterate2", 15,
                                  cl_AtAx_b, cl_Ax_b, cl_DtDx, cl_Dx, cl_FGindices,
                                  cl_gx, cl_gy, cl_gz, cl_ctr, cl_sin2beta, cl_mapx,
                                  cl_mapy, cl_mapz, cl_skernel, cl_mask);
  kernel_zero = cl_get_kernel(cl, "zero", 4, cl_Ax_b, cl_Dx, cl_AtAx_b, cl_DtDx);
  kernel_delta_b = cl_get_kernel(cl, "delta_b", 2, cl_Ax_b, cl_deltab);
  kernel_rms_new_x = cl_get_kernel(cl, "rms_new_x", 4, cl_AtAx_b, cl_DtDx, cl_x, cl_out);
  kernel_collect = cl_get_kernel(cl, "collect", 1, cl_out);
#endif
}

Problem::~Problem() {
  free(Ax_b);
  free(AtAx_b);
  free(Dx);
  free(DtDx);
  free(x);
  free(FGindices);
  free(cylColumns);

#ifdef USE_OPENCL
  //Free memory
  clReleaseMemObject(cl_Ax_b);
  clReleaseMemObject(cl_Dx);
  clReleaseMemObject(cl_AtAx_b);
  clReleaseMemObject(cl_DtDx);
  clReleaseMemObject(cl_x);
  clReleaseMemObject(cl_FGindices);
  clReleaseMemObject(cl_gx);
  clReleaseMemObject(cl_gy);
  clReleaseMemObject(cl_gz);
  clReleaseMemObject(cl_ctr);
  clReleaseMemObject(cl_sin2beta);
  clReleaseMemObject(cl_mapx);
  clReleaseMemObject(cl_mapy);
  clReleaseMemObject(cl_mapz);
  clReleaseMemObject(cl_skernel);
  clReleaseMemObject(cl_mask);
  clReleaseMemObject(cl_deltab);
  clReleaseMemObject(cl_out);
  clReleaseKernel(kernel_iterate1);
  clReleaseKernel(kernel_iterate2);
  clReleaseKernel(kernel_zero);
  clReleaseKernel(kernel_delta_b);
  clReleaseKernel(kernel_rms_new_x);
  clReleaseKernel(kernel_collect);
  cl_delete(cl);
#endif
}
