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
// Author: timothy.daniel.roberts@gmail.com, amanda.ng@gmail.com

/*! \file process.cc
  \brief Process class file.

  Implementation of the Process class.
*/

#include "hybrid/process.h"
#ifdef USE_OPENCL
#include "hybrid/opencl_base.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

#ifdef HPM
extern "C" {
	#include "hpm.h"
}
#endif

#ifdef MPI_PROFILE
#include <mpt.h>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif


using ::util::checkVersion;
using ::util::checkEndianness;
using ::util::byteswap;

Process::Process() {
  rank = 0;
  size = 0;
  mask = NULL;
  deltab = NULL;
  chi = NULL;
}
Process::~Process() {}
bool Process::Init(int argc, char** args) {
  char *filepath;
  MPI_File fptr;
  HandleArgs(argc, args);
  StartMPI(argc, args);
  // Initialize some extra dataspec parameters
  dspec.caxis[0] = 1;
  dspec.caxis[1] = 0;
  dspec.caxis[2] = 0;
  // initilize the output object
  myout.Init(arghandler, rank, size);
  // Load the data, mask, and create the model
  if (!loadDeltaB()) goto exitnow;
  if (!loadMask()) goto exitnow;
  if (!kernel.modelmap.Create(dspec,arghandler)) goto exitnow;
  printroot("ncyls modelmap: %d",kernel.modelmap.ncyls);
  // Adjust the dataspec using orthogonal recursive bisection 
  printroot("rank: %d\n",rank);
  printroot("dspec start: %d\n", dspec.start);
  printroot("dspec end: %d\n", dspec.end);

  myout.DistrArray(deltab, dspec.range, 3, dspec.size, "deltab");
  myout.LocalArray(0, mask, 3, dspec.size, "mask");
  // initialize chi vector
  chi = (Real*) calloc(dspec.N, sizeof(Real));
  if (arghandler.GetArg("-chi", filepath)) {
    if (rank == 0) {
      MPI_File_open(MPI_COMM_SELF, filepath,
                    MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
      MPI_File_read(fptr, chi, dspec.N, MPI_Real, MPI_STATUS_IGNORE);
      MPI_File_close(&fptr);
    }
    MPI_Bcast(chi, dspec.N, MPI_Real, 0, MPI_COMM_WORLD);
    free(filepath);
    filepath = NULL;
  }
  else {
    memset(chi, 0, dspec.N*sizeof(Real));
  }
  // CREATE KERNEL
  if (!kernel.Create(model, dspec, threshold)) goto exitnow;
  return 1;
exitnow:
  return 0;
}

void Process::HandleArgs(int argc, char** args) {
  char *str = NULL;
  arghandler.Init(argc, args);
#ifdef DEBUG
  if (rank==0) printroot("   Argument handler initialised\n");
#endif
  // GET THRESHOLD
  if (!arghandler.GetArg("-threshold", threshold)) {
    threshold = 1e-6;
  }
  if (rank==0) printroot("   Threshold = %0.3e\n", threshold);
  if (!arghandler.GetArg("-model", str)){
    model = MODEL_MIXED;
  } else {
    switch (str[0]) {
      case 's':
        model = MODEL_SPHERICAL;
        break;
      case 'm':
        model = MODEL_MIXED;
        break;
      default:
        if (rank==0) printroot("Unrecognised model");
        break;
    }
  }
  if (rank==0) printroot("   Model: %s\n", model == MODEL_SPHERICAL ?
            "spherical" : "mixed");

}
void Process::StartMPI(int argc, char** args) {
  MPI_Init(&argc, &args);
 #ifdef HPM
   hpmInit();
//   hpmStart("main function");
 #endif
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0) printroot("\n------------------------------------------\n");
  if (rank==0) printroot("MPI Environment\n");
  if (rank==0) printroot("Number of processes: %d\n", size);
  if (rank==0) printroot("This rank: %d\n", rank);
}
bool Process::loadDeltaB() {
  // reads in DeltaB to a vector
  char *filepath;
  int err;
  MPI_File fptr;
  MPI_Status status;
  MPI_Offset offset, disp;
  bool flgByteSwap;
  double buf[11];
  double bmag;
  double *castbuf;

  if (rank==0) printroot("Loading DeltaB ...\n");
  if (!arghandler.GetArg("-DeltaB", filepath)) {
    if (rank==0) printroot("Data file not specified\n");
    return false;
  }
  if (rank==0) printroot("   file: %s\n", filepath);

  // open file for reading
  err = MPI_File_open(MPI_COMM_SELF, filepath,
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
  if (err) {
    if (rank==0) printroot("Could not open file");
    return false;
  }

  // check endianness
  checkEndianness(fptr, flgByteSwap);

  // check version
  checkVersion(fptr, flgByteSwap);

  // read in header
  if (rank==0) printroot("Reading-in MPI header...\n");
  MPI_File_read(fptr, buf, 11, MPI_DOUBLE, &status);

  // byte swap header if required
  if (rank==0) printroot("Byte-swapping header...\n");
  if (flgByteSwap) byteswap((char*)buf,11,sizeof(double));

  // create the dataspec here
  if (rank==0) printroot("Creating dataspec...\n");
  dspec.Create(buf,rank,size);
  //printall("[%d] start = %d end = %d\n", rank, dspec.start, dspec.end);
  if (rank==0) printroot("Allocating deltaB array portion...\n");
  deltab = (Real*) calloc(dspec.range, sizeof(Real));
  if (rank==0) printroot("Finding start of deltab...\n");
  MPI_File_seek(fptr, dspec.start*sizeof(double), MPI_SEEK_CUR);

  if (sizeof(Real) == sizeof(double)) {
    MPI_File_read(fptr, deltab, dspec.range, MPI_DOUBLE, &status);
    if (flgByteSwap) {
      byteswap((char*)deltab, dspec.range, sizeof(double));
    }
  } else {
    if (rank==0) printroot("Casting temporary buffer into double...\n");
    castbuf = (double*) calloc(dspec.range, sizeof(double));
    MPI_File_read(fptr, castbuf, dspec.range, MPI_DOUBLE, &status);
    if (flgByteSwap) {
      byteswap((char*)castbuf, dspec.range, sizeof(double));
    }
    if (rank==0) printroot("Copying temporary buffer into deltab...\n");
    for (int i = 0; i < dspec.range; i++)
      deltab[i] = (float) castbuf[i];
    free(castbuf);
  }

  // Close file pointer
  MPI_File_close(&fptr);

  // set offsets
  /* dspec.yoffset = dspec.size[0]; */
  /* dspec.zoffset = dspec.size[0] * dspec.size[1]; */

  if (rank==0) printroot("   size = %d %d %d\n", dspec.size[0],
            dspec.size[1], dspec.size[2]);
  if (rank==0) printroot("   number of elements = %d\n", dspec.N);

  free(filepath);

  return true;
}

bool Process::loadMask() {
  // Collective function loads mask onto each process

  int err;
  char *filepath;
  MPI_File fptr;
  double buf[4];
  bool flgByteSwap;
  MPI_Offset offset;
  if (rank==0) printroot("Loading mask ...\n");
  if (!arghandler.GetArg("-mask", filepath)) {
    if (rank==0) printroot("Mask file not specified\n");
    return false;
  }
  if (rank==0) printroot("   file: %s\n", filepath);
  mask = (bool*) calloc(dspec.N, sizeof(bool));
  // open file for reading
  if (rank == 0) {
    err = MPI_File_open(MPI_COMM_SELF, filepath,
                        MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
    if (err) {
      if (rank==0) printroot("Could not open mask file");
      return false;
    }

    // check endianness
    checkEndianness(fptr, flgByteSwap);

    // check version
    checkVersion(fptr, flgByteSwap);

    // read in header
    MPI_File_read(fptr, buf, 4, MPI_DOUBLE, MPI_STATUS_IGNORE);

    // byte swap header if required
    if (flgByteSwap) byteswap((char*)buf,4,sizeof(double));

    // Check parameters match
    if (dspec.N != buf[0]) {
      if (rank==0) printroot("Number of elements in Mask does not match DeltaB");
      return false;}
    if (dspec.size[0] != buf[1]) {
      if (rank==0) printroot("Size of first dimension of Mask does not match DeltaB");
      return false;}
    if (dspec.size[1] != buf[2]) {
      if (rank==0) printroot("Size of second dimension of Mask does not match DeltaB");
      return false;}
    if (dspec.size[2] != buf[3]) {
      if (rank==0) printroot("Size of third dimension of Mask does not match DeltaB");
      return false;}

    MPI_File_read(fptr, mask, dspec.N, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fptr);

  }
  MPI_Bcast(mask, dspec.N, MPI_CHAR, 0, MPI_COMM_WORLD);
  //TODO(timseries): do this somewhere else
  dspec.nBG = 0; dspec.nFG = 0;
  for (int i=0; i<dspec.N; i++) {
    if (mask[i])
      dspec.nFG++;
    else
      dspec.nBG++;
  }
  if (rank==0) printroot("   number of foreground voxels = %d\n", dspec.nFG);
  if (rank==0) printroot("   number of background voxels = %d\n", dspec.nBG);
  if (rank==0) printroot("   background percentage = %0.1f%%\n", 100.0*dspec.nBG/dspec.N);

  free(filepath);

  return true;
}

bool Process::FullPass() {
  // Landweber iteration
  // x1_{n+1} = x_{n} - alpha * A' * (A * x_{n} - b) - beta * D' * (D * x_{n})
  // where b     = deltab,
  //       x     = chi,
  //       alpha = 0.5 * tau,
  //       beta  = 0.5 * tau,
  //       tau   = 2/norm(A), and
  //       D     = represents a 3D laplacian filter

  // TODO(timseries): these constants should be specified somewhere else, not as magic nubmers
  Real tau = 0.15, alpha = 0.75, beta = 0.25;
  float Lfactors[3] = {-1, 0.10416667f, 0.03125f};
  Real *cylColumns;
  int *FGindices; // indices corresponding to foreground elements
  int *FGindicesUniform; // indices corresponding to foreground elements, reordered so that the distribution of indices of cylinders and spheres is as uniform as possible along the array
  int *FGindicesCyl; // indices corresponding to foreground elements, reordered so that the distribution of indices of cylinders and spheres is as uniform as possible along the array
  int *FGindicesSphere; // indices corresponding to foreground elements, reordered so that the distribution of indices of cylinders and spheres is as uniform as possible along the array

  Real *new_x;

  int dStart, dEnd;
  int dN;  
  int recvcounts, displs;

  //stopping criteria
  Real relative_threshold = 1e-6;
  Real absolute_threshold = 1e-16;
  int max_iters = 1000;
  double rms_x, rms_diff_x; //, old_rms_diff_x;
  int iteration = 0;

  int o, ox, oy, oz;
  int p, px, py, pz;
  int rx, ry, rz;
  int _rx, _ry, _rz;
  int mix;

  int nthreads,chunk,tid;

  // int roffset;

  MPI_File fptr;
  char *fname;

  time_t tsecs;
  time_t tmins;
  time_t tStart, tEnd, tsave;
  time_t tIterStart1, tIterStart2, tIterEnd1, tIterEnd2, tRms, tOverhead;

  bool flg;

  //  Problem *P;

  if (rank==0) printroot("Full pass ...\n");

  fname = (char*) calloc(256, sizeof(char));

  flg = arghandler.GetArg("-rth", relative_threshold);
  flg = arghandler.GetArg("-ath", absolute_threshold);
  flg = arghandler.GetArg("-maxiters", max_iters);

  if (rank==0) printroot("   Relative threshold = %0.3e\n", relative_threshold);
  if (rank==0) printroot("   Absolute threshold = %0.3e\n", absolute_threshold);
  if (rank==0) printroot("   Maximum iterations = %d\n", max_iters);

  flg = arghandler.GetArg("-tau", tau);
  flg = arghandler.GetArg("-alpha", alpha);
  if (flg)
    beta = 1 - alpha;

  if (rank==0) printroot("   tau = %0.3f\n", tau);
  if (rank==0) printroot("   alpha = %0.3f\n", alpha);
  if (rank==0) printroot("   beta = %0.3f\n", beta);

  //==================================================================================================================
  // Get start and end of local portion of Deltab array

  //==================================================================================================================
  // Create Arrays....
  P = new Problem(kernel, dspec, arghandler, tau, alpha, beta, rank);

  //==================================================================================================================
  // Create foreground indices array and initialise x
  FGindices=P->FGindices;
  FGindicesCyl=P->FGindicesCyl;
  FGindicesSphere=P->FGindicesSphere;
  if (rank==0) printroot("   Creating foreground indices array ...\n");
  //FGindices = (int*) calloc(dspec.N, sizeof(int));
  o = 0;
  int nFGCyls=0;
  int nFGSpheres=0;
  for (p = 0; p < dspec.N; p++) {
    if (mask[p]) {
      FGindices[o] = p;
      if (kernel.modelmap.mask[p]==-1) {
        FGindicesSphere[nFGSpheres]=p;
        nFGSpheres++;
      } else {
        FGindicesCyl[nFGCyls]=p;
        nFGCyls++;
      }
      o++;
    }
  }
  if (rank==0) printroot("Number of spheres: %d\n", nFGSpheres);
  if (rank==0) printroot("Number of cylinders: %d\n", nFGCyls);
  // Create uniform foreground indices array and initialise x
  double ratio=0;
  int large_count=0;
  int small_count=0;
  int target_large_count=0;
  int target_small_count=0;
  double target_ratio=0;
  int* larger_array;
  int* smaller_array;
  if (nFGSpheres>nFGCyls) {
    larger_array=FGindicesSphere;
    smaller_array=FGindicesCyl;
    target_large_count=nFGSpheres;
    target_small_count=nFGCyls;
  } else {
    larger_array=FGindicesCyl;
    smaller_array=FGindicesSphere;
    target_large_count=nFGCyls;
    target_small_count=nFGSpheres;
  }
  target_ratio=((double) target_large_count)/target_small_count;
  FGindicesUniform=P->FGindicesUniform;
  if (rank==0) printroot("   Creating uniform foreground indices array ...\n");
  o = 0;

  for (p = 0; p < dspec.N; p++) {
    if (mask[p]) {
      if ((ratio<target_ratio) && (large_count<target_large_count)) { //add some from the large set
        FGindicesUniform[o] = larger_array[large_count];                   
        large_count++;
      } else if ((ratio>target_ratio) && (small_count<target_small_count)){
        FGindicesUniform[o] = smaller_array[small_count];                   
        small_count++;
      } else if (large_count<target_large_count) {
        FGindicesUniform[o] = larger_array[large_count];                   
        large_count=large_count+(target_large_count-large_count>0);
      } else if (small_count<target_small_count){
        FGindicesUniform[o] = smaller_array[small_count];                   
        small_count=small_count+(target_small_count-small_count>0);
      } else {
        if (rank==0) printroot("error: index out of bounds, sum of small/large=%d, total FG number=%d\n", 
                               small_count+large_count,dspec.nFG);
      }
      ratio = ((double) large_count)/small_count;
      o++;
    }
  }

  //==================================================================================================================
  // Initialise x

  flg = arghandler.GetArg("-x", fname);

  if (flg) {
    int err = MPI_File_open(MPI_COMM_SELF, fname,
                            MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
    if (err) {
      if (rank==0) printroot("Cannot find file %s\n", fname);
      return false;
    }
 
    MPI_File_read(fptr, &iteration, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read(fptr, P->x, dspec.nFG, MPI_Real, MPI_STATUS_IGNORE);
    
    MPI_File_close(&fptr);

    iteration++;
  }
  else {
    o = 0;
    for (p = 0; p < dspec.N; p++) {
      if (mask[p]) {
        P->x[o] = chi[p];
        o++;
      }
    }
  }


  // if (rank==0) printroot("   Creating new_x array ...\n");
  // new_x = (Real*) calloc(dspec.nFG, sizeof(Real));
  // memset(new_x, 0, dspec.nFG * sizeof(Real));

  //==================================================================================================================
  // Create cylColumns array
  cylColumns = P->cylColumns;
  if (cylColumns == NULL) {
    printroot("   Not enough memory to pre-calculate cylinder kernels.\n");
    printroot("   Cylinder kernels will be calculated on-the-fly, which will take longer to process.\n");
    P->PreCalcCylinders = false;
  }
  else {
    
    printroot("   Creating cylColumns array ...\n");
    
    // Initialise cylColumns array
    for (p = dspec.start; p < dspec.end; p++) {
      
      pz = p / dspec.zoffset;
      py = (p - pz * dspec.zoffset) / dspec.yoffset;
      px = p - py * dspec.yoffset - pz * dspec.zoffset;
      
      for (o = 0; o < dspec.nFG; o++) {
        //        mix = kernel.modelmap.mask[FGindices[o]];
               mix = kernel.modelmap.mask[FGindicesUniform[o]];
        if (mix != -1) {
          // oz = FGindices[o] / dspec.zoffset;
          // oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
          // ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;
          oz = FGindicesUniform[o] / dspec.zoffset;
          oy = (FGindicesUniform[o] - oz * dspec.zoffset) / dspec.yoffset;
          ox = FGindicesUniform[o] - oy * dspec.yoffset - oz * dspec.zoffset;
          
          rx = px - ox + kernel.halfsize;
          ry = py - oy + kernel.halfsize;
          rz = pz - oz + kernel.halfsize;          
          
          // cylColumns[mix * dspec.range + p - dspec.start] = kernel.Get(rx, ry, rz, FGindices[o]);
          cylColumns[mix * dspec.range + p - dspec.start] = kernel.Get(rx, ry, rz, FGindicesUniform[o]);
        }
      }
    }
    
    P->PreCalcCylinders = true;
  }

#ifdef USE_OPENCL
//Write buffers
  cl_enqueue_write(P->cl, P->cl_FGindices, sizeof(int) * P->dspec.N, P->FGindices);
  cl_enqueue_write(P->cl, P->cl_gx, P->cl_size_cyl, kernel.gx);
  cl_enqueue_write(P->cl, P->cl_gy, P->cl_size_cyl, kernel.gy);
  cl_enqueue_write(P->cl, P->cl_gz, P->cl_size_cyl, kernel.gz);
  cl_enqueue_write(P->cl, P->cl_ctr, P->cl_size_cyl, kernel.ctr);
  cl_enqueue_write(P->cl, P->cl_sin2beta, P->cl_size_cyl, kernel.sin2beta);
  cl_enqueue_write(P->cl, P->cl_mapx, P->cl_size_cyl, kernel.modelmap.x);
  cl_enqueue_write(P->cl, P->cl_mapy, P->cl_size_cyl, kernel.modelmap.y);
  cl_enqueue_write(P->cl, P->cl_mapz, P->cl_size_cyl, kernel.modelmap.z);
  cl_enqueue_write(P->cl, P->cl_skernel, kernel._nnz * sizeof(Real), kernel.skernel);
  cl_enqueue_write(P->cl, P->cl_mask, sizeof(int) * P->dspec.N, kernel.modelmap.mask);
  cl_enqueue_write(P->cl, P->cl_deltab, P->cl_size_n, deltab);
  // Write initial x buffer
  cl_enqueue_write(P->cl, P->cl_x, P->cl_size_fg, P->x);
#endif
  //==================================================================================================================
  // Start iterations
#ifdef HPM
  hpmStart("iteration");
#endif

  if (rank==0) printroot("   Starting iterations ...\n");
  //old_rms_diff_x = -1;

  tStart = MPI_Wtime();
  tsave = tStart + 3600;
  int dn = -1; //kernel.halfsize - 1;
  int dz = 0; //kernel.halfsize;
  int dp = 1; //kernel.halfsize + 1;

  double wall_time = MPI_Wtime();

  do {

    tIterStart1 = MPI_Wtime();
    if (rank==0) printroot("      Iteration %d: ", iteration);

    // printroot("rank: %d first multadd start :%.3f\n", rank, MPI_Wtime());
    // Ax_b = A * x - b
    // Dx = D * x
    MultAdd(P->Ax_b,P->Dx,P->x,P->x,deltab,true,iteration);
    // printroot("rank: %d first multadd finish:%.3f\n", rank, MPI_Wtime());
    tIterEnd1 = MPI_Wtime();
    tsecs = tIterEnd1 - tIterStart1;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    printroot("first multadd iteration %d, rank:%d (%ldmin %lds) \n", iteration, rank, tmins, tsecs);
    tIterStart2 = MPI_Wtime();

    // AtAx_b = A' * Ax_b
    // DtDx = D' * Dx
    // printroot("rank: %d second multadd start: %.3f\n", rank, MPI_Wtime());
    MultAdd(P->AtAx_b,P->DtDx,P->Ax_b,P->Dx,NULL,false,iteration);
    // printroot("rank: %d second multadd finish: %.3f\n", rank, MPI_Wtime());
          // printroot("after second  multadd....\n",p);
          // printroot("P->Ax_b[8]: %0.3e\n",P->AtAx_b[8]);


    tIterEnd2 = MPI_Wtime();
    tsecs = tIterEnd2 - tIterStart2;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    // if (rank==0) printroot(" (%ldmin %lds)", tmins, tsecs);
    printroot("second multadd iteration %d, rank:%d (%ldmin %lds) \n", iteration, rank, tmins, tsecs);
    // Reduce AtAx_b
    tOverhead = MPI_Wtime();

    //if (rank==0) printroot("      reducing AtAx_b ...\n");
    // printroot("rank: %d first allreduce start: %.3f\n", rank, MPI_Wtime());
    tIterStart2 = MPI_Wtime();
    MPI_Allreduce(MPI_IN_PLACE, P->AtAx_b, P->dspec.nFG, MPI_Real, MPI_SUM, MPI_COMM_WORLD);
    tIterEnd2 = MPI_Wtime();
    tsecs = tIterEnd2 - tIterStart2;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    // if (rank==0) printroot(" (%ldmin %lds)", tmins, tsecs);
    printroot("first allreduce iteration %d, rank:%d (%ldmin %lds) \n", iteration, rank, tmins, tsecs);
    // printroot("rank: %d first allreduce end: %.3f\n", rank, MPI_Wtime());
    // printroot("rank: %d second allreduce start: %.3f\n", rank, MPI_Wtime());
    tIterStart2 = MPI_Wtime();
    MPI_Allreduce(MPI_IN_PLACE, P->DtDx, P->dspec.nFG, MPI_Real, MPI_SUM, MPI_COMM_WORLD);
    tIterEnd2 = MPI_Wtime();
    tsecs = tIterEnd2 - tIterStart2;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    // if (rank==0) printroot(" (%ldmin %lds)", tmins, tsecs);
    printroot("second allreduce iteration %d, rank:%d (%ldmin %lds) \n", iteration, rank, tmins, tsecs);
    // printroot("rank: %d second allreduce end: %.3f\n", rank, MPI_Wtime());
    tsecs = MPI_Wtime() - tOverhead;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    // if (rank==0) printroot(" (%ldmin %lds)", tmins, tsecs);
    // Calculate new x and rms values
    tRms = MPI_Wtime();
#ifdef USE_OPENCL
    //Calculate new_x and rms
    P->profile3.kern_time = 0;
    //Re-Write reduced AtAx_b and DtDx buffers
    cl_enqueue_write(P->cl, P->cl_AtAx_b, P->cl_size_fg, P->AtAx_b);
    cl_enqueue_write(P->cl, P->cl_DtDx, P->cl_size_fg, P->DtDx);
    // Perform the operations
    cl_enqueue_kernel(P->cl, P->kernel_rms_new_x, &P->profile3.event);   //Profile this kernel
    cl_enqueue_kernel(P->cl, P->kernel_collect, NULL);
    // Read the results back
    cl_enqueue_read(P->cl, P->cl_out, sizeof(OutputCL), &P->out);
    //Run queued operations
    cl_run(P->cl);
    //Profile
    cl_profile(P->cl, &P->profile3);
    rms_x = P->out.rms_x;
    rms_diff_x = P->out.rms_diff_x;
#else
    // Calculate new x and rms values
    rms_x = 0;
    rms_diff_x = 0;

    for (o = 0; o < dspec.nFG; o++) {
      Real new_x = P->x[o] - alpha * tau * P->AtAx_b[o] - beta * tau * P->DtDx[o];
      rms_x += new_x * new_x;
      rms_diff_x += (new_x - P->x[o]) * (new_x - P->x[o]);
      P->x[o] = new_x;  //Copy result to x[]
    }
#endif

    rms_x = sqrt(rms_x / dspec.nFG);
    rms_diff_x = sqrt(rms_diff_x / dspec.nFG);

    tIterEnd2 = MPI_Wtime();
    tsecs = tIterEnd2 - tIterStart1;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;

    if (rank==0) {
      printroot(" rms_diff_x = %0.3e rms_x = %0.3e rms_diff_x / rms_x  = %0.3e (%ldmin %lds)\n",
                rms_diff_x, rms_x, rms_diff_x / rms_x, tmins, tsecs);
    }
    // periodic saves just in case the job gets interrupted for any reason...
    tEnd = MPI_Wtime();
    if (rank == 0 && (tEnd > tsave || iteration < 2 || iteration == max_iters - 1)) {
      MPI_File fptr;
      sprintf(fname, "%s/x_iter%06d.bin", myout.outdir, iteration);
      MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
      MPI_File_set_size(fptr, 0);
      MPI_File_write(fptr, &iteration, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fptr, P->x, dspec.nFG, MPI_Real, MPI_STATUS_IGNORE);
      MPI_File_close(&fptr);

      tsave = tEnd + 3600;
    }

    iteration++;}
    while (rms_diff_x / rms_x > relative_threshold && rms_diff_x > absolute_threshold && iteration < max_iters);
#ifdef HPM
  hpmStop("iteration");
#endif

  tEnd = MPI_Wtime();
  tsecs = tEnd - tStart;
  tmins = floor(tsecs/60);

  if (rms_diff_x / rms_x <= relative_threshold)
    if (rank==0) printroot("   Reached relative threshold in %ldmin %lds\n", tmins, tsecs - tmins*60);
  else if (rms_diff_x <= absolute_threshold)
    if (rank==0) printroot("   Reached absolute threshold in %ldmin %lds\n", tmins, tsecs - tmins*60);
  else
    if (rank==0) printroot("   Reached maximum iterations in %ldmin %lds\n", tmins, tsecs - tmins*60);

  if (rank==0) printroot("   Completed solve.\n ");

  // Create chi vector
  if (rank==0) printroot("Compiling chi vector ...\n");
  memset(chi, 0, dspec.N*sizeof(Real));

  for (o = 0; o < dspec.nFG; o++) {
    chi[FGindicesUniform[o]] = P->x[o];
  }

    wall_time = MPI_Wtime() - wall_time;
#ifdef USE_OPENCL
  if (rank == 0)
  {
    printf("\n#procs, threads, per GPU size, global size, problem size, kernel1 time, kernel2 time, kernel3 time, wall time\n");
    printf("%d, %lu, %lu, %d, %d, ", size, P->cl->nthreads, P->cl->global_size, 0, P->cl->problem_size);
    printf("%5.2lf, %5.2lf, %5.2lf, %5.2lf\n", P->profile1.total_time, P->profile2.total_time, P->profile3.total_time, wall_time);
  }
#else
  if (rank == 0)
    fprintf(stderr, "\nProfiling: Wall time %5.2lf seconds\n", wall_time);
#endif
  // Free memory
  // if (rank==0) printroot("   Finished. Cleaning up ...\n");
  // free(x);
  // free(FGindices);
  // free(Ax_b);
  // free(AtAx_b);
  // free(new_x);
  // free(cylColumns);
  // free(Dx);
  // free(DtDx);
  // free(fname);

  return true;
}
//returns an array result which is: multiplier*multipicand+addend. Dir is an indexing option (true->forward,false->backward) 
// Eg.
// Forwards: Ax_b=Ax-b
// Backwards: AtAx_b=A'(Ax-b)
// Calls Mult version of this method, which just does  a multiply when no addend is specified.
  // void Process::MultAdd(Real result_fidelity, Real result_regularizer,
  //                     Real multiplicand_fidelity,Real multiplicand_regularizer, 
  //                     Real addend, bool dir) {
//void Process::MultAdd(Problem *P,bool dir) {
void Process::MultAdd(Real* result_fidelity, Real* result_regularizer, Real* multiplicand_fidelity, Real* multiplicand_regularizer, Real* addend, bool dir, int iteration) {

  float Lfactors[3] = {-1, 0.10416667f, 0.03125f};
  int recvcounts, displs;

  //stopping criteria
  Real relative_threshold = 1e-6;
  Real absolute_threshold = 1e-16;
  int max_iters = 1000;
  double rms_x, rms_diff_x; //, old_rms_diff_x;
  //int iteration = 0;

  int o, ox, oy, oz;
  int p, px, py, pz;
  int rx = 0;
  int ry = 0;
  int rz = 0;
  int _rx = 0;
  int _ry = 0;
  int _rz = 0;
  int mix = 0;

  int nthreads,chunk,tid;

#ifdef USE_OPENCL
  //TODO(timseries): reduce codesize here
  if (dir){
    P->profile1.kern_time = 0;
    cl_size(P->cl, P->dspec.range, 0, P->threads, rank);  //Resize
    cl_set_arg(cl, P->kernel_iterate1, 14, P->dspec.start);
    cl_set_arg(cl, P->kernel_iterate1, 15, P->dspec.end);
    // Perform the operations
    cl_set_arg(cl, P->kernel_zero, 4, P->dspec.range);
    cl_enqueue_kernel(P->cl, P->kernel_zero, NULL);
    //Only have x values on 2nd and subsequent iterations
    if (iteration > 0)
    {
      //Split the job up if requested
      int BLOCK = (dspec.nFG / P->divide + 0.5);
      //printf("iterate1 0 ");
      for (int start=0; start<dspec.nFG; start += BLOCK)
      {
        int end = start + BLOCK;
        if (end > dspec.nFG) end = dspec.nFG;
        //printf("%d ", end);
        cl_set_arg(cl, P->kernel_iterate1, 16, start);
        cl_set_arg(cl, P->kernel_iterate1, 17, end);
        cl_enqueue_kernel(P->cl, P->kernel_iterate1, &P->profile1.event);   //Profile this kernel
      }
      //printf("\n");
    }
    cl_set_arg(cl, P->kernel_delta_b, 2, P->dspec.start);
    cl_set_arg(cl, P->kernel_delta_b, 3, P->dspec.end);
    cl_enqueue_kernel(P->cl, P->kernel_delta_b, NULL);
    cl_run(P->cl);
    //Profile
    if (iteration > 0) cl_profile(P->cl, &P->profile1);
  } else {
        P->profile2.kern_time = 0;
    cl_size(P->cl, P->dspec.nFG, 0, P->threads, rank);  //Resize
    //Split the job up if requested
    int BLOCK = (P->dspec.range / P->divide + 0.5);
    //printf("iterate2 0 ");
    for (int start=P->dspec.start; start<P->dspec.end; start += BLOCK)
    {
      int end = start + BLOCK;
      if (end > P->dspec.end) end = P->dspec.end;
      //printf("%d ", end);
      cl_set_arg(cl, P->kernel_iterate2, 15, start);
      cl_set_arg(cl, P->kernel_iterate2, 16, end);
      cl_set_arg(cl, P->kernel_iterate2, 17, P->dspec.start);
      // Perform the operations
      cl_enqueue_kernel(P->cl, P->kernel_iterate2, &P->profile2.event);   //Profile this kernel
    }
    //printf("\n");
    // Read the results back
    cl_enqueue_read(P->cl, P->cl_AtAx_b, P->cl_size_fg, P->AtAx_b);
    cl_enqueue_read(P->cl, P->cl_DtDx, P->cl_size_fg, P->DtDx);
    //Run queued operations
    cl_run(P->cl);
    //Profile
    cl_profile(P->cl, &P->profile2);
  }
#else //assume we're using CPU on a bluegene or PC
    //TODO(timseries): fix this initialization section to not use P!
    if (dir) {
        memset(P->Ax_b, 0, (P->dspec.range) * sizeof(Real));
        memset(P->Dx, 0, (P->dspec.range) * sizeof(Real));
      }else{
        memset(P->AtAx_b, 0, (P->dspec.nFG) * sizeof(Real));
        memset(P->DtDx, 0, (P->dspec.nFG) * sizeof(Real));
    }
    int index1=0;
    int index2=0;
    int numcyls=0;
    int numspheres=0;
#ifdef USE_OPENMP
    if (rank==0 && iteration==0) printroot("Num of CPU: %d\n", omp_get_num_procs());
    if (rank==0 && iteration==0) printroot("Max threads: %d\n", omp_get_max_threads());
    chunk=P->dspec.nFG/omp_get_max_threads();
    if (rank==0 && iteration==0) printroot("Chunk size: %d\n", chunk);
#endif
#pragma omp parallel shared(nthreads,chunk,dir) private(tid,o,p,ox,oy,oz,px,py,pz,rx,ry,rz,_rx,_ry,_rz,mix,index1,index2,numcyls,numspheres) if (OPENMP)
    {
#ifdef USE_OPENMP
      tid = omp_get_thread_num();

      //      if (tid == 0 && rank==0 && iteration==0)
//      {
        nthreads = omp_get_num_threads();
        // printf("Starting  A*x-b using %d threads\n",nthreads);
        // printf("Initializing matrices on thread %d: \n", tid);
	//      }
        numspheres=0;
        numcyls=0;
#endif
#pragma omp for schedule(static, chunk)
    for (o = 0; o < P->dspec.nFG; o++) {
      mix = kernel.modelmap.mask[P->FGindicesUniform[o]];
      if (mix==-1){
        numspheres++;
      } else {
        numcyls++;
      }
      if ((P->x[o] and dir) or (!dir)) {
        
        // oz = P->FGindices[o] / P->dspec.zoffset;
        // oy = (P->FGindices[o] - oz * P->dspec.zoffset) / P->dspec.yoffset;
        // ox = P->FGindices[o] - oy * P->dspec.yoffset - oz * P->dspec.zoffset;
        oz = P->FGindicesUniform[o] / P->dspec.zoffset;
        oy = (P->FGindicesUniform[o] - oz * P->dspec.zoffset) / P->dspec.yoffset;
        ox = P->FGindicesUniform[o] - oy * P->dspec.yoffset - oz * P->dspec.zoffset;
        pz = P->dspec.start / P->dspec.zoffset;
        py = (P->dspec.start - pz * dspec.zoffset) / P->dspec.yoffset;
        px = P->dspec.start - py * P->dspec.yoffset - pz * P->dspec.zoffset - 1;			    

        rx = px - ox;
        ry = py - oy;
        rz = pz - oz;

        _rx = abs(rx);
        _ry = abs(ry);
        _rz = abs(rz);
//iterate through the elements of this part of the data...convolution
        for (p = P->dspec.start; p < P->dspec.end; p++) {
          index1 = dir ? o : p - P->dspec.start;
          index2 = dir ? p - P->dspec.start : o;
          px++;
          if (px == P->dspec.size[0]) {
            px = 0;
            py++;
            if (py == P->dspec.size[1]) {
              py = 0;
              pz++;
              rz = pz - oz;
              _rz = abs(rz);
            }
            ry = py - oy;
            _ry = abs(ry);
          }        
          rx = px - ox;
          _rx = abs(rx);
          
          if (_rx <= kernel.halfsize && _ry <= kernel.halfsize && _rz <= kernel.halfsize) {
            // Linear system
            // mix = kernel.modelmap.mask[P->FGindices[o]];
            mix = kernel.modelmap.mask[P->FGindicesUniform[o]];
            if (mix == -1) { // spherical kernel
              result_fidelity[index2] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * multiplicand_fidelity[index1];
            }
            else if (P->PreCalcCylinders) {
              result_fidelity[index2] += P->cylColumns[mix*P->dspec.range + p-P->dspec.start] * multiplicand_fidelity[index1];
              // if (dir) {
              //   P->Ax_b[index2] += P->cylColumns[mix*P->dspec.range + index2] * P->x[index1];
              // } else {
              //   P->AtAx_b[index2] += P->cylColumns[mix*P->dspec.range + index2] * P->Ax_b[index1];
              // }
            }
            else if (rx == 0 && ry == 0 && rz == 0) {
              result_fidelity[index2] += kernel.ctr[mix] * multiplicand_fidelity[index1];
            }
            else {
              result_fidelity[index2] += kernel.GetCyl(mix, rx, ry, rz) * multiplicand_fidelity[index1];
            }
            
            //OK 7/13: Optimised Laplacian,
            //This is not faster on CPU (as it is on GPU due to reduced branching)
            //but it is cleaner, using pre-calculated single precision constants also slightly faster
            if (_rx <= 1 && _ry <= 1 && _rz <= 1)
            {
              result_regularizer[index2] += Lfactors[_rx + _ry + _rz] *  multiplicand_regularizer[index1];
            }
          }
        }
      }
    }
    if (dir){
#pragma omp for
      for (p = 0; p < P->dspec.range; p++) {
        result_fidelity[p] -= addend[p];
      }
    }//end if (dir)
    if (rank==0) printroot("tid:%d numcyls:%d numspheres:%d\n", tid, numcyls, numspheres);
    numcyls=0;
    numspheres=0;
    }//end omp parallel section
#endif

}
bool Process::WriteOut(){
  char *filepath;
  MPI_File fptr;  
  myout.LocalArray(0, chi, 3, dspec.size, "chi");
  if (rank == 0) {
    filepath = (char*) calloc(strlen(myout.outdir) + 9, sizeof(char));
    sprintf(filepath, "%s/chi.bin", myout.outdir);
    MPI_File_open(MPI_COMM_SELF, filepath,
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
    MPI_File_write(fptr, chi, dspec.N, MPI_Real, MPI_STATUS_IGNORE);
    MPI_File_close(&fptr);
    free(filepath);
    filepath = NULL;
  }
  // TODO(timseries): include some error checking code as with the other methods in this class. Just return 1 for now....
  return 1;
}

bool Process::CleanUp(){
  if (rank==0) printroot("\n------------------------------------------\n");
  if (rank==0) printroot("CLEANING UP\n");
  // if (rank > 0) {
  //   MPI_Recv(NULL, 0, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // }
  if (rank==0) printroot("Freeing memory ...\n");
  if (deltab != NULL)			free(deltab);
  if (rank==0) printroot("Freed deltab ...\n");
  if (chi != NULL)			free(chi);
  if (rank==0) printroot("Freed chi ...\n");
  if (mask != NULL)			free(mask);
  if (rank==0) printroot("Freed mask ...\n");
  kernel.close();
  if (rank==0) printroot("Finished\n", rank); fflush(stdout);
  myout.Close();
  printf("finished rank=%d\n",rank);
  // if (rank < size -1) {
  //   MPI_Send(NULL, 0, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD);
  //     }
  MPI_Finalize();
#ifdef HPM
  // hpmStop("main function");
  hpmTerminate();
#endif
  // TODO(timhseries): include some error checking code as with the other methods in this class. Just return 1 for now....
  return 1;  
}
