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

/*! \file process.cc
  \brief Process class file.

  Implementation of the Process class.
*/

#include "hybrid/process.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

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
  filepath = NULL;
}
Process::~Process() {}
bool Process::Init(int argc, char** args) {
  HandleArgs(argc, args);
  StartMPI(argc, args);
  // Initialize some extra dataspec parameters
  dspec.caxis[0] = 1;
  dspec.caxis[1] = 0;
  dspec.caxis[2] = 0;
  // initilize the output object
  myout.Init(arghandler, rank, size);
  // Load the data, mask, and model
  if (!loadDeltaB()) goto exitnow;
  myout.DistrArray(deltab, dspec.end - dspec.start, 3, dspec.size, "deltab");
  if (!loadMask()) goto exitnow;
  myout.LocalArray(0, mask, 3, dspec.size, "mask");
  if (!kernel.modelmap.Process(dspec,arghandler)) goto exitnow;
  // initialize chi vector
  chi = (usedtype*) calloc(dspec.N, sizeof(usedtype));
  if (arghandler.GetArg("-chi", filepath)) {
    if (rank == 0) {
      MPI_File_open(MPI_COMM_SELF, filepath,
                    MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
      MPI_File_read(fptr, chi, dspec.N, MPI_USEDTYPE, MPI_STATUS_IGNORE);
      MPI_File_close(&fptr);
    }
    MPI_Bcast(chi, dspec.N, MPI_USEDTYPE, 0, MPI_COMM_WORLD);
    free(filepath);
    filepath = NULL;
  }
  else {
    memset(chi, 0, dspec.N*sizeof(usedtype));
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
  int ElsPerProc;

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
  MPI_File_read(fptr, buf, 11, MPI_DOUBLE, &status);

  // byte swap header if required
  if (flgByteSwap) byteswap((char*)buf,11,sizeof(double));

  // process header into relevant variables
  dspec.N = (int) buf[0];
  dspec.size[0] = (int) buf[1];
  dspec.size[1] = (int) buf[2];
  dspec.size[2] = (int) buf[3];
  dspec.B0 = buf[4];
  memcpy(dspec.bhat, &buf[5], 3*sizeof(double));

  // normalise bhat
  bmag = sqrt(dspec.bhat[0]*dspec.bhat[0] +
              dspec.bhat[1]*dspec.bhat[1] + dspec.bhat[2]*dspec.bhat[2]);
  dspec.bhat[0] /=bmag;
  dspec.bhat[1] /=bmag;
  dspec.bhat[2] /=bmag;

  // Create vector and read in DeltaB
  ElsPerProc = dspec.N / size;
  dspec.start = 0;
  for (int p = 0; p < rank; p++) {
    dspec.start += ElsPerProc +
        ((p < dspec.N - ElsPerProc * size) ? 1 : 0);
  }
  dspec.end = dspec.start + ElsPerProc +
      ((rank < dspec.N - ElsPerProc * size) ? 1 : 0);

  //printall("[%d] start = %d end = %d\n", rank, dspec.start, dspec.end);

  deltab = (usedtype*) calloc(dspec.end - dspec.start, sizeof(usedtype));

  MPI_File_seek(fptr, dspec.start*sizeof(double), MPI_SEEK_CUR);

  if (sizeof(usedtype) == sizeof(double)) {
    MPI_File_read(fptr, deltab, dspec.end-dspec.start, MPI_DOUBLE, &status);
    if (flgByteSwap) {
      byteswap((char*)deltab, dspec.end-dspec.start, sizeof(double));
    }
  } else {
    castbuf = (double*) calloc(dspec.end-dspec.start, sizeof(double));
    MPI_File_read(fptr, castbuf, dspec.end-dspec.start, MPI_DOUBLE, &status);
    if (flgByteSwap) {
      byteswap((char*)castbuf, dspec.end-dspec.start, sizeof(double));
    }
    for (int i = 0; i < dspec.end - dspec.start; i++)
      deltab[i] = (float) castbuf[i];
    free(castbuf);
  }

  // Close file pointer
  MPI_File_close(&fptr);

  // set offsets
  dspec.yoffset = dspec.size[0];
  dspec.zoffset = dspec.size[0] * dspec.size[1];

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


  usedtype tau = 0.15, alpha = 0.75, beta = 0.25;

  usedtype *cylColumns;
  int *FGindices; // indices corresponding to foreground elements

  usedtype *new_x;

  int dStart, dEnd;
  int dN;  
  int recvcounts, displs;

  //stopping criteria
  usedtype relative_threshold = 1e-6;
  usedtype absolute_threshold = 1e-16;
  int max_iters = 1000;
  double rms_x, rms_diff_x; //, old_rms_diff_x;
  int iteration = 0;

  //int o, ox, oy, oz;
  //int p, px, py, pz;
  //int rx, ry, rz;

  // int roffset;

  MPI_File fptr;
  char *fname;

  time_t tsecs;
  time_t tmins;
  time_t tStart, tEnd, tsave;
  time_t tIterStart1, tIterStart2, tIterEnd1, tIterEnd2, tRms, tOverhead;

  bool flg;

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
  dStart = dspec.start;
  dEnd = dspec.end;
  dN = dEnd - dStart;

  //==================================================================================================================
  // Create Arrays....
  P = new Problem(kernel, dspec, tau, alpha, beta);

  //==================================================================================================================
  // Create foreground indices array and initialise x
  if (rank==0) printroot("   Creating foreground indices array ...\n");
  FGindices = (int*) calloc(dspec.N, sizeof(int));
  int o = 0;
  for (int p = 0; p < dspec.N; p++) {
    if (mask[p]) {
      FGindices[o] = p;
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
    MPI_File_read(fptr, P->x, dspec.nFG, MPI_USEDTYPE, MPI_STATUS_IGNORE);

    MPI_File_close(&fptr);

    iteration++;
  }
  else {
    int o = 0;
    for (int p = 0; p < dspec.N; p++) {
      if (mask[p]) {
        P->x[o] = chi[p];
        o++;
      }
    }
  }

  // Create new_x array, OK implementation reallocates on every iteration for some reason
  if (rank==0) printroot("   Creating new_x array ...\n");
  new_x = (usedtype*) calloc(dspec.nFG, sizeof(usedtype));
  memset(new_x, 0, dspec.nFG * sizeof(usedtype));

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
    for (int p = dStart; p < dEnd; p++) {
      
      int pz = p / dspec.zoffset;
      int py = (p - pz * dspec.zoffset) / dspec.yoffset;
      int px = p - py * dspec.yoffset - pz * dspec.zoffset;
      
      for (o = 0; o < dspec.nFG; o++) {
        int mix = kernel.modelmap.mask[FGindices[o]];
        if (mix != -1) {
          int oz = FGindices[o] / dspec.zoffset;
          int oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
          int ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;
          
          int rx = px - ox + kernel.halfsize;
          int ry = py - oy + kernel.halfsize;
          int rz = pz - oz + kernel.halfsize;          
          
          cylColumns[mix * dN + p - dStart] = kernel.Get(rx, ry, rz, FGindices[o]);
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
  cl_enqueue_write(P->cl, P->cl_skernel, kernel._nnz * sizeof(usedtype), kernel.skernel);
  cl_enqueue_write(P->cl, P->cl_mask, sizeof(int) * P->dspec.N, kernel.modelmap.mask);
  cl_enqueue_write(P->cl, P->cl_deltab, P->cl_size_n, DeltaB);
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
  bool first = true;

  do {

    tIterStart1 = MPI_Wtime();
    if (rank==0) printroot("      Iteration %d: ", iteration);

    // Ax_b = A * x - b

#ifdef USE_OPENCL
    P->profile1.kern_time = 0;
    cl_size(P->cl, P->dN, 0, P->threads, rank);  //Resize
    if (first)
    {
      // Write initial x buffer
      cl_enqueue_write(P->cl, P->cl_x, P->cl_size_fg, P->x);
      first = false;
    }
    cl_set_arg(cl, P->kernel_iterate1, 14, P->dspec.start);
    cl_set_arg(cl, P->kernel_iterate1, 15, P->dspec.end);
    // Perform the operations
    cl_set_arg(cl, P->kernel_zero, 4, P->dN);
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
#elif defined(USE_OPENMP)
    // if (omp_get_dynamic()) {
    //   printroot ("dynamic threads enabled in OPENMP\n");
    // }

    memset(P->Ax_b, 0, (P->dN) * sizeof(usedtype));
    memset(P->Dx, 0, (P->dN) * sizeof(usedtype));
#pragma omp parallel for schedule(dynamic) num_threads(2)
    {
      omp_set_num_threads(2);
    printroot("num threads = %d\n",omp_get_num_threads());
    for (o = 0; o < dspec.nFG; o++) {
      int ox, oy, oz, px, py, pz, rx, ry, rz;
      
      if (P->x[o]) {
        
        oz = FGindices[o] / dspec.zoffset;
        oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
        ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;

        // TODO(timseries): find out why these were missing from OK's implementation.
        pz = dStart / dspec.zoffset;
        py = (dStart - pz * dspec.zoffset) / dspec.yoffset;
        px = dStart - py * dspec.yoffset - pz * dspec.zoffset - 1;			    

        ry = py - oy;
        rz = pz - oz;
        
        for (int p = dStart; p < dEnd; p++) {
          int size = dspec.size[1] * dspec.size[0];
          pz = p / size;
          py = (p - pz*size) / dspec.size[0];
          px = p - pz*size - py * dspec.size[0];

          rx = px - ox;
          ry = py - oy;
          rz = pz - oz;
          
          if (rx >= -kernel.halfsize && rx <= kernel.halfsize && ry >= -kernel.halfsize && ry <= kernel.halfsize && rz >= -kernel.halfsize && rz <= kernel.halfsize) {
            // Linear system
            int mix = kernel.modelmap.mask[FGindices[o]];
            if (mix == -1) { // spherical kernel
              P->Ax_b[p - dStart] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * P->x[o];
            }
            else if (P->PreCalcCylinders) {
              P->Ax_b[p - dStart] += cylColumns[mix*dN + p-dStart] * P->x[o];
            }
            else if (rx == 0 && ry == 0 && rz == 0) {
              P->Ax_b[p - dStart] += kernel.ctr[mix] * P->x[o];
            }
            else {
              P->Ax_b[p - dStart] += kernel.GetCyl(mix, rx, ry, rz) * P->x[o];
            }
            
            // Laplacian
            Laplacian(rx, ry, rz, P->Dx, P->x, p-dStart, o);
          }
        }
      }
    }
    }
#pragma omp parallel for schedule(dynamic)
    {
    for (int p = 0; p < P->dN; p++) {
      P->Ax_b[p] -= deltab[p];
    }
    }
#else
    memset(P->Ax_b, 0, (P->dN) * sizeof(usedtype));
    memset(P->Dx, 0, (P->dN) * sizeof(usedtype));

    //LWIterate(kernel, P, P->Ax_b, P->x, -1, P->Dx, P->x);
#if 1
    for (o = 0; o < dspec.nFG; o++) {
      int ox, oy, oz, px, py, pz, rx, ry, rz;
      
      if (P->x[o]) {
        
        oz = FGindices[o] / dspec.zoffset;
        oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
        ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;

        // TODO(timseries): find out why these were missing from OK's implementation.
        pz = dStart / dspec.zoffset;
        py = (dStart - pz * dspec.zoffset) / dspec.yoffset;
        px = dStart - py * dspec.yoffset - pz * dspec.zoffset - 1;			    

        ry = py - oy;
        rz = pz - oz;
        
        for (int p = dStart; p < dEnd; p++) {
          int size = dspec.size[1] * dspec.size[0];
          pz = p / size;
          py = (p - pz*size) / dspec.size[0];
          px = p - pz*size - py * dspec.size[0];

          rx = px - ox;
          ry = py - oy;
          rz = pz - oz;
          
          if (rx >= -kernel.halfsize && rx <= kernel.halfsize && ry >= -kernel.halfsize && ry <= kernel.halfsize && rz >= -kernel.halfsize && rz <= kernel.halfsize) {
            // Linear system
            int mix = kernel.modelmap.mask[FGindices[o]];
            if (mix == -1) { // spherical kernel
              P->Ax_b[p - dStart] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * P->x[o];
            }
            else if (P->PreCalcCylinders) {
              P->Ax_b[p - dStart] += cylColumns[mix*dN + p-dStart] * P->x[o];
            }
            else if (rx == 0 && ry == 0 && rz == 0) {
              P->Ax_b[p - dStart] += kernel.ctr[mix] * P->x[o];
            }
            else {
              P->Ax_b[p - dStart] += kernel.GetCyl(mix, rx, ry, rz) * P->x[o];
            }
            
            // Laplacian
            Laplacian(rx, ry, rz, P->Dx, P->x, p-dStart, o);
          }
        }
      }
    }
#endif
 
    for (int p = 0; p < P->dN; p++) {
      P->Ax_b[p] -= deltab[p];
    }
#endif

    tIterEnd1 = MPI_Wtime();
    tsecs = tIterEnd1 - tIterStart1;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    if (rank==0) printroot("(%ldmin %lds)", tmins, tsecs);
    tIterStart2 = MPI_Wtime();

    // AtAx_b = A' * Ax_b
    // DtDx = D' * Dx
#ifdef USE_OPENCL
    P->profile2.kern_time = 0;
    cl_size(P->cl, P->dspec.nFG, 0, P->threads, rank);  //Resize
    //Split the job up if requested
    int BLOCK = (dN / P->divide + 0.5);
    //printf("iterate2 0 ");
    for (int start=dStart; start<P->dEnd; start += BLOCK)
    {
      int end = start + BLOCK;
      if (end > dEnd) end = dEnd;
      //printf("%d ", end);
      cl_set_arg(cl, P->kernel_iterate2, 15, start);
      cl_set_arg(cl, P->kernel_iterate2, 16, end);
      cl_set_arg(cl, P->kernel_iterate2, 17, dStart);
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
#else
    memset(P->AtAx_b, 0, P->dspec.nFG*sizeof(usedtype));
    memset(P->DtDx, 0, P->dspec.nFG*sizeof(usedtype));

    LWIterate(P->AtAx_b, P->Ax_b, 1, P->DtDx, P->Dx);
#if 0
    for (o = 0; o < dspec.nFG; o++) {
      int ox, oy, oz, px, py, pz, rx, ry, rz;
      
      oz = FGindices[o] / dspec.zoffset;
      oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
      ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;
      
      pz = dStart / dspec.zoffset;
      py = (dStart - pz * dspec.zoffset) / dspec.yoffset;
      px = dStart - py * dspec.yoffset - pz * dspec.zoffset - 1;      
      
      ry = py - oy;
      rz = pz - oz;
      
      //roffset = (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset;
      
      for (int p = dStart; p < dEnd; p++) {
        
        px++;
        if (px == dspec.size[0]) {
          px = 0;
          py++;
          if (py == dspec.size[1]) {
            py = 0;
            pz++;
            rz = pz - oz;
          }
          ry = py - oy;
        }        
        rx = px - ox;
     
        if (rx >= -kernel.halfsize && rx <= kernel.halfsize && ry >= -kernel.halfsize && ry <= kernel.halfsize && rz >= -kernel.halfsize && rz <= kernel.halfsize) {
     
          // Linear system
          int mix = kernel.modelmap.mask[FGindices[o]];

          if (mix == -1) { // spherical kernel
            P->AtAx_b[o] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * P->Ax_b[p - dStart];
          }
          else if (P->PreCalcCylinders) {
            P->AtAx_b[o] += cylColumns[mix*dN + p-dStart] * P->Ax_b[p - dStart];
          }
          else if (rx == 0 && ry == 0 && rz == 0) {
            P->AtAx_b[o] += kernel.ctr[mix] * P->Ax_b[p - dStart];
          }
          else {
            P->AtAx_b[o] += kernel.GetCyl(mix, rx, ry, rz) * P->Ax_b[p - dStart];
          }
  
          // Laplacian
          Laplacian(rx, ry, rz, P->DtDx, P->Dx, o, p-dStart);
        }
      }
    }
#endif
#endif

    tIterEnd2 = MPI_Wtime();
    tsecs = tIterEnd2 - tIterStart2;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    if (rank==0) printroot(" (%ldmin %lds)", tmins, tsecs);

    // Reduce AtAx_b
    tOverhead = MPI_Wtime();

    //if (rank==0) printroot("      reducing AtAx_b ...\n");
    MPI_Allreduce(MPI_IN_PLACE, P->AtAx_b, dspec.nFG, MPI_USEDTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, P->DtDx, dspec.nFG, MPI_USEDTYPE, MPI_SUM, MPI_COMM_WORLD);

    tsecs = MPI_Wtime() - tOverhead;
    tmins = floor(tsecs/60);
    tsecs -= tmins*60;
    if (rank==0) printroot(" (%ldmin %lds)", tmins, tsecs);


    //		if (rank == 0) {
    //			sprintf(fname, "%s/DtDx.bin", outdir);
    //			MPI_File_open(PETSC_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
    //			MPI_File_write_ordered(fptr, DtDx, dspec.nFG, MPIU_SCALAR, MPI_STATUS_IGNORE);
    //			MPI_File_close(&fptr);
    //
    //			sprintf(fname, "%s/AtAx_b.bin", outdir);
    //			MPI_File_open(PETSC_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
    //			MPI_File_write_ordered(fptr, AtAx_b, dspec.nFG, MPIU_SCALAR, MPI_STATUS_IGNORE);
    //			MPI_File_close(&fptr);
    //		}

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
      usedtype new_x = P->x[o] - alpha * tau * P->AtAx_b[o] - beta * tau * P->DtDx[o];
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
    //old_rms_diff_x = rms_diff_x;

    // Copy results to x
    //    memcpy(x, new_x, dspec.nFG*sizeof(usedtype));

    // periodic saves just in case the job gets interrupted for any reason...
    tEnd = MPI_Wtime();
    if (rank == 0 && (tEnd > tsave || iteration < 2 || iteration == max_iters - 1)) {
      MPI_File fptr;
      sprintf(fname, "%s/x_iter%06d.bin", myout.outdir, iteration);
      MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
      MPI_File_set_size(fptr, 0);
      MPI_File_write(fptr, &iteration, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fptr, P->x, dspec.nFG, MPI_USEDTYPE, MPI_STATUS_IGNORE);
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
  memset(chi, 0, dspec.N*sizeof(usedtype));

  for (o = 0; o < dspec.nFG; o++) {
    chi[FGindices[o]] = P->x[o];
  }

    wall_time = MPI_Wtime() - wall_time;
#ifdef USE_OPENCL
  if (rank == 0)
  {
    printf("\n#procs, threads, per GPU size, global size, problem size, kernel1 time, kernel2 time, kernel3 time, wall time\n");
    printf("%d, %d, %d, %d, %d, ", size, P->cl->nthreads, P->cl->global_size, 0, P->cl->problem_size);
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
void Process::LWIterate(usedtype* LHS, usedtype* RHS, int dir, usedtype* LapLHS, usedtype* LapRHS) {
//Outer loop
for (int o = 0; o < P->dspec.nFG; o++) {
  if (dir < 0 && !P->x[o]) continue;  //Skip zero

  int oz = P->FGindices[o] / P->dspec.zoffset;
  int oy = (P->FGindices[o] - oz * P->dspec.zoffset) / P->dspec.yoffset;
  int ox = P->FGindices[o] - oy * P->dspec.yoffset - oz * P->dspec.zoffset;

  int pz0 = P->dspec.start / P->dspec.zoffset;
  int py0 = (P->dspec.start - pz0 * P->dspec.zoffset) / P->dspec.yoffset;
  int px0 = P->dspec.start - py0 * P->dspec.yoffset - pz0 * P->dspec.zoffset - 1;      
#if 1
  int px = px0, py = py0, pz = pz0;
  int ry = py0 - oy;
  int rz = pz0 - oz;
  //(Removed(OK)int roffset = (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset;
#endif

  //Inner loop
  for (int p = P->dspec.start; p < P->dspec.end; p++) {
#if 1
    px++;
    if (px == P->dspec.size[0]) {
      px = 0;
      py++;
      if (py == P->dspec.size[1]) {
        py = 0;
        pz++;
        rz = pz - oz;
      }
      ry = py - oy;
    }
    int rx = px - ox;
#else
    //Directly calculating from index is nearly twice as slow (divisions)
    //Implemented for testing as necessary for part of GPU calc
    int R = P->dspec.size[1];
    int C = P->dspec.size[0];
    int RC = R * C;
    int pz = p / RC;
    int pzRC = pz*RC;
    int py = (p - pzRC) / R;
    int px = p - pzRC - py * R;

    int rx = px - ox;
    int ry = py - oy;
    int rz = pz - oz;
#endif

    if (rx >= -kernel.halfsize && rx <= kernel.halfsize && ry >= -kernel.halfsize && ry <= kernel.halfsize && rz >= -kernel.halfsize && rz <= kernel.halfsize) {
      // Linear system
      int mix = kernel.modelmap.mask[P->FGindices[o]];
      int oo = o;
      int pp = p-P->dspec.start;
      if (dir<0) {
        oo = pp;
        pp = o;
      }

      if (mix == -1) { // spherical kernel
        LHS[oo] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * RHS[pp];
      }
      else if (P->PreCalcCylinders) {
        LHS[oo] += P->cylColumns[mix*P->dN + p-P->dspec.start] * RHS[pp];
      }
      else if (rx == 0 && ry == 0 && rz == 0) {
        LHS[oo] += kernel.ctr[mix] * RHS[pp];
      }
      else {
        LHS[oo] += kernel.GetCyl(mix, rx, ry, rz) * RHS[pp];
      }

      // Laplacian
      Laplacian(rx, ry, rz, LapLHS, LapRHS, oo, pp);
    }
  }
}
}
void Process::Laplacian(int rx, int ry, int rz, usedtype* LHS, usedtype* RHS, int o, int p) {
  static int dn = -1; //kernel.halfsize - 1;
  static int dz = 0; //kernel.halfsize;
  static int dp = 1; //kernel.halfsize + 1;
  static usedtype D3_96 = 3.0/96.0;
  static usedtype D10_96 = 10.0/96.0;

  // Laplacian
  if (rx == dn){
    if (ry == dn) {
      if (rz == dz)
        LHS[o] += D3_96 * RHS[p];
    }
    else if (ry == dz) {
      if (rz == dn)
        LHS[o] += D3_96 * RHS[p];
      else if (rz == dz)
        LHS[o] += D10_96 * RHS[p];
      else if (rz == dp)
        LHS[o] += D3_96 * RHS[p];
    }
    else if (ry == dp) {
      if (rz == dz)
        LHS[o] += D3_96 * RHS[p];
    }
  }
  else if (rx == dz) {
    if (ry == dn) {
      if (rz == dn)
        LHS[o] += D3_96 * RHS[p];
      else if (rz == dz)
        LHS[o] += D10_96 * RHS[p];
      else if (rz == dp)
        LHS[o] += D3_96 * RHS[p];
    }
    else if (ry == dz) {
      if (rz == dn)
        LHS[o] += D10_96 * RHS[p];
      else if (rz == dz)
        LHS[o] -= RHS[p];
      else if (rz == dp)
        LHS[o] += D10_96 * RHS[p];
    }
    else if (ry == dp) {
      if (rz == dn)
        LHS[o] += D3_96 * RHS[p];
      else if (rz == dz)
        LHS[o] += D10_96 * RHS[p];
      else if (rz == dp)
        LHS[o] += D3_96 * RHS[p];
    }
    
  }
  else if (rx == dp) {
    if (ry == dn) {
      if (rz == dz)
        LHS[o] += D3_96 * RHS[p];
    }
    else if (ry == dz) {
      if (rz == dn)
        LHS[o] += D3_96 * RHS[p];
      else if (rz == dz)
        LHS[o] += D10_96 * RHS[p];
      else if (rz == dp)
        LHS[o] += D3_96 * RHS[p];
    }
    else if (ry == dp) {
      if (rz == dz)
        LHS[o] += D3_96 * RHS[p];
    }            
  }
}

bool Process::WriteOut(){
  myout.LocalArray(0, chi, 3, dspec.size, "chi");
  if (rank == 0) {
    filepath = (char*) calloc(strlen(myout.outdir) + 9, sizeof(char));
    sprintf(filepath, "%s/chi.bin", myout.outdir);
    MPI_File_open(MPI_COMM_SELF, filepath,
                  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
    MPI_File_write(fptr, chi, dspec.N, MPI_USEDTYPE, MPI_STATUS_IGNORE);
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
