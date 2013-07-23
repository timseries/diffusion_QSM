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

/*! \file modelmap.cc
  \brief ModelMap class file.

  Implementation of the ModelMap class.
*/

#include "hybrid/modelmap.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

#include "hybrid/util.h"

using ::util::checkEndianness;
using ::util::checkVersion;
using ::util::byteswap;

ModelMap::ModelMap() {
  mask = NULL;
  x = NULL;
  y = NULL;
  z = NULL;
  ncyls = 0;
}
ModelMap::~ModelMap() {}
// TODO(timseries): The sizeof(type)'s would be best replaced with a constant...
bool ModelMap::Process(DataSpec &dspec, const ArgHandler &arghandler) {
  int err;
  char *filepath;
  MPI_File fptr;
  double buf[5];
  bool flgByteSwap;
  MPI_Offset offset, disp;
  bool flg;
  Real threshold;
  double *mx, *my, *mz;

  //  if (rank==0) printroot("Loading model map ...\n");
  flg = arghandler.GetArg("-modelmap", filepath);
  if (!flg) {
    //    if (rank==0) printroot("Model map file not specified\n");
    return false;
  }
  //  if (rank==0) printroot("   file: %s\n", filepath);

  // open file for reading
  err = MPI_File_open(MPI_COMM_SELF,
                      filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
  //if (err) {
  //  if (rank==0) printroot("Could not open model mask file"); return false;}

  // check endianness
  checkEndianness(fptr, flgByteSwap);

  // check version
  checkVersion(fptr, flgByteSwap);

  // read in header
  MPI_File_read(fptr, buf, 5, MPI_DOUBLE, MPI_STATUS_IGNORE);

  // Check parameters match
  if (flgByteSwap) byteswap(reinterpret_cast<char*>(buf), 5, sizeof(double));

  if (dspec.size[0] != buf[1]) {
    //    if (rank==0) printroot("Size of first dimension of ModelMask does not match DeltaB");
    return false; }
  if (dspec.size[1] != buf[2]) {
    //    if (rank==0) printroot("Size of second dimension of ModelMask does not match DeltaB");
    return false; }
  if (dspec.size[2] != buf[3]) {
    //    if (rank==0) printroot("Size of third dimension of ModelMask does not match DeltaB");
    return false; }
  if (buf[4] != 3) {
    //    if (rank==0) printroot("Size of fourth dimension of ModelMask is not 3");
    return false; }
  if (buf[0] != buf[1] * buf[2] * buf[3] * buf[4]) {
    //    if (rank==0) printroot("Number of elements ModelMask does not match array sizes");
    return false; }

  // Load model data
  mx = reinterpret_cast<double*>(calloc(dspec.N, sizeof(double)));
  my = reinterpret_cast<double*>(calloc(dspec.N, sizeof(double)));
  mz = reinterpret_cast<double*>(calloc(dspec.N, sizeof(double)));

  if (mx == NULL || my == NULL || mz == NULL) {
    //    if (rank==0) printroot("Not enough memory available.");
    return false;
  }

  MPI_File_read(fptr, mx, dspec.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_read(fptr, my, dspec.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_read(fptr, mz, dspec.N, MPI_DOUBLE, MPI_STATUS_IGNORE);

  if (flgByteSwap) byteswap(reinterpret_cast<char*>(mx),
                            dspec.N, sizeof(double));
  if (flgByteSwap) byteswap(reinterpret_cast<char*>(my),
                            dspec.N, sizeof(double));
  if (flgByteSwap) byteswap(reinterpret_cast<char*>(mz),
                            dspec.N, sizeof(double));
  MPI_File_close(&fptr);

  // GET THRESHOLD
  flg = arghandler.GetArg("-mt", threshold);
  if (!flg) threshold = 0.2;

  //  if (rank==0) printroot("   threshold = %0.3f\n", threshold);

  // CREATE MASK
  mask = reinterpret_cast<int*>(calloc(dspec.N, sizeof(int)));

  int ix = 0;
  for (int i = 0; i < dspec.N; i++) {
    mask[i] = (sqrt(mx[i]*mx[i] + my[i]*my[i] + mz[i]*mz[i]) >
               threshold) ? ix++ : -1;
  }

  // GATHER X Y AND Z ARRAYS
  ncyls = ix;
  //  if (rank==0) printroot("   number cylinders = %d\n", ncyls);
  x = reinterpret_cast<Real*>(calloc(ncyls, sizeof(Real)));
  y = reinterpret_cast<Real*>(calloc(ncyls, sizeof(Real)));
  z = reinterpret_cast<Real*>(calloc(ncyls, sizeof(Real)));
  for (int i = 0; i < dspec.N; i++) {
    if (mask[i] >= 0) {
      x[mask[i]] = mx[i];
      y[mask[i]] = my[i];
      z[mask[i]] = mz[i];
    }
  }
  //  if (rank==0) printroot("   ix = %d\n", ix);
  free(mx);
  free(my);
  free(mz);

  return true;
}

void ModelMap::close() {
  if (mask != NULL) free(mask);
  if (x != NULL) free(x);
  if (y != NULL) free(y);
  if (z != NULL) free(z);
}
