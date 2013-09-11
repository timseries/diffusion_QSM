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

/*! \file dataspec.h
  \brief Dataspec class declarations file.

*/
#include "hybrid/dataspec.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

using ::util::checkVersion;
using ::util::checkEndianness;
using ::util::byteswap;

DataSpec::DataSpec(){
  size[0]=0;
  size[1]=0;
  size[2]=0;
  N=0;
  yoffset=0;
  zoffset=0;
  orb_divisions=NULL;
  workmatrix=NULL;
  orb_divisions_size=0;
  B0=0;
  bhat[0]=0;
  bhat[1]=0;
  bhat[2]=0;
  caxis[0]=0;
  caxis[1]=0;
  caxis[2]=0;
  nBG=0;
  nFG=0;
  start=0;
  end=0;
  range=0;
  mpi_world_size=0;
  rank=0;
}
DataSpec::~DataSpec() {}
bool DataSpec::Create(const ArgHandler &arghandler, int rank, int mpi_world_size) {
  // process header into relevant variables
  //info about datas (size and number of voxels) determined from dltab
  double bmag=0;
  char *filepath;
  int err=0;
  MPI_File fptr;
  MPI_Status status;
  MPI_Offset offset, disp;
  bool flgByteSwap;
  double buf[11];
  double *castbuf;

  if (rank==0) printroot("Loading DeltaB ...\n");
  if (!arghandler.GetArg("-DeltaB", filepath)) {
    if (rank==0) printroot("Data file not specified\n");
    return false;
  }
  if (rank==0) printroot("   file: %s\n", filepath);
  err = MPI_File_open(MPI_COMM_SELF, filepath,
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
  if (err) {
    printroot("rank:%d Could not open file\n", rank);
    return false;
  }
  if (rank==0) printroot("error status from file open on proc 0: %d\n",err);
  checkEndianness(fptr, flgByteSwap);
  checkVersion(fptr, flgByteSwap);
  if (rank==0) printroot("Reading-in MPI header...\n");
  MPI_File_read(fptr, buf, 11, MPI_DOUBLE, &status);
  if (rank==0) printroot("Byte-swapping header...\n");
  if (flgByteSwap) byteswap((char*)buf,11,sizeof(double));
  // create the dataspec here
  if (rank==0) printroot("Creating dataspec...\n");

  this->mpi_world_size=mpi_world_size;
  this->rank=rank;
  N = (int) buf[0];

  size[0] = (int) buf[1];
  size[1] = (int) buf[2];
  size[2] = (int) buf[3];

  caxis[0] = 1;
  caxis[1] = 0;
  caxis[2] = 0;

  yoffset = size[0]; 
  zoffset = size[0] * size[1]; 

  B0 = buf[4];
  memcpy(bhat, &buf[5], 3*sizeof(double));

  // normalise bhat
  bmag = sqrt(bhat[0]*bhat[0] +
	      bhat[1]*bhat[1] + bhat[2]*bhat[2]);
  bhat[0] /=bmag;
  bhat[1] /=bmag;
  bhat[2] /=bmag;
  err=MPI_File_close(&fptr);
  free(filepath);

  // Create workmatrix
  if (rank==0) printroot("   Creating workmatrix array ...\n");
  workmatrix = (int*) calloc(N, sizeof(int));
  return true;
}

void DataSpec::AllocatePartitions(bool orb_flag) {
  int ElsPerProc=0;
  workmatrix = (int*) calloc(N, sizeof(int));
  orb_divisions_size=mpi_world_size-2;
  orb_divisions=(int*) calloc(orb_divisions_size, sizeof(int));

  if (~orb_flag) {
    
    //set the start and finish using continguous, even-sized blocks
  ElsPerProc = N / mpi_world_size;
  start = 0;
  for (int p = 0; p < rank; p++) {
    start += ElsPerProc +
      ((p < N - ElsPerProc * mpi_world_size) ? 1 : 0);
  }
  end = start + ElsPerProc +
    ((rank < N - ElsPerProc * mpi_world_size) ? 1 : 0);
 } else {
    int prevworkmatrix=0;
    //calculate the cumulative work required in place
    for (int p=0; p < N; p++){
      workmatrix[p]=workmatrix[p]+prevworkmatrix;
      prevworkmatrix=workmatrix[p];
    }
  PartitionByORB();
  }
  range=end-start;
}

void DataSpec::PartitionByORB() {
  //compute average amount of work
  int average_process_work=(int) (workmatrix[N-1]-workmatrix[0])/(orb_divisions_size+2);
  //iterate through the workmatrix and get appropriate bounds for each of the processes
  int orb_start=0;
  int orb_divisions_index=0;
  for (int orb_end=0; orb_end < N; orb_end++) {
    if ((workmatrix[orb_end]-workmatrix[orb_start])>=average_process_work) {
      if (orb_divisions_index>=orb_divisions_size) {
        // finished, last process will get what's leftover...
        break;
      }  
      orb_divisions[orb_divisions_index]=orb_end;
      orb_divisions_index++;
      orb_start=orb_end;
    }
  }
  if (rank==0) {
    start=0;
  } else {
    start=orb_divisions[rank-1];
  }
  if (rank==(mpi_world_size-1)) {
    end=N;
  } else {
    end=orb_divisions[rank];
  }

}

