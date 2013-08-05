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


DataSpec::DataSpec(){
  size[0]=0;
  size[1]=0;
  size[2]=0;
  N=0;
  yoffset=0;
  zoffset=0;
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
}
DataSpec::~DataSpec() {}
void DataSpec::Create(double* buf,int rank, int mpi_world_size) {
  // process header into relevant variables
  int ElsPerProc=0;
  double bmag=0;
  N = (int) buf[0];
  size[0] = (int) buf[1];
  size[1] = (int) buf[2];
  size[2] = (int) buf[3];
  B0 = buf[4];
  memcpy(bhat, &buf[5], 3*sizeof(double));

  // normalise bhat
  bmag = sqrt(bhat[0]*bhat[0] +
	      bhat[1]*bhat[1] + bhat[2]*bhat[2]);
  bhat[0] /=bmag;
  bhat[1] /=bmag;
  bhat[2] /=bmag;

  // Create vector and read in DeltaB
  ElsPerProc = N / mpi_world_size;
  start = 0;
  for (int p = 0; p < rank; p++) {
    start += ElsPerProc +
      ((p < N - ElsPerProc * mpi_world_size) ? 1 : 0);
  }
  end = start + ElsPerProc +
    ((rank < N - ElsPerProc * mpi_world_size) ? 1 : 0);
  range=end-start;
}
