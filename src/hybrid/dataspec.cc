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
  orb_divisions=NULL;
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
void DataSpec::Create(double* buf,int rank, int mpi_world_size) {
  // process header into relevant variables
  //buf is a file buffer created by loadDeltaB
  int ElsPerProc=0;
  double bmag=0;
  int *workmatrix = NULL;
  int px,pypz;
  this->mpi_world_size=mpi_world_size;
  this->rank=rank;
  N = (int) buf[0];
  size[0] = (int) buf[1];
  size[1] = (int) buf[2];
  size[2] = (int) buf[3];
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
  //orb arrray divisions based on the intialized world size
  orb_divisions=(int*) calloc(mpi_world_size-2, sizeof(int));

  if (0) {
    
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
    workmatrix=(int*) calloc(N, sizeof(int));
    int rx=0;
    int ry=0;
    int rz=0;
    int px=0;
    int py=0;
    int pz=0;

    int distance_to_center=0;
    int reference=0;
    int prevworkmatrix=0;
    int center = N/2;
    int centerz = center / zoffset;
    int centery = (center - centerz * zoffset) / yoffset;
    int centerx = center - centery * yoffset - centerz * zoffset - 1;			    
    //set the start and finish using blocks proportional to the amount of work...
    //calculate the total work required, which is invsersely related to the distance to the center
    for (int p=0; p < N; p++){
      pz = p / zoffset;
      py = (p - pz * zoffset) / yoffset;
      px = p - py * yoffset - pz * zoffset - 1;			    
      rx=px-centerx;
      ry=py-centery;
      rz=pz-centerz;      
      distance_to_center=ceil(pow((pow(rx,2)+pow(ry,2)+pow(rz,2)),1.0/2.0));
      // printroot("distance_to_center: %d\n",distance_to_center);
      if (p==0) {
        reference=distance_to_center;
      }
      workmatrix[p]=prevworkmatrix+abs(reference-distance_to_center);
      prevworkmatrix=workmatrix[p];
      // printroot("workmatrix[]]: %d\n",workmatrix[p]);
    }

//do orthognonal recursive bisection
  PartitionByORB(workmatrix);
  }
  range=end-start;
}
void DataSpec::PartitionByORB(int* workmatrix) {
  // printroot("workend: %d, workstart: %d\n", workmatrix[N],workmatrix[0]);
  ORB(0,N-1,0,mpi_world_size-2,workmatrix);
  if (rank==0) {
    start=0;
  } else {
    start=orb_divisions[rank-1]+1;
  }
  if (rank==(mpi_world_size-1)) {
    end=N;
  } else {
    end=orb_divisions[rank];
  }
}
void DataSpec::ORB(int data_start, int data_end, int orb_start,int orb_end,int* workmatrix) {

  //recursively partition the data
  int workdiv=(workmatrix[data_end]-workmatrix[data_start])/2;
  //find the index where the work is divided in half
  int data_div=data_start;
  // printrootl("dataend: %d, datastart: %d, workend: %d, workstart: %d\n", data_end, data_start, workmatrix[data_end],workmatrix[data_start]);
  printroot("data_start: %d, data_end: %d, orb_start: %d, orb_end: %d, workdiv: %d\n", data_start, data_end, orb_start, orb_end, workdiv);
  while (workmatrix[data_div]-workmatrix[data_start] < workdiv) {
    data_div++;
  }
  if (orb_start==orb_end) {
    orb_divisions[orb_start]=data_div;
    printroot("orb_divisions[orb_start]: %d, data_div: %d\n", orb_divisions[orb_start], data_div);
  } else {
    int orb_div=(orb_start+orb_end)/2;
    orb_divisions[orb_div]=data_div;
    ORB(data_start,data_div,orb_start,orb_div-1,workmatrix);
    ORB(data_div,data_end,orb_div+1,orb_end,workmatrix);
  }
}
