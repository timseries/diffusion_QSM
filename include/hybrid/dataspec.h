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

#ifndef INCLUDE_DATASPEC_H_
#define INCLUDE_DATASPEC_H_

#include "hybrid/common.h"

class DataSpec {
 public:
  DataSpec();
  virtual ~DataSpec();

  void Create(double* buf,int rank, int size);
  void PartitionByORB(int* workmatrix);
  void ORB(int data_start, int data_end, int orb_start,int orb_end,int* workmatrix);
  int size[3];  // data size
  int N; // number of foxelsdelmap
  int yoffset;
  int zoffset;
  int mpi_world_size; //!< Brief number of processes, used for splitting up data
  int rank; //!< Brief the rank of the process which has this dataspec
  int* orb_divisions; //!< Brief array of divisions from orb algorithm
  double B0;
  double bhat[3];
  double caxis[3];
  int nBG;
  int nFG;
  unsigned start;
  unsigned end;
  unsigned range;
};
#endif  // INCLUDE_DATASPEC_H_
