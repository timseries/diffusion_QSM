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

/*! \file process.h
    \brief Process class definitions file.

*/

#ifndef INCLUDE_PROCESS_H_
#define INCLUDE_PROCESS_H_

#include "hybrid/basictypes.h"
#include "hybrid/dataspec.h"
#include "hybrid/arghandler.h"
#include "hybrid/kernel.h"
#include "hybrid/output.h"
#include "hybrid/problem.h"
#include "hybrid/util.h"


class Process {
 public:
  Process();
  virtual ~Process();
  bool Init(int argc, char** args);
  void HandleArgs(int argc, char** args);
  void StartMPI(int argc, char** args);
  bool loadDeltaB();
  bool loadMask();
  bool FullPass();
  inline void MultAdd(Real* result_fidelity, Real* result_reguliarizer,
                      Real* multiplicand_fidelity,Real* multiplicand_regularizer, 
                      Real* addend, bool dir);
  Real Mult();
  bool WriteOut();
  bool CleanUp();
  models model;
  bool *mask;
  Kernel kernel;
  DataSpec dspec;
  Problem *P;
  
  Real threshold;
  Real *deltab;
  Real *chi;
  
  int rank;
  int size;

  Output myout;
  ArgHandler arghandler;
  char *filepath;
  MPI_File fptr;
};
#endif  // INCLUDE_PROCESS_H_
