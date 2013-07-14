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

#include "include/process.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

Process::Process() {
  rank = 0;
  size = 0;
  *mask = NULL;
  *deltab = NULL;
  *chi = NULL;
  *filepath = NULL
      }

Process::Init(int argc, char** args) {
  HandleArgs(int argc, char** args);
  StartMPI(int argc, char** args);
}
Process::HandleArgs(int argc, char** args) {
  char *str = NULL;
  arghandler.Init(argc, args);
#ifdef DEBUG
  printroot("   Argument handler initialised\n");
#endif
  // GET THRESHOLD
  if (!arghandler.GetArg("-threshold", threshold)) {
    threshold = 1e-6;
  }
  printroot("   Threshold = %0.3e\n", threshold);
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
        printroot("Unrecognised model");
        return 1;
        break;
    }
  }
  printroot("   Model: %s\n", model == MODEL_SPHERICAL ?
            "spherical" : "mixed");
      
}
Process::StartMPI(int argc, char** args) {
  MPI_Init(&argc, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printroot("\n------------------------------------------\n");
  printroot("MPI Environment\n");
  printroot("Number of processes: %d\n", size);
  printroot("This rank: %d\n", rank);
}        
