// Copyright (c) 2013, Timothy Roberts, Amanda Ng
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

/*! \file output.h
    \brief Output class definitions file.

*/

#ifndef INCLUDE_OUTPUT_H_
  #define INCLUDE_OUTPUT_H_

#include <mpi.h>

#include "hybrid/basictypes.h"
#include "hybrid/arghandler.h"

class Output {
 public:
  Output();
  virtual ~Output();
  void Init(const ArgHandler &arghandler, int rank, int size);
  void LocalArray(int onproc, usedtype* array,
                  int ndims, int* dims, const char* arrayname);
  void LocalArray(int onproc, bool* array,
                  int ndims, int* dims, const char* arrayname);
  void DistrArray(usedtype* array, int localsize,
                  int ndims, int* dims, const char* arrayname);
  void Close();
  char *outdir;
  int rank;
  int size;
 private:
        MPI_File binfile;
        MPI_File matfile;
        FILE *errfile;
        bool initialised;
        char *tmpstr;
};
#endif  // INCLUDE_OUTPUT_H_