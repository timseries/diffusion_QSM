// Copyright (c) 2013, Amanda Ng, Timothy Roberts 
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
// Author: amanda.ng@gmail.com, timothy.daniel.roberts@gmail.com

/*! \file output.h
    \brief Output class definitions file.

*/

#ifndef INCLUDE_OUTPUT_H_
  #define INCLUDE_OUTPUT_H_

#include <mpi.h>

#include "hybrid/common.h"
#include "hybrid/arghandler.h"

/**
* Output class. Used to store store the binary output data \f$\Chi\f$ as well as the binary mask showing locations of foreground voxels (1) or background voxels (0).
*/
class Output {
 public:
/**
* Output constructor.
*/
  Output();
/**
* Output destructor.
*/
  virtual ~Output();
/**
* Method for intiializing the class. Opens matlab and binary files for writing, and redirects standard error (stderr) to an output error file errfile.
* @param &arghandler Object pointer of class Arghandler.
* @param rank Rank of this process.
* @param size Number of MPI processes.
* @see Arghandler
*/
  void Init(const ArgHandler &arghandler, int rank, int size);
/**
* Method for writing a local Real array to a binary file.
* @param onproc The process which performs the write.
* @param array Real array to write to binary file. 
* @param ndims Number of dimensions of array.
* @param dims array of length ndims, each element being the size of array along that dim.
* @param arrayname Name of matlab variable that matlab code will read binary array into.
*/
  void LocalArray(int onproc, Real* array,
                  int ndims, int* dims, const char* arrayname);
/**
* Method for writing a local bool array to a binary file.
* @param onproc The process which performs the write.
* @param array bool array to write to binary file. 
* @param ndims Number of dimensions of array.
* @param dims array of length ndims, each element being the size of array along that dim.
* @param arrayname Name of matlab variable that matlab code will read binary array into.
*/
  void LocalArray(int onproc, bool* array,
                  int ndims, int* dims, const char* arrayname);
/**
* Method for writing a distributed (MPI)Real array to a binary file and matlab code to read this binary file (into matlab).
* @param array Real array to write to binary file. 
* @param ndims Number of dimensions of array.
* @param dims array of length ndims, each element being the size of array along that dim.
* @param arrayname Name of matlab variable that matlab code will read binary array into.
*/
  void DistrArray(Real* array, int localsize,
                  int ndims, int* dims, const char* arrayname);
/**
* Method for closing the binary and matlab files.
*/
  void Close();
  char *outdir; ///< Output directory where binary output files are written.
  int rank; ///< Array size N, 1 where there is a foreground voxel, 0 background
  int size; ///< Number of processes (in the MPI communication world).
 private:
        MPI_File binfile; ///< Temporary binary file handle.
        MPI_File matfile; ///< Temporary matlab file handle.
        FILE *errfile; ///< File which sterr are redirected to.
        bool initialised; ///< Set to true when Init is run successfully (files in open state). Set to false after Close finishes.
        char *tmpstr; ///< Character array used to buffer writes to local/distributed arrays.
};
#endif  // INCLUDE_OUTPUT_H_
