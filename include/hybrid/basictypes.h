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

/*! \file basictypes.h
    \brief Define types used in this project. Some are set durin based on command line arguments.

    Define types and constants and macros used in this project. Some are set during compile based on command line arguments.
    Real: the datatype used for numerical computations throughout (MPI, OpenMP). This is specified using a -DUSEDTYPE command line argument.
    models: this typedef specifies different model assumptions.

*/

#ifndef INCLUDE_BASICTYPES_H_
  #define INCLUDE_BASICTYPES_H_

//Types
#ifndef USEDTYPE
  #define Real float
  #define MPI_Real MPI_FLOAT
#else
  #if USEDTYPE == DOUBLE
    #define Real double
    #define MPI_Real MPI_DOUBLE
  #else
    #define Real float
    #define MPI_Real MPI_FLOAT
  #endif  // USEDTYPE == DOUBLE
#endif  // INCLUDE_BASICTYPES_H_

//Constants
#define SIZEOF_CHAR 1
#define LOADDATAVERSION 1

#ifdef USE_OPENMP
  #define OPENMP 1
#else
  #define OPENMP 0
#endif

//Macros
//TODO(timseries): consider inlining these or including these or something else, as per google coding standards on preprocessor macros.

#define printall printf
#define printroot printf

//
// The models typedef used for computing the kernel basis functions.
// MODEL_SPHERICAL: use the spherical model only.
// MODLE_MIXED: use the spherical+cylindrical model.
typedef enum {
  MODEL_SPHERICAL, MODEL_MIXED
} models;
#endif  // INCLUDE_BASICTYPES_H_
