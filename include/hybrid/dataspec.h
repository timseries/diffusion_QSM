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
#include "hybrid/arghandler.h"
#include "hybrid/util.h"

/**
* DataSpec class. Used to store the members and methods used in keeping track of local (and global) allocations of  \f$\mathbf{\Delta B}\f$ voxels.
*/
class DataSpec {
 public:
/**
* DataSpec constructor.
*/
  DataSpec();
/**
* DataSpec destructor.
*/
  virtual ~DataSpec();
/**
* Class initialization method. 
* @param &arghandler pointer to an object of Arghandler.
* @param rank the rank of this process.
* @param mpi_world_size the number of MPI processes.
* @return True if successful, false otherwise.
*/
  bool Create(const ArgHandler &arghandler, int rank, int mpi_world_size);
/**
* Method to allocate partitions to processes in one of two modes. Should only be called once, and calculates the cumulative work of the workmatrix (initialized by the first MultAdd in process).
* @param orb_flag True if using ORB to allocate portions of \f$\mathbf{\Delta B}\f$ to processes, false otherwise.
* @see Process::MultAdd()
*/
  void AllocatePartitions(bool orb_flag);
/**
* Method to calculate partitions based on the size of an N-sized workmatrix using orthogonal recursive bisection concept.
*/
  void PartitionByORB(); 
  int size[3];  ///< Dimensions of the signal (x,y,z).
  int N; ///< Total number of voxels in the \f$\mathbf{\Delta B}\f$ dataset.
#ifdef USE_FOURIER_SPHERES
  int N_fft; ///> Number of elements in a Fourier transformed signal of length \f$ N \f$.
#endif
  int yoffset; ///< Distance between succesive y indices in a column-major 3d array.
  int zoffset; ///< Distance between succesive z indices in a column-major 3d array.
  int mpi_world_size; ///< Number of processes, used for splitting up data
  int rank; ///< The rank of the process which has this dataspec.
  int *orb_divisions; ///< Array of divisions from orb algorithm.
  int orb_divisions_size; ///< Size of the array orb_dvisions
  int *workmatrix; ///<  Array (size N) with comulative work required to process each voxel in \f$\mathbf{\Delta B}\f$.
  double B0; ///< External static magnetic field \f$B_0\f$. Scalar.
  double bhat[3]; ///< Magnetic field flux density \f$\mathbf{B}\f$. Vector.
  double caxis[3]; ///< Reference cylinder orientation \f$\mathbf{\hat c}\f$. Vector.
  int nBG; ///< Number of background voxels. Scalar.
  int nFG; ///< Number of foreground voxels. Scalar.
  unsigned start; ///< Starting index of a local process' portion of \f$\mathbf{\Delta B}\f$.
  unsigned end; ///< Ending index of a local process' portion of \f$\mathbf{\Delta B}\f$.
  unsigned range; ///< Size of a local process' portion of \f$\mathbf{\Delta B}\f$. End-start.
};
#endif  // INCLUDE_DATASPEC_H_
