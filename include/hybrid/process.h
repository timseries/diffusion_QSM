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

#include "hybrid/common.h"
#include "hybrid/dataspec.h"
#include "hybrid/arghandler.h"
#include "hybrid/kernel.h"
#include "hybrid/output.h"
#include "hybrid/problem.h"
#include "hybrid/util.h"

/**
* Process class. Contains the data members and methods for solving the linear system \f$ \mathbf{A \Delta \chi \!= \! \Delta B} \f$ where \f$ \mathbf{A} \f$ is the convolution kernel of spherical and cylinderical voxels.
*/
class Process {
 public:
/**
* Process constructor.
*/
  Process();
/**
* Process destructor.
*/
  virtual ~Process();
/**
* Method for intiializing the class. Initializes the member arrays and MPI functions and processes the command line arguments.
* @param argc Number of standard command line arguments.
* @param args Command line arguments.
* @see Arghandler
* @see Process::HandleArgs
* @see Process::StartMPI
*/
  bool Init(int argc, char** args);
/**
* Method for passing the command line arguments to the Arghandler arghanlder and setting some of the solver parameters.
* @param argc Number of standard command line arguments.
* @param args Command line arguments.
* @see Arghandler
*/
  void HandleArgs(int argc, char** args);
/**
* Method for initializing MPI, using special command line arguments if provided.
* @param argc Number of standard command line arguments.
* @param args Command line arguments.
*/
  void StartMPI(int argc, char** args);
/**
* Loads \f$\mathbf{\Delta B}\f$ from file into class member deltab.
* @return False if this method failed.
*/
  bool loadDeltaB();
/**
* Loads foreground/background binary mask from file into class member mask.
* @return False if this method failed.
*/
  bool loadMask();
/**
* The main solve method. Calls the load methods and runs the Landweber iterations via Multadd
* Landweber iteration: 
* \f$ \mathbf{x}_{n+1} \!=\! \mathbf{x}_{n} - \alpha\mathbf{A'(A*x}_{n} - \mathbf{b)} - \beta * \mathbf{D'(Dx}_{n}\mathbf{)} \f$
* where \f$ \mathbf{b \!=\! \Delta B} \f$ (deltab),
* \f$ \mathbf{x\!=\! \Delta \chi} \f$ (x), 
* \f$ \alpha \!=\! 0.5\tau \f$,
* \f$ \beta \!=\! 0.5\tau \f$,
* \f$ \tau \!=\! \frac{2}{\|\mathbf{A}\|} \f$, and
* \f$ \mathbf{D} \f$ represents a 3D laplacian filter.
*
* @return False if this method failed.
* @see Process::Multadd
*/
  bool FullPass();
/**
* Performs a mult add operation. This is where most of the computatation takes place.
* Returns an array result  (result_fidelity) which is: multiplier*multipicand [+addend].
* Dir determines the direction (true->forward,false->backward), where
* Forward means compute \f$ \mathbf{Ax-b} \f$ (Ax_b) and 
* Backward means compute \f$ \mathbf{A'(Ax-b)} \f$ (AtAx_b).
* This method supports three mutually exclusive modes:
* 0. No optimization.
* 1. Threaded (USE_OPENMP) for Bluegene hardware
* 2. Vectorized (USE_OPENCL) for GPU hardware (e.g. Kepler). Threaded (US
* Additionally, there is a Fourier optimization (FOURIER_SPHERES) which can be used with the mode #1.
*/
  void MultAdd(Real* result_fidelity, Real* result_reguliarizer,
                      Real* multiplicand_fidelity, Real* multiplicand_regularizer, 
                      Real* addend, int* workmatrix, bool dir, int iteration);
/**
* Write the data \f$ \mathbf{\Delta \chi} \f$ to output data file chi.bin.
* @return False if this method failed.
*/
  bool WriteOut();
/**
* Free the memory allocated for the data members of this class.
* @return False if this method failed.
*/
  bool CleanUp();
  models model; ///< Enum of type models for the mode (spherical or spherical+cylindrical).
  bool *mask; ///< Array size N, 1 where there is a foreground voxel, 0 background
  Kernel kernel; ///< Object of class kernel.
  DataSpec dspec; ///< Object of class dspec.
  
  Real threshold; ///< Threshold used to determine the kernel extent.
  Real *deltab; ///< The N-sized array \f$\mathbf{\Delta B}\f$.
  Real *chi; ///< The nFG-sized array \f$ \mathbf{\Delta \chi} \f$.
  
  int rank; ///< The rank of this process.
  int size; ///< The number of processes in the MPI communication world.

  Output myout; ///< Object of class Output for wrtiing \f$ \mathbf{\Delta \chi} \f$ (chi) and mask to binary files.
  ArgHandler arghandler; ///< Object of class Arghandler.
  Problem *P; ///< Object of class Problem.
};
#endif  // INCLUDE_PROCESS_H_
