// Copyright (c) 2013, Amanda Ng and Timothy Roberts
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

/*! \file modelmap.h
  \brief ModelMap class definitions file.

*/

#ifndef INCLUDE_MODELMAP_H_
#define INCLUDE_MODELMAP_H_

#include "hybrid/common.h"
#include "hybrid/dataspec.h"
#include "hybrid/arghandler.h"

/**
* ModelMap class. Used to compute and store the binary mask of foreground voxels (mask) and determine the cylinder orientations.
*/
class ModelMap {
 public:
/**
* DataSpec constructor.
*/
  ModelMap();
/**
* DataSpec destructor.
*/
  virtual ~ModelMap();
/**
* Class initialization method. 
* @param &dspec pointer to an object of DataSpec.
* @param &arghandler pointer to an object of Arghandler.
* @return True if successful, false otherwise.
*/
  bool Create(DataSpec &dspec, const ArgHandler &arghandler);
  void close(void);
  int *mask; ///< Array of size N, entries corresponding to spherical voxel (-1) or cylyindrical voxel (n>0, n is the index of the cylinder)
  int ncyls; ///< Number of cylinders in the dataset
  Real *x; ///< Array of size ncyls, the vector of x components of the cylindrical voxel direction
  Real *y; ///< Array of size ncyls, the y component of the cylindrical voxel direction
  Real *z; ///< Array of size ncyls, the z component of the cylindrical voxel direction
};
#endif  // INCLUDE_MODELMAP_H_
