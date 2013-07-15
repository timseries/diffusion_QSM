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

/*! \file kernel.cc
  \brief Kernel class file.

  Implementation of the Kernel class.
*/

#include "hybrid/util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

namespace util {
  // This function checks the first flag in the data file for whether an endian change has occurred.
  // It compares the endianness against the local process, not the server process.
  // Therefore subsequent reads of data from the file should either checked and converted before sending
  // to other processes, or sent and then checked and converted.
  // INPUTS:	MPI_File	&fptr
  // OUTPUTS:	bool		&flgByteSwap
void checkEndianness(MPI_File &fptr, bool &flgByteSwap) {

  MPI_Status	status;
  unsigned	endianCheckValue;

  MPI_File_read(fptr, &endianCheckValue, 1, MPI_UNSIGNED, &status);
  if (endianCheckValue == 0x12345678) //Hex value 0x12345678
    flgByteSwap = false;
  else if (endianCheckValue == 0x78563412) //Hex value 0x78563412
    flgByteSwap = true;
  else {
    printroot("Could not recognise endian check flag\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void checkVersion(MPI_File &fptr, bool flgByteSwap) {
  MPI_Status	status;
  unsigned	versionCheckValue;
  bool		check;
  MPI_File_read(fptr, &versionCheckValue, 1, MPI_UNSIGNED, &status);
  if (flgByteSwap) byteswap((char*)&versionCheckValue, 1, sizeof(unsigned));

  if (versionCheckValue != LOADDATAVERSION) {
    printroot("Data file version is not compatible\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void byteswap(char* buf, int buflength, int dataTypeBytes) {
  // buf			 : array
  // buflength	 : number of elements in array
  // dataTypeBytes : number of bytes in datatype
  char swappee[dataTypeBytes];
  char swapped[dataTypeBytes];

  for (int i=0; i<buflength; i++) {
    memcpy(swappee, &buf[i*dataTypeBytes] , dataTypeBytes);
    for (int j=0; j<dataTypeBytes/2; j++) {
      swapped[j] = swappee[dataTypeBytes-1-j];
      swapped[dataTypeBytes-1-j] = swappee[j];
    }
    memcpy(&buf[i*dataTypeBytes], swapped, dataTypeBytes);
  }
}
}
