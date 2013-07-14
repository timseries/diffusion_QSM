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

/*! \file output.cc
  \brief Output class file.

  Implementation of the Output class.
*/

#include "include/output.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <mpi.h>

Output::Output() {
  initialised = false;
}

Output::Init() {
  bool flg;
  // TODO(timseries): replace sprintf with snprintf which is more secure.
  flg = arghandler.GetArg("-out", outdir);
  if (!flg) {
    outdir = reinterpret_cast<char*>(calloc(7, SIZEOF_CHAR));
    sprintf(outdir, "output");
  }
  tmpstr = reinterpret_cast<char*>(calloc(1024, SIZEOF_CHAR));
  // redirect stderr to error files
  sprintf(tmpstr, "%s/err%03d.txt", outdir, rank);
  if ((errfile = freopen(tmpstr, "w", stderr)) == NULL) {
    printall("[%d] Could not redirect stderr\n", rank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#ifdef DEBUG
  printroot("   Output::Init : redirected stderr\n");
#endif
  // Open binary file
  if (rank == 0) {
    sprintf(tmpstr, "%s/out.bin", outdir);
    MPI_File_open(MPI_COMM_SELF, tmpstr, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &binfile);
#ifdef DEBUG
    printroot("   Output::Init : opened binary file\n");
#endif
    MPI_File_set_size(binfile, 0);
    // Open matlab file on process 0
    sprintf(tmpstr, "%s/out.m", outdir);
    MPI_File_open(MPI_COMM_SELF, tmpstr, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &matfile);
    MPI_File_set_size(matfile, 0);
#ifdef DEBUG
    printroot("   Output::Init : opened matlab file\n");
#endif
    int endchkval = 1;
    // Write precision type to mat file
    if (sizeof(usedtype) == sizeof(float)) {
      sprintf(tmpstr, "TmpChiMapOut.precision = 'single';\n");
    } else {
      sprintf(tmpstr, "TmpChiMapOut.precision = 'double';\n");
    }
    MPI_File_write(matfile, tmpstr, strlen(tmpstr),
                   MPI_CHAR, MPI_STATUS_IGNORE);

    // Write endian check value to binary file
    MPI_File_write(binfile, &endchkval, 1, MPI_INT, MPI_STATUS_IGNORE);

    // Write file open command to matlab
    sprintf(tmpstr, "TmpChiMapOut.fid = fopen('out.bin', 'r', 'b');\n");
    sprintf(tmpstr, "%sif fread(TmpChiMapOut.fid, 1, 'int32') ~= 1\n",
            tmpstr);
    sprintf(tmpstr, "%s\tfclose(TmpChiMapOut.fid);\n", tmpstr);
    sprintf(tmpstr, "%s\tTmpChiMapOut.fid = fopen('out.bin', 'r', 'l');\n",
            tmpstr);
    sprintf(tmpstr, "%s\tfread(TmpChiMapOut.fid, 1, 'int32');\n", tmpstr);
    sprintf(tmpstr, "%send\n", tmpstr);
    MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR,
                   MPI_STATUS_IGNORE);
  }

  initialised = true;
}

