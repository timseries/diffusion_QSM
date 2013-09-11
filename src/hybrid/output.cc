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
// THIS SOFTWARE IS PROVIDED BY Amanda Ng and Timothy Roberts ''AS IS'' AND ANY
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

/*! \file output.cc
  \brief Output class file.

  Implementation of the Output class.
*/

#include "hybrid/output.h"

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
Output::~Output() {}
void Output::Init(const ArgHandler &arghandler, int rank, int size) {
  bool flg;
  // TODO(timseries): replace sprintf with snprintf which is more secure.
  this->rank=rank;
  this->size=size;
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
  if (rank==0) printroot("   Output::Init : redirected stderr\n");
#endif
  // Open binary file
  if (rank == 0) {
    sprintf(tmpstr, "%s/out.bin", outdir);
    MPI_File_open(MPI_COMM_SELF, tmpstr, MPI_MODE_WRONLY | MPI_MODE_CREATE,
s                  MPI_INFO_NULL, &binfile);
#ifdef DEBUG
    if (rank==0) printroot("   Output::Init : opened binary file\n");
#endif
    MPI_File_set_size(binfile, 0);
    // Open matlab file on process 0
    sprintf(tmpstr, "%s/out.m", outdir);
    MPI_File_open(MPI_COMM_SELF, tmpstr, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &matfile);
    MPI_File_set_size(matfile, 0);
#ifdef DEBUG
    if (rank==0) printroot("   Output::Init : opened matlab file\n");
#endif
    int endchkval = 1;
    // Write precision type to mat file
    if (sizeof(Real) == sizeof(float)) {
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

void Output::LocalArray(int onproc, Real* array, int ndims, int* dims, const char* arrayname) {
  int			ndims0;
  int			*dims0;
  int			arraynamelength;
  char		*arrayname0;
  int			n;
  Real	*array0;

  if (rank != onproc && rank != 0) return;

  // if data is not on process 0, send data to process 0
  if (onproc != 0) {

    // send data
    if (rank == onproc) {
      n = 1;
      for (int i = 0; i < ndims; i++)
        n *= dims[i];

      arraynamelength = strlen(arrayname);

      MPI_Send(&ndims,1,MPI_INT,0,0,MPI_COMM_WORLD);
      MPI_Send(dims,ndims,MPI_INT,0,1,MPI_COMM_WORLD);
      MPI_Send(&arraynamelength,1,MPI_INT,0,2,MPI_COMM_WORLD);
      MPI_Send((char*)arrayname,arraynamelength,MPI_CHAR,0,3,MPI_COMM_WORLD);
      MPI_Send(array,n,MPI_Real,0, 4, MPI_COMM_WORLD);
    }

    // receive data
    if (rank == 0) {
      MPI_Recv(&ndims0,1,MPI_INT,onproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      dims0 = (int*) calloc(ndims0, sizeof(int));
      MPI_Recv(dims0,ndims0,MPI_INT,onproc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&arraynamelength,1,MPI_INT,onproc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      arrayname0 = (char*) calloc(arraynamelength, sizeof(char));
      MPI_Recv(arrayname0,arraynamelength,MPI_CHAR,onproc,3,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      n = 1;
      for (int i = 0; i < ndims0; i++)
        n *= dims0[i];
      array0 = (Real*) calloc(n, sizeof(Real));
      MPI_Recv(array0,n,MPI_Real,onproc,4,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
  // otherwise, set values and pointers as though process 0 had sent data to itself
  else {
    ndims0 = ndims;
    dims0 = dims;
    arrayname0 = (char*) arrayname;
    array0 = array;
    n = 1;
    for (int i = 0; i < ndims0; i++)
      n *= dims0[i];
  }

  // On process 0, write the data to file
  if (rank == 0) {

    // write array to binary file
    MPI_File_write(binfile, array0, n, MPI_Real, MPI_STATUS_IGNORE);

    // write read-in matlab code to matlab file
    if (ndims > 1) {
      sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname0, n, dims0[0]);
      for (int i = 1; i < ndims0; i++) {
        sprintf(tmpstr, "%s %d", tmpstr, dims0[i]);
      }
      sprintf(tmpstr, "%s]);\n", tmpstr);
    }
    else {
      sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname0, n);
    }

    MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);

    // if data was sent from another process, free the memory that was allocated
    if (onproc != 0) {
      free(dims0);
      free(arrayname0);
      free(array0);
    }

  }

  // THE FOLLOWING IS THE TRUE MPI IMPLEMENTATION. DUE TO A BUG ON AVOCA, THIS HAS BEEN REPLACED BY THE ABOVE
  /*
     if (rank != onproc) return;

     int n = 1;
     for (int i = 0; i < ndims; i++)
     n *= dims[i];

     MPI_File_write_shared(binfile, array, n, MPI_Real, MPI_STATUS_IGNORE);

     if (ndims > 1) {
     sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname, n, dims[0]);
     for (int i = 1; i < ndims; i++) {
     sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
     }
     sprintf(tmpstr, "%s]);\n", tmpstr);
     }
     else {
     sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname, n);
     }

     MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
  */
}

void Output::LocalArray(int onproc, bool* array, int ndims, int* dims, const char* arrayname) {
  int			ndims0;
  int			*dims0;
  int			arraynamelength;
  char		*arrayname0;
  int			n;
  bool		*array0;

  if (rank != onproc && rank != 0) return;

  // if data is not on process 0, send data to process 0
  if (onproc != 0) {

    // send data
    if (rank == onproc) {
      n = 1;
      for (int i = 0; i < ndims; i++)
        n *= dims[i];

      arraynamelength = strlen(arrayname);
      MPI_Send(&ndims, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(dims, ndims, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&arraynamelength, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
      MPI_Send((char*)arrayname, arraynamelength, MPI_CHAR, 0, 3, MPI_COMM_WORLD);
      MPI_Send(array, n, MPI_CHAR, 0, 4, MPI_COMM_WORLD);
    }

    // receive data
    if (rank == 0) {
      MPI_Recv(&ndims0, 1, MPI_INT, onproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      dims0 = (int*) calloc(ndims0, sizeof(int));
      MPI_Recv(dims0, ndims0, MPI_INT, onproc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&arraynamelength, 1, MPI_INT, onproc, 
               2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      arrayname0 = (char*) calloc(arraynamelength, sizeof(char));
      MPI_Recv(arrayname0, arraynamelength, MPI_CHAR,
               onproc, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      n = 1;
      for (int i = 0; i < ndims0; i++)
        n *= dims0[i];
      array0 = (bool*) calloc(n, sizeof(bool));
      MPI_Recv(array0, n, MPI_CHAR, onproc, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // otherwise, set values and pointers as though process 0 had sent data to itself
  else {
    ndims0 = ndims;
    dims0 = dims;
    arrayname0 = (char*) arrayname;
    array0 = array;
    n = 1;
    for (int i = 0; i < ndims0; i++)
      n *= dims0[i];
  }

  // On process 0, write the data to file
  if (rank == 0) {

    // write array to binary file
    MPI_File_write(binfile, array0, n, MPI_CHAR, MPI_STATUS_IGNORE);

    // write read-in matlab code to matlab file
    if (ndims > 1) {
      sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, 'uint8'), [%d", arrayname0, n, dims0[0]);
      for (int i = 1; i < ndims0; i++) {
        sprintf(tmpstr, "%s %d", tmpstr, dims0[i]);
      }
      sprintf(tmpstr, "%s]);\n", tmpstr);
    }
    else {
      sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, 'uint8');\n", arrayname0, n);
    }

    MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);

    // if data was sent from another process, free the memory that was allocated
    if (onproc != 0) {
      free(dims0);
      free(arrayname0);
      free(array0);
    }
  }
}

void Output::DistrArray(Real* array, int localsize, int ndims, int* dims, const char* arrayname) {
  int			n;
  Real	*array0 = NULL;
  MPI_Request req[2];
  // Send data to process 0
  if (rank>0){
    MPI_Isend(&localsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req[0]);
    MPI_Isend(array, localsize, MPI_Real, 0, 1, MPI_COMM_WORLD, &req[1]);
  }

  // On process 0, receive and write data to file.
  if (rank == 0) {

    n = 0; // current length of array0

    // Receive array length and data from each process in order, and write to binary file
    for (int p = 0; p < size; p++) {
      if (rank==0) printroot( "receiving 1 of 2 ...\n");
      if (rank==0) printroot("local size: %d\n", localsize);
      if (rank==0) printroot("mpi_int: %d\n", MPI_INT);
      if (rank==0) printroot("n: %d\n", n);
      if (rank==0) printroot("p: %d\n", p);
      if (rank==0) printroot("MPI_COMM_WORLD: %d\n", MPI_COMM_WORLD);
      if (p>0){
        MPI_Recv(&localsize, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        // realloc memory if the localsize has increased
        if (n < localsize) {
          if (rank==0) printroot( "reallocating  ...\n");
          array0 = (Real*) realloc(array0, localsize*sizeof(Real));
          n = localsize;
        }
        if (rank==0) printroot( "receiving 2 of 2 ...\n");

        MPI_Recv(array0, localsize, MPI_Real, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      else{
        array0=array;
      }
      if (rank==0) printroot( "writing 1 of 2 ...\n");

      MPI_File_write(binfile, array0, localsize, MPI_Real, MPI_STATUS_IGNORE);
      if (p==0){
        array0=NULL;
      }
    }
    // Write read-in matlab code to the matlab file
    int n = 1;
    for (int i = 0; i < ndims; i++)
      n *= dims[i];

    if (ndims > 1) {
      sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname, n, dims[0]);
      for (int i = 1; i < ndims; i++) {
        sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
      }
      sprintf(tmpstr, "%s]);\n", tmpstr);
    }
    else {
      sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname, n);
    }
    if (rank==0) printroot( "writing 2 of 2 ...\n");

    MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  // Free array
  if (array0 != NULL)	free(array0);

  // Wait for requests to complete
  if (rank==0) printroot( "waiting ...\n");
  if (rank>0){
    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
  }
  if (rank==0) printroot( "done waiting ...\n");
  // THE FOLLOWING IS THE TRUE MPI IMPLEMENTATION. DUE TO A BUG ON AVOCA, THIS HAS BEEN REPLACED BY THE ABOVE
  /*
    int n = 1;
    for (int i = 0; i < ndims; i++)
    n *= dims[i];

    MPI_File_write_ordered(binfile, array, localsize, MPI_Real, MPI_STATUS_IGNORE);

    if (rank == 0) {
    if (ndims > 1) {
    sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname, n, dims[0]);
    for (int i = 1; i < ndims; i++) {
    sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
    }
    sprintf(tmpstr, "%s]);\n", tmpstr);
    }
    else {
    sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname, n);
    }

    MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
    }
  */
}

void Output::Close() {

  MPI_Offset offset, disp;

  if (rank == 0) {
    sprintf(tmpstr, "fclose(TmpChiMapOut.fid);\n");
    sprintf(tmpstr, "%sclear TmpChiMapOut\n", tmpstr);
    //MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Close binary file
  //	MPI_File_get_position_shared(binfile, &offset);
  //	MPI_File_get_byte_offset(binfile, offset, &disp);
  //	MPI_File_set_size(binfile, disp);
  if (rank == 0)
    MPI_File_close(&binfile);

  // Close matlab file
  //	MPI_File_get_position_shared(matfile, &offset);
  //	MPI_File_get_byte_offset(matfile, offset, &disp);
  //	MPI_File_set_size(matfile, disp);
  if (rank == 0)
    MPI_File_close(&matfile);

  // Close error log file and remove if empty
  if (ftell(errfile) == 0) {
    fclose(errfile);
    sprintf(tmpstr, "%s/err%03d.txt", outdir, rank);
    remove(tmpstr);
  }
  else {
    fclose(errfile);
  }

  // Close files
  initialised = false;
  free(tmpstr);
}
