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

/*! \file arghandler.cc
  \brief ArgHandler class file.

  Implementation of the ArgHandler class.
*/

#include "include/arghandler.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

ArgHandler::Init(int argc, char** args) {
  this->argc = argc;
  this->args = args;
}

ArgHandler::GetArg(const char* option, char *&str) {
  for (int i = 0; i < argc; i++) {
    if (strcmp(args[i], option) == 0) {
      if (i == argc-1) {
        return false;
      }
      str = reinterpret_cast<char*>(calloc(strlen(args[i+1])+1, 
                                           SIZEOF_CHAR));
      // TODO(timseries): replace strcpy with snprintf which is more secure.
      strcpy(str, args[i+1]);
      return true;
    }
  }
  return false;
}

ArgHandler::GetArg(const char* option, const usedtype &value) {
  for (int i = 0; i < argc; i++) {
    if (strcmp(args[i], option) == 0) {
      if (i == argc-1) {
        return false;
      }
      value = atof(args[i+1]);
      return true;
    }
  }
  return false;
}

ArgHandler::GetArg(const char* option, int &value) {
  for (int i = 0; i < argc; i++) {
    if (strcmp(args[i], option) == 0) {
      if (i == argc-1) {
        return false;
      }
      value = atoi(args[i+1]);
      return true;
    }
  }
  return false;
}
