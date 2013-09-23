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

/*! \file arghandler.h
    \brief ArgHandler class definitions file.

*/

#ifndef INCLUDE_HYBRID_ARGHANDLER_H_
#define INCLUDE_HYBRID_ARGHANDLER_H_

#include "hybrid/common.h"

/**
* ArHandler class. Used to parse the command line arguments.
*/
class ArgHandler {
 public:
/**
* ArgHandler constructor.
*/
  ArgHandler();
/**
* ArgHandler destructor.
*/
  virtual ~ArgHandler();
/**
* Class initialization method. 
* @param argc Number of command line arguments
* @param args The command line arguments
*/
  void Init(int argc, char** args);
/**
* Getter method (string version).
* @param *option The command line option we wish to fetch. 
* @param *&str A pointer to a string to which the command line option value is copied.
* @return True if the command line option is found. False otherwise.
*/
  bool GetArg(const char* option, char *&str) const ;
/**
* Getter method (Real version).
* @param *option The command line option we wish to fetch. 
* @param &value A pointer to a Real to which the command line option value is copied.
* @return True if the command line option is found. False otherwise.
*/
  bool GetArg(const char* option, Real &value) const ;
/**
* Getter method (int version).
* @param *option The command line option we wish to fetch. 
* @param &value A pointer to an int to which the command line option value is copied.
* @return True if the command line option is found. False otherwise.
*/
  bool GetArg(const char* option, int &value) const;
 private:
  int argc; ///< Number of command line arguments
  char** args; ///< The command line arguments
};
#endif  // INCLUDE_HYBRID_ARGHANDLER_H_
