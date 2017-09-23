/* Copyright (C) 2017 IBM Corp.
 * Licensed under the Apache License, Version 2.0 (the "License"); 
 * you may not use this file except in compliance with the License. 
 * You may obtain a copy of the License at
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, 
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either express or implied. See the License for the specific
 * language governing permissions and limitations under the License. 
 */
#include "argmap.h"

// Code for parsing command line

bool parseArgs(int argc,  char *argv[], argmap_t& argmap)
{
  for (long i = 1; i < argc; i++) {
    char *x = argv[i];
    long j = 0;
    while (x[j] != '=' && x[j] != '\0') j++;
    if (x[j] == '\0') return false;
    std::string arg(x, j);
    if (argmap[arg] == NULL) return false;
    argmap[arg] = x+j+1;
  }

  return true;
}

bool doArgProcessing(std::string *value, const char *s)
{
  *value = std::string(s);
  return true;
}
