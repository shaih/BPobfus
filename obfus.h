#ifndef _OBFUS_H_
#define _OBFUS_H_
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
/**************************************************************************
 obfus.h - obfuscating (read-once) branching programs
 *
 * Usage pattern:
 *
 * The BP is given in the 'trans' argument, represented in an NTL matrix with
 * one row per step of the BP, and one column per symbol of the alphabet. Each
 * entry trans[i][j] is a transition matrix of the BP (represented as mat_l),
 * and we assume that all the matrices are square and have the same dimension.
 *
 * The obfuscated BP is written to the subdirectory 'dirName' of the current
 * directory, each matrix in a separate file. In the obfuscated program, each
 * transition matrix is represented by two randomized matrices, one on the
 * "main chain" and one on the "dummy chain". Each matrix is written in a file
 * whose name is <from>_<to>_<symbol>.dat
 *
 * Evaluating this BP on an input string is done via the function
 *
 *   bool funcVal = evalNFA(dirName, symbolString, params=NULL);
 *
 * This function expects to find the obfuscated branching program in the
 * subdirectory 'dirName' of the current directory, it evaluates it on the
 * given input string and returns a boolean true/false. If the params argument
 * is pecified then it is used to check that the BP is consistent with the
 * given paameters (else an exception is raised).
 *
 *******************************************************************************/

#include "GGH15.h"

// Obfuscate a (read one) BP using the given nodes
void obfuscateNFASavedNodes(const std::string& dirName,
			    const NTL::Mat< NTL::Mat<long> >& trans,
			    TDMatrixParams* params, bool pipes=false,
			    long nThreads=1);

// Obfuscate a (read one) BP on a new single0use mmaps instance
void initObfuscateNFA(const std::string& dirName,
                      const NTL::Mat< NTL::Mat<long> >& trans, long sec,
		      TDMatrixParams* params, bool bPipe=false, long threads=1);

// Evaluate an obfuscated (read once) BP on a given input string
bool evalNFA(const std::string& name, const Vec<long>& symbolString,
	     TDMatrixParams *params=NULL, long nThreads=1);

//pipedInit forks two threads, one "writer" thread and the "producer" thread, which can be further multi-threaded
bool pipedInit(const std::string& dirName, TDMatrixParams& p, long n,
	       long nThreads);

// The following two functions are low-level functions that are used
// by obfuscateNFA, they are externalized for debugging purposes

// Encode a transition matrix via a pair of matrices
void randomizeTransitions(GGH15ptxt& p1, GGH15ptxt& p2,
			  const NTL::Mat<long>& transition,
			  long i, long j, long nSteps, long dim);

// Create a new directory and write the parameters in that directory
long saveParamstoNewDir(const std::string& dirName, TDMatrixParams& p);

long makeRandomBP(std::string dirName, std::string name, long dim, long sig, long L);
//long makeRandomBP(std::string name, long dim, long sig, long L);

void writeNodesToFile(SynchronizedPipe<nodePair>& pipe);
typedef NTL::Pair<FILE*,GGH15encoding*> encPair;

// Obfuscate a (read one) BP using the given mmaps instance
inline void obfuscateNFA(const std::string& name,
                  const NTL::Mat< NTL::Mat<long> >& trans, GGH15mmaps &maps)
{
  throw std::logic_error("obfuscateNFA is a deprecated function");
}
#endif //_OBFUS_H_
