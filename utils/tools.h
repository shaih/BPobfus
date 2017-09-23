#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED
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
#include <sstream>
#include <thread>
#include <NTL/mat_lzz_p.h>
#include <NTL/pair.h>

// VJS: I defined this because std::to_string is C++11,
// which I would prefer avoiding
// SH: ifdefed it to use std::to_string if available
#if (__cplusplus<201103L)
template<class T>
std::string ToString(const T& x)
{
   std::stringstream ss;
   ss << x;
   return ss.str();
}
#else
#define ToString std::to_string
#endif


// Find the number of cores on the machine. This is part fo the C++11 standard,
// for older *Nix machines we use sysconf, for older Win32 we use GetSystemInfo
#if (__cplusplus<201103L)
#if defined(__unix__) || defined(__unix) || defined(unix) // Unix
#include <unistd.h>
inline long numCores() { return sysconf( _SC_NPROCESSORS_ONLN ); }
#elif defined(__MINGW32__) || defined(WIN32)              // Windows
#include <windows.h>
inline long numCores() { SYSTEM_INFO sysinfo;
                         GetSystemInfo( &sysinfo );
                         return sysinfo.dwNumberOfProcessors; }
#else                                                     // somethign else??
inline long numCores() { return 1; } // FIXME: anything else to do here?
#endif
#else // C++11
inline long numCores() { return std::thread::hardware_concurrency(); }
#endif


#if 1

#define GEXEC_RANGE_END(seq)\
NTL_GEXEC_RANGE_END

#define GEXEC_RANGE(seq,cat,first,last)\
  NTL_GEXEC_RANGE(seq,cat,first,last)

  #define EXEC_RANGE(cat,first,last)\
   NTL_EXEC_RANGE(cat,first,last)

   #define EXEC_RANGE_END \
NTL_EXEC_RANGE_END

#define EXEC_INDEX(nt,index)\
NTL_EXEC_INDEX(nt,index)

#define EXEC_INDEX_END \
NTL_EXEC_INDEX_END

#else

#define GEXEC_RANGE_END(seq)\
NTL_GEXEC_RANGE_END\
    if (!seq) {\
    fprintf(stderr,"end loop %s %d\n", __FILE__, __LINE__);}

#define GEXEC_RANGE(seq,cat,first,last)\
    if (!seq) {\
    fprintf(stderr,"start loop %s %d\n", __FILE__, __LINE__);}\
    NTL_GEXEC_RANGE(seq,cat,first,last)

#define EXEC_RANGE(cat,first,last)\
    fprintf(stderr,"start loop %s %d\n", __FILE__, __LINE__);\
    NTL_EXEC_RANGE(cat,first,last)

#define EXEC_RANGE_END \
NTL_EXEC_RANGE_END \
    fprintf(stderr,"end loop %s %d\n"  , __FILE__, __LINE__);

#define EXEC_INDEX(nt,index)\
fprintf(stderr,"start loop %s %d\n", __FILE__, __LINE__);\
NTL_EXEC_INDEX(nt,index)

#define EXEC_INDEX_END \
NTL_EXEC_INDEX_END \
fprintf(stderr,"end loop %s %d\n", __FILE__, __LINE__);\

#endif


// Matrix inversion modulo the current modulus
bool ppInvert(NTL::Mat<NTL::zz_p>& X, const NTL::Mat<NTL::zz_p>& A, long p, long r);


// Preconditions:
//   * MM and M have same dimensions
//   0 <= MM[i][j] < P
//   gcd(p, P) = 1, where p = zz_p::modulus()
//   p*P can be represented as a long
// Postcondition:
//   0 <= MM'[i][j] < P*p
//   MM'[i][j] = MM[i][j] mod P and MM'[i][j] = M[i][j] mod p
//
//  Here, MM' is the new value of MM and MM is the original value
void CRT(NTL::Mat<long>& MM, long P, const NTL::Mat<NTL::zz_p>& M);

// Returns a list of prime factors and their multiplicity,
// N = \prod_i (factors[i].a)^{factors[i].b}
void factorize(NTL::Vec< NTL::Pair<long,long> > &factors, long N);


// VJS: I changed it to an assignment of 0,
// rather than T(0).  Should be more portable and more efficient.
template<class T> void clear(NTL::Vec<T>& v)
{
  for (long i=0; i < v.length(); i++) v[i] = 0;
}

template<class T> void clear(NTL::Mat<T>& X)
{
  for (long i=0; i<X.NumRows(); i++)
    for (long j=0; j<X.NumCols(); j++) X[i][j] = 0;
}

void RandomFill(NTL::Mat<NTL::zz_p>& X);

// Convert to zz_p and invert, assumes that p is easy to factor
// Returns true on success, false on failure
bool invMod(NTL::Mat<NTL::zz_p>& X, const NTL::Mat<long>& Y);

// generates a random pair of nxn matrices with X*Y = I
// assumes zz_p::modulus() is easy to factor
void GenerateMatrixPair(NTL::Mat<NTL::zz_p>& X, NTL::Mat<NTL::zz_p>& Y, long n);

void seekPastChar(std::istream& str, int cc);

#endif // ifdef TOOLS_H_INCLUDED
