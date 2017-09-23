#ifndef MAT_L__H
#define MAT_L__H
/* Copyright (C) 2017 IBM Corp.
 *  Licensed under the Apache License, Version 2.0 (the "License"); 
 * you may not use this file except in compliance with the License. 
 * You may obtain a copy of the License at
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, 
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either express or implied. See the License for the specific
 * language governing permissions and limitations under the License. 
 */
/*
 * mat_l - functions for matrices of type long
 */
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/mat_ZZ_p.h>
#include "vec_l.h"

typedef NTL::Mat<long> mat_l;

void mul(mat_l& X, const mat_l& A, const mat_l& B);
void mul(mat_l& X, const mat_l& A, const mat_l& B);
void mul(vec_l& x, const mat_l& A, const vec_l& b);
void mul(vec_l& x, const vec_l& a, const mat_l& B);
void mul(mat_l& X, const mat_l& A, long b);
void square(mat_l& X, const mat_l& A, bool transpose=false);
void ident(mat_l& X, long n);
bool isZero(const mat_l& M);

long writeToFile(mat_l& X, FILE* handle);

long readFromFile(mat_l& X, FILE* handle);

inline mat_l& operator-=(mat_l& X, const mat_l& Y)
{
  for (long i=0; i<X.NumRows(); i++) for (long j=0; j<X.NumCols(); j++) {
      X[i][j] -= Y[i][j];
    }
  return X;
}

inline mat_l& operator*=(mat_l& x, long a)
{
   mul(x, x, a);
   return x;
}

inline mat_l& operator*=(mat_l& X, const mat_l& Y)
{
  mul(X,X,Y);
  return X;
}

inline bool operator==(mat_l& X, const mat_l& Y)
{
  for (long i=0; i<X.NumRows(); i++) for (long j=0; j<X.NumCols(); j++) {
      if (!(X[i][j] == Y[i][j]))
        return false;
    }
  return true;
}

#endif // MATRIXEXT_H_INCLUDED
