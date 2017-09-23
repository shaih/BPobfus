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
/****************************************************************************
mat_l: a module for handling matrices of type long
Used mainly in TDMatrix.cpp as trapDoorMatR is of type mat_l
*****************************************************************************/
#include <cassert>
#include <stdexcept>

#include <NTL/matrix.h>
#include "utils/tools.h"
#include "mat_l.h"
#include <cstdio>
#include <cstdlib>

NTL_CLIENT

//#define DEBUG

//X = I * n, set the dimension to n, clear and then set to value
void ident(mat_l& X, long n)
{
    X.SetDims(n,n);
    clear(X);
    for (long i = 0; i < n; i++)
        X[i][i] = 1L;
}

bool isZero(const mat_l& M)
{
  for (long i=0; i<M.NumRows(); i++) for (long j=0; j<M.NumCols(); j++) {
      if (M[i][j] != 0)
	return false;
    }
  return true;
}

/* square - compute X := A*A^t or X := A^t*A, depending on transpose flag.
   Does not reset the matrix X dimension as long as it is large enough
*/
void square(mat_l& X, const mat_l& A, bool transpose)
{
    long n = A.NumRows();
    long m = A.NumCols();

    // If X is not large enough, set its dimension. If X is large enough then
    // its dimension is not reset, and only the top-left n-by-m is affected.

    if (transpose) swap(n,m);
    if (X.NumRows() < n || X.NumCols() < m)
        X.SetDims(n,n);

    for (long i=0; i<n; i++)
    {
        for (long j=0; j<n; j++)
        {
            X[i][j] = 0;
            if (transpose)
            {
                for (long k=0; k<m; k++)
                    X[i][j] += A[k][i] * A[k][j];
            }
            else
            {
                for (long k=0; k<m; k++)
                    X[i][j] += A[i][k] * A[j][k];
            }
        }
    }
}



//mul - does not reset the matrix dimension as long as it is large enough
//multiplies a vector and a matrix

void mul(vec_l& x, const vec_l& a, const mat_l& B)
{
    long l = a.length();
    long m = B.NumCols();

    if (l != B.NumRows())
    {
        cerr << "mul(vec_l, mat_l) dimension mismatch, vec["<< l
             << "] x mat[" << B.NumRows() << "][" << B.NumCols() << "]\n";
        exit(0);
    }

    int oRow,iIndex, val;

    x.SetLength(l);
    for (oRow = 0; oRow < m; oRow++)
    {
        val=0;
        for (iIndex = 0; iIndex<l; iIndex++)
        {
            val += (a.at(iIndex) * B[iIndex][oRow]);
        }
        x[oRow] = val;
    }
}

//compute x = A*b, A is a matrix, x and b are vectors,
//function sets length of x to be equal to A.NumRows()
void mul(vec_l& x, const mat_l& A, const vec_l& b)
{
    int nRows, nCols, nVecSize;
    int i,j;
    long val;

    nCols = A.NumCols();
    nVecSize = b.length();
    nRows = A.NumRows();
    x.SetLength(nRows);

    if (nCols != nVecSize)
    {
        throw std::logic_error("Matrix dimensions do not match");
        return;
    }

    for (i = 0; i < nRows; i++)
    {
        val=0;
        for (j = 0; j < nVecSize; j++)
            val+= A[i][j] * b.at(j);
        x[i] = val;
    }
}

//matrix multiplication by a number
void mul(mat_l& X, const mat_l& A, long b)
{
    long n = A.NumRows();
    long m = A.NumCols();

    X.SetDims(n, m);
    long i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            X[i][j] = A[i][j] * b;
}


// matrix multiplication, the most naive implementation
void mul(mat_l& X, const mat_l& A, const mat_l& B)
{
    assert(A.NumCols() == B.NumRows());
    mat_l tmp(INIT_SIZE, A.NumRows(), B.NumCols());
    for (long i=0; i<tmp.NumRows(); i++) for (long j=0; j<tmp.NumCols(); j++)
        {
            tmp[i][j] = 0;
            for (long k=0; k<A.NumCols(); k++)
                tmp[i][j] += A[i][k] * B[k][j];
        }
    X = tmp;
}

//Reads and sets the size of the matrix and its values from an open file. The function gets the handle of the file to read from.
long readFromFile(mat_l& X, FILE* handle)
{
    long numRows,numCols;
    long count = fread(&numRows,sizeof(long),1,handle);
    count += fread(&numCols,sizeof(long),1,handle);
    X.SetDims(numRows,numCols);
    for (long i = 0; i < numRows; i++)
        count+= fread(X[i].elts(), numCols*sizeof(long), 1, handle);
    return count;
}

//Writes the size of the matrix and its values to an open file. Gets the handle of the file to write to.
long writeToFile(mat_l& X, FILE* handle)
{
    long numRows = X.NumRows();
    long numCols = X.NumCols();
    long count = fwrite(&numRows,sizeof(long),1,handle);
    count += fwrite(&numCols,sizeof(long),1,handle);
    for (long i = 0; i < numRows; i++)
        count+= fwrite(X[i].elts(), numCols*sizeof(long), 1, handle);
    return count;
}
