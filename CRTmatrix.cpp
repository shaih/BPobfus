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
/***********************************************************************
CRTmatrix: implements A matrix in CRT representation
************************************************************************/

#include <vector>
#include <cassert>
#include <cstdlib>
//#include <mpfr.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <NTL/FFT.h>
#include <NTL/SmartPtr.h>

NTL_CLIENT
#include "utils/tools.h"
#include "DGaussSampler.h"
#include "CRTmatrix.h"
#include "TDMatrixParams.h"

//#define DEBUGPRINT

CRTmatrix& CRTmatrix::leftMultBy(const CRTmatrix& M) // A := M*A
{
    FHE_TIMER_START;
    assert(length()==params->kFactors);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    // multi-threaded implementation, threads handle different moduli

    EXEC_RANGE(length(), first, last)
    //EXEC_RANGE(length(), first, last)
    for (long i=first; i<last; i++)
    {
        params->zzp_context[i].restore();
        mat_zz_p& matMod = (*this)[i];
        const mat_zz_p& Mmod = M[i];
        mul(matMod, Mmod, matMod);
    }
    //EXEC_RANGE_END
    EXEC_RANGE_END
    return *this;
}

CRTmatrix& CRTmatrix::leftMultBy(const mat_l& M) // A := M*A
{
    FHE_TIMER_START;
    assert(length()==params->kFactors);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    // multi-threaded implementation, threads handle different moduli
    EXEC_RANGE(length(), first, last)
    for (long i=first; i<last; i++)
    {
        params->zzp_context[i].restore();
        mat_zz_p& matMod = (*this)[i];
        mat_zz_p Mmod = conv<mat_zz_p>(M);
        mul(matMod, Mmod, matMod);
    }
    EXEC_RANGE_END

    return *this;
}

CRTmatrix& CRTmatrix::leftMultBy(const vec_l& u) // A := u*A
{
    FHE_TIMER_START;
    assert(length()==params->kFactors);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    // multi-threaded implementation, threads handle different moduli
    EXEC_RANGE(length(), first, last)

    for (long i=first; i<last; i++)
    {
        params->zzp_context[i].restore();
        mat_zz_p& matMod = (*this)[i]; // A
        const vec_zz_p umod = conv<vec_zz_p>(u);
        vec_zz_p vTmp;
        mul(vTmp, umod, matMod); // tmp = u*A

        // copy vTmp back to matMod
        matMod.SetDims(1, vTmp.length());
        matMod[0] = vTmp;
    }
    EXEC_RANGE_END

    return *this;
}

// Convert Mat<long> to CRT format, we assume that params are already set
CRTmatrix& CRTmatrix::operator=(const mat_l& M) // A := M
{
    FHE_TIMER_START;
    assert(params != NULL);

    NTL::zz_pPush ppush; // backup NTL's current modulus
    SetLength(params->kFactors);
    for (long i=0; i<length(); i++)
    {
        params->zzp_context[i].restore();
        conv((*this)[i], M);
    }

    return *this;
}


// Default assignment should work: NTL now guarantees
// that assignment and default constructors work, even
// "out of context"
CRTmatrix& CRTmatrix::operator=(const CRTmatrix& other) // A := M
{
    FHE_TIMER_START;
    params = other.params;
    NTL::zz_pPush ppush; // backup NTL's current modulus
    SetLength(other.length());
    for (long i=0; i<params->kFactors; i++)
    {
        params->zzp_context[i].restore();
        (*this)[i] = other[i];
    }
    return *this;
}

CRTmatrix& CRTmatrix::operator*=(const mat_l& M) // A := M*A
{
    FHE_TIMER_START;
    assert(length()==params->kFactors);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    // multi-threaded implementation, threads handle different moduli
    EXEC_RANGE(length(), first, last)
    mat_zz_p Mmod;
    for (long i=first; i<last; i++)
    {
        params->zzp_context[i].restore();
        mat_zz_p& Amod = (*this)[i];
        Mmod = conv<mat_zz_p>(M);
        mul(Amod, Amod, Mmod);
    }
    EXEC_RANGE_END

    return *this;
}

CRTmatrix& CRTmatrix::operator*=(const CRTmatrix& M) // A := A*M
{
    FHE_TIMER_START;
    assert(length()==params->kFactors
           && M.length()==M.params->kFactors
           && params==M.params);
    NTL::zz_pPush ppush; // dibackup NTL's current modulus

    // multi-threaded implementation, threads handle different moduli
    EXEC_RANGE(length(), first, last)
    for (long i=first; i<last; i++)
    {
        params->zzp_context[i].restore();
        mat_zz_p& Amod = (*this)[i];
        const mat_zz_p& Mmod = M[i];
        mul(Amod, Amod, Mmod);
    }
    EXEC_RANGE_END

    return *this;
}

CRTmatrix& CRTmatrix::operator+=(const mat_l& M) // A := A+M
{
    FHE_TIMER_START;
    assert(length()==params->kFactors);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    // multi-threaded implementation, threads handle different moduli
    EXEC_RANGE(length(), first, last)
    //for (long i=0; i<length(); i++)
    for (long i=first; i<last; i++)
    {
        mat_zz_p& matMod = (*this)[i];
        assert(matMod.NumRows()==M.NumRows() && matMod.NumCols()==M.NumCols());

        params->zzp_context[i].restore();
        for (long j=0; j<matMod.NumRows(); j++)
            for (long k=0; k<matMod.NumCols(); k++)
                matMod[j][k] += M[j][k];
    }
    EXEC_RANGE_END
    return *this;
}
CRTmatrix& CRTmatrix::operator-=(const mat_l& M) // A := A+M
{
    FHE_TIMER_START;
    assert(length()==params->kFactors);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    EXEC_RANGE(length(), first, last)
    //for (long i=0; i<length(); i++)
    for (long i=first; i<last; i++)
    {
        mat_zz_p& matMod = (*this)[i];
        assert(matMod.NumRows()==M.NumRows() && matMod.NumCols()==M.NumCols());

        params->zzp_context[i].restore();
        for (long j=0; j<matMod.NumRows(); j++)
            for (long k=0; k<matMod.NumCols(); k++)
                matMod[j][k] -= M[j][k];
    }
    EXEC_RANGE_END
    return *this;
}

CRTmatrix& CRTmatrix::operator+=(const CRTmatrix& M) // A := A+M
{
    FHE_TIMER_START;
    assert(length()==params->kFactors
           && M.length()==M.params->kFactors
           && params==M.params);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    EXEC_RANGE(length(), first, last)
    //for (long i=0; i<length(); i++)
    for (long i=first; i<last; i++)
    {
        mat_zz_p& Amod = (*this)[i];
        const mat_zz_p& Mmod = M[i];

        assert(Amod.NumRows()==Mmod.NumRows() && Amod.NumCols()==Mmod.NumCols());
        params->zzp_context[i].restore();
        Amod += Mmod;
    }
    EXEC_RANGE_END
    return *this;
}

CRTmatrix& CRTmatrix::operator-=(const CRTmatrix& M) // A := A+M
{
    FHE_TIMER_START;
    assert(length()==params->kFactors
           && M.length()==M.params->kFactors
           && params==M.params);
    NTL::zz_pPush ppush; // backup NTL's current modulus

    EXEC_RANGE(length(), first, last)
    //for (long i=0; i<length(); i++)
    for (long i=first; i<last; i++)
    {
        mat_zz_p& Amod = (*this)[i];
        const mat_zz_p& Mmod = M[i];

        assert(Amod.NumRows()==Mmod.NumRows() && Amod.NumCols()==Mmod.NumCols());
        params->zzp_context[i].restore();
        Amod -= Mmod;
    }
    EXEC_RANGE_END
    return *this;
}

// Get a CRT column: crtCol[i] is the j'th column of this[i]
void CRTmatrix::getColumn(Vec<vec_zz_p>& crtCol, long j)
{
    FHE_TIMER_START;
    crtCol.SetLength(length());
    for (long n=0; n<length(); n++)
    {
        vec_zz_p& colMod = crtCol[n];
        mat_zz_p& matMod = (*this)[n];

        colMod.SetLength(matMod.NumRows());
        // NOTE: this is OK, but beware that this would
        // not work for ZZ_p's without restoring context

        for (long i=0; i<matMod.NumRows(); i++)
            colMod[i] = matMod[i][j];
    }
}

//checks if all values of the matrix are smaller than the modolus
bool CRTmatrix::isSmall(long bitGap) const
{
    FHE_TIMER_START;
    mat_ZZ M; // reconstruct the matrix in ZZ format
    ZZ q;     // the modulus

    convert(M, q, *this); // cnovert to ZZ format

    q >>= bitGap;         // q / 2^{bitGap}

    for (long i=0; i<M.NumRows(); i++) for (long j=0; j<M.NumCols(); j++)
        {
            if (abs(M[i][j]) > q)
            {
                return false;
            }
        }
    return true;
}

// Choose a random n-by-m matrix in CRT representation
void CRTmatrix::randomFill(long n, long m)
{
    FHE_TIMER_START;

    const TDMatrixParams& prms = *params;
    SetLength(prms.kFactors);

    zz_pPush ppush; // backup NTL's current modulus

    EXEC_RANGE(length(), first, last)
    //for (long i=0; i< prms.kFactors; i++)
    for (long i=first; i< last; i++)
    {

        prms.zzp_context[i].restore();
        (*this)[i].SetDims(n, m);
        RandomFill((*this)[i]);
    }
    EXEC_RANGE_END
}

// Invert M modulo the product of the factors
bool CRTmatrix::invert(const mat_l& M)
{
    FHE_TIMER_START;
    SetLength(params->kFactors);

    zz_pPush ppush; // backup NTL's current modulus

    std::atomic<bool> result(true);

    EXEC_RANGE(params->kFactors, first, last)

    // Invert M wrt each factor separately
    //for (long i=0; i< params->kFactors; i++)
     for (long i=first; i< last; i++)
    {
        params->zzp_context[i].restore();
        // NOTE: Alternatively, can call relaxed_inv directly now
        if (!invMod((*this)[i], M))   // not invertible
        {
            SetLength(0);
            result = false;
            //return false;
        }
    }
    EXEC_RANGE_END

    return result;
    //return true;
}


void CRTmatrix::identity(long n)
{
    FHE_TIMER_START;
    zz_pPush ppush; // backup NTL's current modulus
    SetLength(params->kFactors);

    // Invert M wrt each factor separately
    for (long i=0; i< params->kFactors; i++)
    {
        params->zzp_context[i].restore();
        ident( (*this)[i], n ); // set as the identity mod f_i
    }
}

bool operator==(const CRTmatrix& A, const CRTmatrix& B)
{
    FHE_TIMER_START;
    if (A.params != B.params) return false;

    zz_pPush push; // backup the NTL current modulus
    for (long i=0; i < A.params->kFactors; i++)
    {
        (A.params->zzp_context)[i].restore();
#ifdef DEBUGPRINT
        cout << "A[i]=" << A[i] << ", B[i]=" << B[i] << endl;
#endif // DEBUGPRINT
        if (A[i] != B[i]) return false;
    }
    return true;
}

//invert a CRT matrix
bool CRTmatrix::invert(const CRTmatrix& A)
{
    FHE_TIMER_START;
    params = A.params;
    SetLength(params->kFactors);

    zz_pPush push; // backup the NTL current modulus

    for (long i=0; i < params->kFactors; i++)
    {
        const mat_zz_p& Amod = A[i];
        params->zzp_context[i].restore();

        zz_p d;
        //    if (isInvertible(Amod, params->factors[i]))
        //inv((*this)[i],Amod); // raises error if Amod is singular
        relaxed_inv(d,(*this)[i],Amod); //raises error if Amod is singular

        if (IsZero(d))
            return false; //non invertible
    }
    return true;
}


// reconstruct the matrix in ZZ fromat, returning in q the modulus
void convert(mat_ZZ& to, ZZ& q, const CRTmatrix& from)
{
    FHE_TIMER_START;
    q = to_ZZ(1L); // initalize to 1
    if (from.length() <=0)
    {
        to.kill(); // set as 0-by-0 matrix
        return;
    }

    to.SetDims(from[0].NumRows(), from[0].NumCols());
    clear(to);     // initialize to zero

    NTL::zz_pPush ppush; // backup NTL's current modulus
    for (long k=0; k<from.length(); k++)   // go over the factors one at a time
    {
        from.params->zzp_context[k].restore();
        NTL::CRT(to, q, from[k]); // incremental CRT
    }
}

// reconstruct the matrix in CRT fromat
// FIXME: can be optimized for many factors
void convert(CRTmatrix& to, const mat_ZZ& from, const TDMatrixParams& prms)
{
    FHE_TIMER_START;
    to.params = (TDMatrixParams*) &prms;
    to.SetLength(prms.kFactors);
    if (prms.kFactors==0) return; // nothing to do

    NTL::zz_pPush ppush; // backup NTL's current modulus

    EXEC_RANGE(to.length(), first, last)

    for (long k=first; k<last; k++)   // go over the factors one at a time
    {
        prms.zzp_context[k].restore();
        conv(to[k], from); // compute M mod current factor
    }

    EXEC_RANGE_END
}


void generateMultiPair(CRTmatrix& A, CRTmatrix& Ainv, long n)
{
    FHE_TIMER_START;

    NTL::Vec<CRTmatrix> Amats(INIT_SIZE, 12, A);
    NTL::Vec<CRTmatrix> Ainvmats(INIT_SIZE, 12, Ainv);

    EXEC_RANGE(12, first, last)
    for (long i=first; i<last; i++)
    {

        CRTmatrix &Apair = Amats[i];
        CRTmatrix &AinvPair = Ainvmats[i];
        generateMatrixPair(Apair, AinvPair, n);
    }

    EXEC_RANGE_END
}


void generateMatrixPair(CRTmatrix& A, CRTmatrix& Ainv, long n)
{
    FHE_TIMER_START;
    const TDMatrixParams& prms = *(A.params);

    A.SetLength(prms.kFactors);
    Ainv.SetLength(prms.kFactors);

    zz_pPush push; // backup the NTL current modulus

    // It's possible to do parallilizing at this level
    // need to call SetNumThreads(nt) in main.

    // NOTE: the zz_pPush in the "main" thread is sufficient,
    // as the "worker" threads do not (and should not) assume
    // any particular contextual state

#ifdef DEBUGPRINT
        cerr << "*** " << prms.kFactors << " " << AvailableThreads() << "\n";
#endif // DEBUGPRINT

    EXEC_RANGE(prms.kFactors, first, last)
    for (long i=first; i<last; i++)
    {

        mat_zz_p& Amod = A[i];
        mat_zz_p& AinvMod = Ainv[i];
        prms.zzp_context[i].restore(); // set i'th factor as NTL's current modulus
        GenerateMatrixPair(Amod, AinvMod, n);
    }
    EXEC_RANGE_END

    FHE_TIMER_STOP;
};

// binary I/O
#undef CRTMAT_CONVtoZZ
long CRTmatrix::writeToFile(FILE* handle) const
{
    FHE_TIMER_START;

    // write the parameters first
    long count = params->writeToFile(handle);
    assert (this->length() == params->kFactors);

#ifndef CRTMAT_CONVtoZZ
    long n,m;
    if (this->length()<=0)
        n = m = 0;
    else
    {
        n = (*this)[0].NumRows();
        m = (*this)[0].NumCols();
    }
    count += fwrite(&n, sizeof(n), 1, handle); // how many rows
    count += fwrite(&m, sizeof(m), 1, handle); // how many cols

    unsigned char buf[params->e];
    NTL::zz_pPush ppush; // backup NTL's current modulus
    for (long f=0; f<this->length(); f++)   // write one factor at a time
    {
        params->zzp_context[f].restore();
        const Mat<zz_p>& mat = (*this)[f];
        for (long i=0; i<n; i++) for (long j=0; j<m; j++)   // write every entry
            {
                for (long ie=0; ie < params->e; ie++) // write every byte
                    buf[ie] = (rep(mat[i][j]) >> (ie*8)) & 0xff; // extract one byte
                count += fwrite(buf, params->e, 1, handle); // do the actual writing
            }
#ifdef DEBUGPRINT
        //cout << "*this[f]="  << (*this)[f] << endl;
#endif // DEBUGPRINT
    }


#else
    ZZ q;
    mat_ZZ zM;
    convert(zM, q, *this); // convert to ZZ representation, then write it

    long qBytes = NumBytes(q); // # of bytes to represent q
    count += fwrite(&qBytes, sizeof(qBytes),1, handle); // how many bytes per int

    long n = zM.NumRows();
    long m = zM.NumCols();
    count += fwrite(&n, sizeof(n), 1, handle); // how many rows
    count += fwrite(&m, sizeof(m), 1, handle); // how many cols

    unsigned char buf[qBytes];
    for (long i=0; i<n; i++)
        for (long j=0; j<m; j++)
        {
            if (sign(zM[i][j])<0) zM[i][j] += q; // map to interval [0,q-1]
            BytesFromZZ(buf, zM[i][j], qBytes);      // get binary representation
            count += fwrite(buf, qBytes, 1, handle); // write it to file
        }
#endif
    return count;
}

long CRTmatrix::readFromFile(FILE* handle, TDMatrixParams* prmBuf)
{
    FHE_TIMER_START;

#ifdef DEBUGPRINT
    cout << "reading" << endl;
#endif

    assert(params != NULL || prmBuf != NULL); // some pointer must be provided
    long count;

    // If buffer is given, make params point to it and don't overwrite from iput
    if (prmBuf != NULL)
    {
        TDMatrixParams p;
        count = p.readFromFile(handle);
        assert(p == *prmBuf); // sanity check
        params = prmBuf;      // point to given params
    }
    else
        count = params->readFromFile(handle); // overwrite params from input

#ifndef CRTMAT_CONVtoZZ
    long n,m;
    count += fread(&n, sizeof(n), 1, handle); // how many rows
    count += fread(&m, sizeof(m), 1, handle); // how many column
    if (count == 0) return count; //nothing is read

    this->SetLength(params->kFactors);
    if (params->kFactors==0) return count; // nothing to do

    unsigned char buf[params->e];
    NTL::zz_pPush ppush; // backup NTL's current modulus
    for (long f=0; f<this->length(); f++)   // write one factor at a time
    {
        params->zzp_context[f].restore();
        (*this)[f].SetDims(n,m);
        for (long i=0; i<n; i++) for (long j=0; j<m; j++)   // write every entry
            {
                count += fread(buf, params->e, 1, handle); // read from file



                long zp = buf[(params->e)-1];  // assemble the integer, byte by byte
                //for (long ie=1; ie < params->e; ie++)
                for (long ie=((params->e)-2); ie >=0 ; ie--)
                {
                    zp <<= 8;
                    zp += buf[ie];
                }
                (*this)[f][i][j].LoopHole() = zp; // avoid calling rem()
            }
    }
#else
    long qBytes, n, m;
    count += fread(&qBytes, sizeof(qBytes), 1, handle); // how many bytes per int
    count += fread(&n, sizeof(n), 1, handle); // how many rows
    count += fread(&m, sizeof(m), 1, handle); // how many column

    if (count == 0) return count; //nothing is read

    unsigned char buf[qBytes];
    mat_ZZ zM(INIT_SIZE, n, m);
    for (long i=0; i<n; i++) for (long j=0; j<m; j++)
        {
            count += fread(buf, qBytes, 1, handle); // read from file
            ZZFromBytes(zM[i][j], buf, qBytes);     // make a ZZ object
        }

    convert(*this, zM, *params); // convert to CRT representation
#endif
    return count;
}

