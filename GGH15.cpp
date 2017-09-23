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
/*****************************************************************************
 * GGH15.h - an implementation of the graph-based multilinear maps from
 *           [Gentry-Gorbunov-Halevi, TCC 2015]
 *****************************************************************************/
#include <stdexcept>
#include <unistd.h>
#include <sys/stat.h>

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>

#include "mat_l.h"
#include "GGH15.h"
#include "utils/tools.h"

bool bGGH15threading=1;

//#define DEBUGPRINT
//#define PRINTDOT
#ifdef PRINTDOT
#define printDot cerr << "." << std::flush
#else
#define printDot
#endif

// A lower bound on |q|/2^10 for an L-edge chain. We need to bound the size
// of a product of L matrices of dim m-by-m, where L-1 have entries drawn
// with parameter sqrt{sigmaX) and one has entries with parameter sec.
// The formula we use is q/2^10 > largest-singular-val(product), which
// translates to
//    q/2^10 > (2*sqrt(m)*sqrt(sigmaX))^{L-1} * lambda*sqrt(m+n)
//      \approx lambda * sigmaX^{(L-1)/2} * m^{L/2} * 2^L
long qBound(long sigmaX, long m, long L)
{
  // cerr << "@ L="<<L<<", sigmaX="<<sigmaX
  //      <<",log_2(sqrt(sigmaX))="<<(log(sigmaX)/(log(2)*2))
  //      <<", log_2(m)="<<(log(m)/log(2)) << endl;

  double logQ = ((L-1)*log(sigmaX) + L*log(m))/2.0;
  return ceil(logQ/log(2.0)) + L + 16;
}

// Computing basic parameters for length-'diam' MMAPs. In particular,
// assuming that "fresh" sampled encodings use "varince" sigmaX, we need
// to ensure that
//  (a) log q > qBound(sigmaX, m, L, sec)      (qBound from above)
//  (b)     m > mBound(log q, errBits=7, sec)  (mBound from TDMatrixParams)
//
// Since the size of q is determined by the number of small factors, i.e.,
// q = \prod_{i=0}^{k-1} p_i^e, we just try to find the smallest k for which
// the above constraints are satisfied.
//
// In addition, we also need to ensure that m >= n*k*e + mBar, where
// n = nStates + ceil(sqrt(sec)/2) and mBar = mbarBound(log q, n, sec)
// (from TDMatrixParams).

GGH15basicParams::GGH15basicParams(COMP_PARAMS_TYPE comp, //=COMP_GGH15PRMS
                                   long nStates, long L, long sec, long ee)
{
    FHE_TIMER_START;

    e = ee;
    n = nStates + ceil(sqrt(sec)/2);

    // k will never be smaller than L
    for (k=L; k<TDMATRIX_NUM_SMALL_FACTORS; k++) {
      long qBits = get_qBits(k, e);
      m = mBound(qBits, /*errBits=*/7, sec);

      long mBar = mbarBound(qBits, n, sec);
      long wLen = n*k*e;
      if (m < mBar+wLen) m = mBar+wLen;
      else               mBar = m-wLen;

      long sigmaX = getsigmaXVal(wLen, mBar, /*maxFactor=*/181);
      long qbound = qBound(sigmaX,m,L);

      if (qBits >= qbound) // all constraints are satisfied
	break;

      // otherwise we try a larger value of k
    }
    if (k >= TDMATRIX_NUM_SMALL_FACTORS) { // could not find a solution
      throw std::logic_error("Cannot set GGH15 parameters");
    }
}


// Initialize a single node, choose A w/ trapdoor and T, Tinv
void GGH15node::initGGH15node(TDMatrixParams& prms, long idx)
{
    #ifdef DEBUGPRINT
        fprintf(stderr,"initGGH15node %s %d\n", __FILE__, __LINE__);
    #endif // DEBUGPRINT

    FHE_TIMER_START;
    index = idx;
    params = &prms;

    A.initTDmatrix(prms); // Choose matrix A with trapdoor

    #ifdef DEBUGPRINT
        fprintf(stderr,"initGGH15node %s %d\n", __FILE__, __LINE__);
    #endif // DEBUGPRINT

    // Choose a small n-by-n matrix P, set Pinv as its inverse mod q
    Pinv.params = &prms;
    do {
      setSmall(P, prms.n, prms.n, /*howSmall=*/prms.maxFactor);
    #ifdef DEBUGPRINT
        fprintf(stderr,"initGGH15node %s %d\n", __FILE__, __LINE__);
    #endif // DEBUGPRINT
    }
    while (!Pinv.invert(P)); // Repeat until you get an invertible P mod q

    #ifdef DEBUGPRINT
        fprintf(stderr,"initGGH15node %s %d\n", __FILE__, __LINE__);
    #endif // DEBUGPRINT

    // choose T, Tinv as random inverses of each other mod q
    T.params = Tinv.params = &prms;
    generateMatrixPair(T, Tinv, prms.m);
}

// Initilize an instance with n nodes
GGH15mmaps::GGH15mmaps(TDMatrixParams& p, long n)
{
    FHE_TIMER_START;
    assert (n>0);

    params = &p;

    nodes.SetLength(n-2);

    //the flag allows parallelizing when it is set to false
    NTL_GEXEC_RANGE(!bGGH15threading, n-2, first, last)

    for (long i=first; i< last; i++)
    {
#ifdef DEBUGPRINT
        cout << "initializing node:" << i << endl;
#if (defined(__unix__) || defined(__unix) || defined(unix))
        struct rusage rusage;
        getrusage( RUSAGE_SELF, &rusage );
        cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
#endif // DEBUGPRINT

        nodes[i].initGGH15node(p, i+1);

#ifdef DEBUGPRINT
        cout << "after initializing node:" << i << endl;
#endif // DEBUGPRINT
    }
    NTL_GEXEC_RANGE_END

#ifdef DEBUGPRINT
    cout << "before randomfill:"  << endl;
#endif // DEBUGPRINT

    // The right bracket is a random non-small n-by-1 matrix
    randomFill(An, /*rows=*/p.n, /*cols=*/1, p);

#ifdef DEBUGPRINT
    cout << "after randomfill:"  << endl;
#endif // DEBUGPRINT
}

long GGH15mmaps::writeToFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = params->writeToFile(handle);
    count += An.writeToFile(handle);
    long nNodes = numNodes()-2;
    count += fwrite(&nNodes,sizeof(long),1,handle);
    for (long i = 0; i < nNodes; i++)
    {
        count+= nodes[i].writeToFile(handle);
    }
    return count;
}

long GGH15mmaps::readFromFile(FILE* handle, TDMatrixParams* prmBuf)
{
    FHE_TIMER_START;
    assert(params != NULL || prmBuf != NULL); // some pointer must be provided

    long count=0;

    if (prmBuf != NULL)
    {
        TDMatrixParams p;
        count = p.readFromFile(handle);
        assert(p == *prmBuf); // sanity check
        params = prmBuf;      // point to given params
    }
    else
        count = params->readFromFile(handle); // overwrite params from input



    count += An.readFromFile(handle,params);

    // Read the actual nodes from disk, one at a time

    long nNodes;
    count += fread(&nNodes,sizeof(nNodes),1,handle);


    nodes.SetLength(nNodes);
    for (long i = 0; i < (nNodes); i++)
    {
        count+= nodes[i].readFromFile(handle, params);
    }
    return count;
}

// Write one GGH15 node to output
long GGH15node::writeToActualFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = params->writeToFile(handle);

    count += fwrite(&index, sizeof(index),1, handle);
    count += A.writeToFile(handle);

    count += ::writeToFile(P,handle);
    count += T.writeToFile(handle);
    count += Tinv.writeToFile(handle);
    count += Pinv.writeToFile(handle);

    return count;
}

long GGH15node::readFromFile(FILE* handle, TDMatrixParams* prmBuf)
{
    FHE_TIMER_START;
    assert(params != NULL || prmBuf != NULL); // some pointer must be provided

    long count = 0;
    if (prmBuf != NULL)
    {
        TDMatrixParams p;
        count = p.readFromFile(handle);
        assert(p == *prmBuf); // sanity check
        params = prmBuf;     // point to given params
    }
    else
        count = params->readFromFile(handle); // overwrite params from input

    count += fread(&index, sizeof(index),1, handle);
    count += A.readFromFile(handle, params);

    count += ::readFromFile(P,handle);
    count += T.readFromFile(handle, params);

    count += Tinv.readFromFile(handle, params);
    count += Pinv.readFromFile(handle, params);

    return count;
}

long GGH15encoding::writeToFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = fwrite(&from, sizeof(from),1, handle);
    count += fwrite(&to, sizeof(to),1, handle);

    count += data.writeToFile(handle);

    return count;
}

long GGH15encoding::readFromFile(FILE* handle, TDMatrixParams& prms)
{
    FHE_TIMER_START;
    long count = fread(&from, sizeof(from),1, handle);
    count += fread(&to, sizeof(to),1, handle);
    return count + data.readFromFile(handle, &prms);
}

/*************************************************************/
/** Encding routines: there are four such routines,         **/
/** depending on whether or not i=0 and whether or not j=n. **/
/*************************************************************/

/* GGH15 encoding for 0<i,j<n. Set C = Tinv_i * CC * T_j, where
 * CC satisfies A_i*CC = (Pinv_i*S*P_j)*A_j+E (mod q). We assume
 * that toA was already multiplied by P_j.
 */
static void encodeGGH15(CRTmatrix& C, const mat_l& S,
                        const GGH15node& from, const GGH15node& to,
                        const CRTmatrix& toA)
{
    FHE_TIMER_START;
    const TDMatrixParams* prms = from.getPrms();

    // Set B = (Pinv_i x S x P_j) x A + E

    CRTmatrix B = from.getPinv(); // Pinv_i
    B *= S;                       // Pinv_i x S
    B *= toA;                     // (Pinv_i x S x P_j) x A

    mat_l E;  // Choose a small n-by-m noise matrix E
    setSmall(E, S.NumRows(), toA.NumCols(), /*sigma=*/prms-> r);

    B += E;  // (Pinv_i x S x P_j) x A  + E

    // Sample small C s.t. fromA x CC = B
    mat_l CC(INIT_SIZE, prms->m, B.NumCols());

    #ifdef DEBUGPRINT
        cout << "encoding" << endl;
    #endif
    //NTL_GEXEC_RANGE(!bGGH15threading, B.NumCols(), first, last)
    GEXEC_RANGE(!bGGH15threading, B.NumCols(), first, last)
    for (long j=first; j< last; j++)
    {
        // Sample vector x_j s.t. A*x_j = j'th column of B
        Vec<vec_zz_p> colMod;
        B.getColumn(colMod, j); // The j'th column of B in CRT format

        vec_l sampledColumn;    // this is where we put the vector x_j

        from.getA().sampleWithTrapdoor(sampledColumn, colMod);// sample a column

        for (long i=0; i<CC.NumRows(); i++)                    // store it in CC
            CC[i][j] = sampledColumn[i];
    }
    //NTL_GEXEC_RANGE_END
    GEXEC_RANGE_END(!bGGH15threading)

    // Compute C = Tinv_i * CC * T_j
    C = from.getTinv();
    C *= CC;
    C *= to.getT();
}
//------------------------------------------------------------------//

/* GGH15 encoding for i=0<j<n. Set C =(S*P_j*A_j +E) *T_j (mod q).
 * We assume that toA it was already multiplied by P_j.
 */
static void encodeGGH15first(CRTmatrix& C, const mat_l& S,
                             const GGH15node& to, const CRTmatrix& toA)
{
    FHE_TIMER_START;
    const TDMatrixParams* prms = to.getPrms();

    // Set C = S x P_j x A + E
    const mat_l& tmpS = S;

    mat_l E;  // Choose a small 1-by-m noise matrix E
    setSmall(E, tmpS.NumRows(), toA.NumCols(), /*sigma=*/prms->r);

    C = toA;            // P_j x A
    C.leftMultBy(tmpS); // S x P_j x A
    C += E;             // S x P_j x A + E

    C *= to.getT();
}
//------------------------------------------------------------------//

/* GGH15 encoding for 0<i<j=n. Set C = Tinv_i * CC, where CC
 * satisfies A_i*CC = (Pinv_i*S)*A_n + E (mod q).
 */
static void encodeGGH15last(CRTmatrix& C, const mat_l& S,
                            const GGH15node& from, const CRTmatrix& toA)
{
    FHE_TIMER_START;
    const TDMatrixParams* prms = from.getPrms();

    // Set B = (Pinv_i x S) x A_n + E

    CRTmatrix B = from.getPinv(); // Pinv_i
    B *= S;                       // Pinv_i x S
    B *= toA;                     // Pinv_i x S x A_n

    mat_l E;  // Choose a small n-by-m noise matrix E
    setSmall(E, B.NumRows(), B.NumCols(), /*sigma=*/prms->r);

    B += E;                       // Pinv_i x S x A_n + E

    // Sample small CC s.t. fromA x CC = B
    mat_l CC(INIT_SIZE, prms->m, B.NumCols());
    for (long j=0; j<B.NumCols(); j++)
    {
        // Sample vector x_j s.t. A*x_j = j'th column of B

        Vec<vec_zz_p> colMod;
        B.getColumn(colMod, j); // The j'th column of B in CRT format

        vec_l sampledColumn;    // this is where we put the vector x_j

        from.getA().sampleWithTrapdoor(sampledColumn, colMod);// sample a column
        for (long i=0; i<CC.NumRows(); i++)                    // store it in CC
            CC[i][j] = sampledColumn[i];
    }

    // Compute C = Tinv_i * CC
    C = from.getTinv();
    C *= CC;
}
//------------------------------------------------------------------//

/* GGH15 encoding for i=0<j=n, set C =(S*A_n +E) (mod q).
 */
static void encodeGGH15_0_n(CRTmatrix& C, const mat_l& S, const CRTmatrix& toA)
{
    FHE_TIMER_START;
    const mat_l& tmpS = S;

    C = toA;            // A_n
    C.leftMultBy(tmpS); // u x S x A

    mat_l E;  // Choose a small 1-by-m noise matrix E
    setSmall(E, tmpS.NumRows(), toA.NumCols(), /*sigma=*/toA.params->r);

    C += E;             // u x S x A + E
}

bool encodeMatrix(const std::string& dirName, GGH15encoding& C, const mat_l& S,
		  long i, long j, GGH15node* iNode, GGH15node* jNode,
		  CRTmatrix& An, TaggedCRTmatrix* toAptr,
		  TDMatrixParams *params, int nNodes)

{
#ifdef DEBUGPRINT
    cout << "begin encode matrix" << endl;
#endif // DEBUGPRINT

    FHE_TIMER_START;

    #ifdef DEBUGPRINT
        cout << "j=" << ToString(j) << endl;
    #endif

    // FIXME: should we enforce i<j ??
    assert(i>=0 && i<nNodes-1);
    assert(j>0  && j<nNodes);

    // Get the destination matrix A in explicit CRT form
    TaggedCRTmatrix toA(params);
    if (toAptr==NULL) toAptr = &toA;

    bool isLast = (nNodes-1 == j);

#ifdef DEBUGPRINT
    cout << "after reading jNode" << endl;
    //getcwd(cwd, sizeof(cwd));
    //cout << "read node info, current directory = " << cwd << endl;
#endif

    // If the caller did not supply a good toA, then compute it now
    if (toAptr->tag != j && !isLast)   // If target is NOT the last node
    {
        jNode->getCRTA(toAptr->data);
        toAptr->data.leftMultBy( jNode->getP() ); // set toA = P_j x A_j
        toAptr->tag = j;
    }

    // Do the actual encoding
    C.setFrom(i);
    C.setTo(j);
    printDot;
    if (i == 0)    // encoding wrt 0 -> ...
    {
      if (isLast)  // encoding wrt 0 -> n
      {
#ifdef DEBUGPRINT
	cout << "encoding encodeGGH15_0_n" << endl;
#endif
	encodeGGH15_0_n(C.getData(), S, An);
      }
      else         // encoding wrt 0 -> j<n
      {
#ifdef DEBUGPRINT
	cout << "encoding first" << endl;
#endif
	encodeGGH15first(C.getData(), S, *jNode, toAptr->data);
      }
    }
    else           // encoding wrt 0<i -> ...
    {

      if (isLast)  // encoding wrt 0<i -> n
      {
#ifdef DEBUGPRINT
	cout << "Encoding last" << endl;
#endif
	encodeGGH15last(C.getData(), S, *iNode, An);
      }
      else         // encoding wrt 0<i -> j<n
      {
	encodeGGH15(C.getData(), S, *iNode, *jNode, toAptr->data);
      }
    }
    return true;
}

bool readEncodeMatrix(const std::string& dirName, GGH15encoding& C, const mat_l& S, long i, long j,
                  TaggedCRTmatrix* toAptr, TDMatrixParams *params, int nNodes)
{
#ifdef DEBUGPRINT
    cout << "begin encode matrix" << endl;
#endif // DEBUGPRINT

    FHE_TIMER_START;
    // FIXME: should we enforce i<j ??
    assert(i>=0 && i<nNodes-1);
    assert(j>0  && j<nNodes);

    // Get the destination matrix A in explicit CRT form
    TaggedCRTmatrix toA(params);
    if (toAptr==NULL) toAptr = &toA;

    bool isLast = (nNodes-1 == j);

    GGH15node jNode, iNode;

#ifdef DEBUGPRINT
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    cout << "read node info, current directory = " << cwd << endl;
    cout << "reading jNode" << endl;
#endif

//  load the correct node
    std::string fileName = "node"+ToString(j)+ ".dat";
    FILE* handle = fopen(fileName.c_str(), "rb");
    if (handle != 0)
    {
        jNode.readFromFile(handle, params);
        fclose(handle);
    }

#ifdef DEBUGPRINT
    cout << "after reading jNode" << endl;
    getcwd(cwd, sizeof(cwd));
    cout << "read node info, current directory = " << cwd << endl;
#endif
    //load An
    CRTmatrix An;// The An matrix (actually n-by-1 "matrix")
    fileName = "nAn.dat";
    handle = fopen(fileName.c_str(), "rb");
    if (handle != 0)
    {
        long n;
        long count = fread(&n, sizeof(n),1, handle);
        count+=An.readFromFile(handle, params);
        fclose(handle);
    }
    else NTL::Error("Cannot read An");

    // If the caller did not supply a good toA, then compute it now
    if (toAptr->tag != j && !isLast)   // If target is NOT the last node
    {
        jNode.getCRTA(toAptr->data);
        toAptr->data.leftMultBy( jNode.getP() ); // set toA = P_j x A_j
        toAptr->tag = j;
    }

    // Do the actual encoding
    C.setFrom(i);
    C.setTo(j);
    printDot;
    if (i == 0)    // encoding wrt 0 -> ...
    {
      if (isLast)  // encoding wrt 0 -> n
      {
#ifdef DEBUGPRINT
	cout << "encoding encodeGGH15_0_n" << endl;
#endif
	encodeGGH15_0_n(C.getData(), S, An);
      }
      else         // encoding wrt 0 -> j<n
      {
#ifdef DEBUGPRINT
	cout << "encoding first" << endl;
#endif
	encodeGGH15first(C.getData(), S, jNode, toAptr->data);
      }
    }
    else           // encoding wrt 0<i -> ...
    {
      //get iNode - we are already in the correct directory
      fileName = "node"+ToString(i)+ ".dat";
      handle = fopen(fileName.c_str(), "rb");
      if (handle != 0)
      {
	iNode.readFromFile(handle , params);
	fclose(handle);
      }
      if (isLast)  // encoding wrt 0<i -> n
      {
#ifdef DEBUGPRINT
	cout << "Encoding last" << endl;
#endif
	encodeGGH15last(C.getData(), S, iNode, An);
      }
      else         // encoding wrt 0<i -> j<n
      {
	encodeGGH15(C.getData(), S, iNode, jNode, toAptr->data);
      }
    }
    return true;
}

// Encode plaintext matrix wrt the path i -> j. The toAptr argument is
// an optimization:
//  + If toAptr is non-NULL and toAptr->tag==j, then it is assumed
//    that to toAptr->data is the matrx P_j x A_j.
//  + If toAptr is non-NULL and toAptr->tag!=j, then toA->data is set
//    P_j x A_j and toAptr->tag is set to j.
bool GGH15mmaps::encode(GGH15encoding& C, const mat_l& S, long i, long j,
                        TaggedCRTmatrix* toAptr)
{
    FHE_TIMER_START;
    // FIXME: should we enforce i<j ??
    assert(i>=0 && i<numNodes()-1);
    assert(j>0  && j<numNodes());

    // Get the destination matrix A in explicit CRT form
    TaggedCRTmatrix toA(params);
    if (toAptr==NULL) toAptr = &toA;

    bool isLast = (numNodes()-1 == j);

    // If the caller did not supply a good toA, then compute it now
    if (toAptr->tag != j && !isLast)   // If target is NOT the last node
    {
        nodes[j-1].getCRTA(toAptr->data);
        toAptr->data.leftMultBy( nodes[j-1].getP() ); // set toA = P_j x A_j
        toAptr->tag = j;
    }

    // Do the actual encoding
    C.from = i;
    C.to = j;

    printDot;
    if (i == 0)    // encoding wrt 0 -> ...
    {
        if (isLast)  // encoding wrt 0 -> n
            encodeGGH15_0_n(C.data, S, An);

        else         // encoding wrt 0 -> j<n
            encodeGGH15first(C.data, S, nodes[j-1], toAptr->data);
    }
    else           // encoding wrt 0<i -> ...
    {

        if (isLast)  // encoding wrt 0<i -> n
            encodeGGH15last(C.data, S, nodes[i-1], An);

        else         // encoding wrt 0<i -> j<n
            encodeGGH15(C.data, S, nodes[i-1], nodes[j-1], toAptr->data);

    }
    return true;
}


// A debugging method to check that an encoding was done right
bool GGH15mmaps::verifyEncoding(const GGH15encoding& C,	const mat_l& S) const
{
    FHE_TIMER_START;
    long i = C.fromNode();
    long j = C.toNode();

    CRTmatrix CC(params), SS(params), toA(params), fromA(params);
    if (i!=0)
    {
        nodes[i-1].getCRTA(fromA);

        if (j!=(numNodes()-1))   // 0<i -> j<n
        {
            nodes[j-1].getCRTA(toA);

            CC = nodes[i-1].getT();
            CC *= C.data;
            CC *= nodes[j-1].getTinv();

            SS = nodes[i-1].getPinv();
            SS *= S;
            SS *= nodes[j-1].getP();
        }
        else                 // 0<i -> j=n
        {
            toA = An;

            CC = nodes[i-1].getT();
            CC *= C.data;

            SS = nodes[i-1].getPinv();
            SS *= S;
        }
    }
    else // i==0
    {
        if (j!=(numNodes()-1))   // 0=i -> j<n
        {
            nodes[j-1].getCRTA(toA);

            CC = C.data;
            CC *= nodes[j-1].getTinv();

            SS = nodes[j-1].getP();
            SS.leftMultBy(S);
        }
        else                 // 0=i -> j=n
        {
            toA = An;
            CC = C.data;
            SS = S;
        }
    }

    if (i>0)   // Check that CC and EE = fromA*CC - SS*toA are small
    {
        fromA *= CC;
        SS *= toA;
        fromA -= SS;
        return (CC.isSmall() && fromA.isSmall());
    }
    else   // Check that EE = CC - SS*toA is small
    {
        SS *= toA;
        CC -= SS;
        return CC.isSmall();
    }
}

bool operator==(const GGH15mmaps& A, const GGH15mmaps& B)
{
    FHE_TIMER_START;

    if (A.getParams()!=B.getParams())
        return false;

    if (A.numNodes()!=B.numNodes())
        return false;

    for (long iNode = 0; iNode < (A.numNodes()-2); iNode++)
    {
        if (A.getNode(iNode)!=B.getNode(iNode))
            return false;
    }

    if (A.getAn()!=B.getAn())
        return false;

    return true;
}

// checks whether two nodes are identical, i.e., have the exact same
// variable values in them. if so, returns true.

bool operator==(const GGH15node& A, const GGH15node& B)
{
    if (A.getIndex()!=B.getIndex())
        return false;
    if (A.getPrms()!=B.getPrms())
        return false;
    if (A.getA()!=B.getA())
        return false;
    if (A.getP()!=B.getP())
        return false;
    if (A.getPinv()!=B.getPinv())
        return false;
    if (A.getT()!=B.getT())
        return false;
    if (A.getTinv()!=B.getTinv())
        return false;

    return true;
}
