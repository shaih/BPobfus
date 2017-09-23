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
/*****************************************************************************
 obfus.cpp - implementing NFA obfuscation
 *****************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#include <cmath>
#include "DGaussSampler.h"
#include "GGH15.h"
#include "obfus.h"
#include "CRTmatrix.h"
#include <NTL/mat_ZZ.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/BasicThreadPool.h>

NTL_CLIENT

//#define PRINTDOT
#ifdef PRINTDOT
#define printDot cerr << "." << std::flush
#else
#define printDot
#endif

//#define DEBUG
//#define DEBUGPRINT

static void prepare4encoding(const NTL::Mat< NTL::Mat<long> >& trans,
			     Mat<GGH15ptxt>& mainGGH15ptxt,
			     Mat<GGH15ptxt>& dummyGGH15ptxt,
			     CRTmatrix& An, TDMatrixParams* params,
			     long nThreads);
static void encodeBranch(const std::string& dirName,
			 const Mat<GGH15ptxt>& branch,
			 CRTmatrix& An, TDMatrixParams* params,
			 SynchronizedPipe<GGH15node*>* inPipe=NULL,
			 SynchronizedPipe<encPair>* outPipe=NULL);
static void readNodesPipe(SynchronizedPipe<GGH15node*>& pipe,
			  const Mat<GGH15ptxt>& mainBranch,
			  const Mat<GGH15ptxt>& dummyBranch,
			  TDMatrixParams* params);
static void writeEncodingPipe(SynchronizedPipe<encPair>& pipe);
static GGH15node* getNode(long idx, TDMatrixParams* params);
static void evalBranch(GGH15encoding& mainPathEnc,
		       const NTL::Vec<int>& indexes,
		       const NTL::Vec<long>& symbolString,
		       TDMatrixParams* params);

// Obfuscate the given NFA and store the result in the directory 'name'
void initObfuscateNFA(const std::string& dirName,
                      const NTL::Mat< NTL::Mat<long> >& trans, long sec,
		      TDMatrixParams* params, bool bPipe, long threads)
{
    FHE_TIMER_START;

    long L = trans.NumRows(); // BP-length, # of edges
    if (L==0) return; // sanity check, nothing to do

    // compute the parameters if they are not provided by caller
    TDMatrixParams p;
    if (params==NULL)
    {
        long e = 3;
        long dim = trans[0][0].NumRows();

#ifdef DEBUG //small parameters in DEBUG mode
        long k = L;
        p.init(dim+2, k, e);

#else // "real" parameters in non-DEBUG mode
        GGH15basicParams bp(COMP_GGH15PRMS, dim, L, sec, e);
        p.init(bp.n, bp.k, bp.e, bp.m);
#endif // DEBUG

        params = &p; // pointer points to local p

	cout << "m=" << p.m << ", kFactors=" << p.kFactors
	 << ", n=" << p.n << ", e=" << p.e
	 << ", sigmaX=" << p.sigmaX << endl;
    }
    if (bPipe) {
      pipedInit(dirName, *params, 2*L, threads);
    } else {
      SetNumThreads(threads);
      initAllNodes(dirName, *params, 2*L);
    }

    // do actual obfuscation
    obfuscateNFASavedNodes(dirName, trans, params, bPipe, threads);
}


// This is called after initAllNodes, so all the nodes are already
// stored in dirName
void obfuscateNFASavedNodes(const std::string& dirName,
                            const NTL::Mat< NTL::Mat<long> >& trans,
			    TDMatrixParams* params, bool pipes, long nThreads)
{
    FHE_TIMER_START;

    cout << "obfuscating\n" << std::flush;

    char cwdOrig[1024];
    if (!getcwd(cwdOrig, sizeof(cwdOrig)))
      throw std::logic_error("Cannot get current directory");

    long ret = chdir(dirName.c_str());  // move into the new directory
    if (ret != 0) {
      cout << "   original directory was "<<cwdOrig<<endl;
      throw std::logic_error("Cannot change directory");
    }

    // Encode each tranition matrix via two encoded randomized matrices

    Mat<GGH15ptxt> mainGGH15ptxt, dummyGGH15ptxt; // main & dummy branches
    CRTmatrix An;             // The A_n matrix (actually n-by-1 "matrix")

    SetNumThreads(3); // For prepare we need 2, for encode we may need 3

    // First randomize the transition matrices and read A_n
    prepare4encoding(trans,mainGGH15ptxt,dummyGGH15ptxt,An,params,nThreads);

    //Next encode the randomized matrices, one branch at a time

    // FIXME: The piped implemtation below is fragile, relying on the reader
    // to read the nodes in exactly the right order that the worker expect
    // them. A more robust implementation would have the worker ask for a
    // node and then the reader read that node, but it is a lot more work
    // to implement this more robust solution (with pre-fetching).

    if (pipes) {
      SynchronizedPipe<GGH15node*> inPipe;
      SynchronizedPipe<encPair> outPipe;
      // "fork" three threads, reader, worker, and writer
      EXEC_INDEX(3, index)
      switch (index) {

      case 0: // The "reader thread, no multi-threading here
	readNodesPipe(inPipe,mainGGH15ptxt,dummyGGH15ptxt,params);
	inPipe.end();
	break;

      case 1: // The "worker" thread, can be multi-threaded
	SetNumThreads((nThreads>2)? (nThreads-2) : 1);
	encodeBranch(dirName, mainGGH15ptxt, An, params, &inPipe, &outPipe);
	encodeBranch(dirName, dummyGGH15ptxt, An, params,&inPipe, &outPipe);
	outPipe.end();
	{GGH15node *tmp; inPipe.receive(tmp);} // clear end-of-pipe signal
	break;

      case 2: // The "writer thread
	writeEncodingPipe(outPipe);
	break;
      }
      EXEC_INDEX_END
    }
    else { // non-piped implementation
      SetNumThreads(nThreads);
      encodeBranch(dirName, mainGGH15ptxt, An, params);
      encodeBranch(dirName, dummyGGH15ptxt, An, params);
    }
    // Restore working directory
    if (!chdir(cwdOrig)) return; // ignore this value
}


//////
//evalNFA gets a bitstring that defines the path
//////
bool evalNFA(const std::string& dirName, const Vec<long>& symbolString,
             TDMatrixParams* params, long nThreads)
{
    FHE_TIMER_START;
    long nEdges = symbolString.length();
    assert(nEdges>0); // sanity check

    char cwd[1024];
    if (!getcwd(cwd, sizeof(cwd))) // store cwd before switching
      throw std::logic_error("Cannot get current directory");
    if (chdir(dirName.c_str()) != 0)
      throw std::logic_error("Cannot change directory " + dirName);

    // read the parameters from file, if not given
    TDMatrixParams p;
    if (params==NULL)
    {
        params = &p;
        FILE* handle = fopen("params.dat", "rb");
        params->readFromFile(handle);
        fclose(handle);
    }

    // The "main" path is 0,1,3,5,...,2*nSymbols-1,
    // the "dummy" path is 0,2,4,..., 2*nSymbols-2, 2*nSymbols-1.
    // First and last nodes are the same on both paths.

    // Along each path we read the matrices that are determined by the
    // symbols in symbolString. For example if symbolString is (0,2,1,...)
    // then on the main path we read 0_1_0, 1_3_2, 3_5_1,...
    // and on the dummy path we read 0_2_0, 2_4_2, 4_6_1,...

    // Then we multiply the matrices along each path (getting two matrices
    // relative to 0-to-(2*nSymbols-1)), subtract these matrices off of each
    // other, and test if the result is a zero-encoding.

    NTL::Vec<int> realIdx(INIT_SIZE, nEdges+1);
    NTL::Vec<int> dummyIdx(INIT_SIZE, nEdges+1);
    for (long i = nEdges-1; i>=0; --i) {
      realIdx[i+1] = 2*i+1;
      dummyIdx[i+1] = (i==nEdges-1)? (2*i +1) : (2*i +2);
    }
    realIdx[0] = dummyIdx[0] = 0;

    // Multiply all matrices, from right to left. Note that rightmost
    // matrices are m-by-1, middle matrices are m-by-m, and leftmost
    // matrices are 7-by-m (?), so right-to-left product will always
    // have a products of m-by-m X m-by-1.

    GGH15encoding mainPathEnc;
    GGH15encoding dummyPathEnc;

    // No. of threads to allocate each branch
    if (nThreads>2) nThreads /= 2;
    else nThreads = 1;

    SetNumThreads(3);
    EXEC_INDEX(3, index)
    switch (index) {
    case 0: break; // Cannot multi-thread theard 0
    case 1:
      SetNumThreads(nThreads);
      evalBranch(mainPathEnc, realIdx, symbolString, params);
      break;
    case 2:
      SetNumThreads(nThreads);
      evalBranch(dummyPathEnc, dummyIdx, symbolString, params);
      break;
    }
    EXEC_INDEX_END

    // subtract and check for zero encoding
    mainPathEnc -= dummyPathEnc;

    bool success = (mainPathEnc.fromNode()==0)
      && (mainPathEnc.toNode()==2*nEdges-1)
      && mainPathEnc.getData().isSmall();

    if (chdir(cwd) != 0) { // restore working directory
      cout << "evalNFA: attempted return to "<<cwd<<endl;
      throw std::logic_error("Cannot change directory");
    }
    return success;
}

//pipedInit forks two threads, one "writer" thread and the "producer" thread, which can be further multi-threaded
bool pipedInit(const std::string& dirName, TDMatrixParams& p, long n,
	       long nThreads)
{

 FHE_TIMER_START;

  bool bsuccess=false;
  SynchronizedPipe<nodePair> pipe;

  // "fork" two thread
  SetNumThreads(2);
  EXEC_INDEX(2, index)
    switch (index) {

    case 0: // The "writer thread, no multi-threading here
      writeNodesToFile(pipe);
      break;

      // NOTE: Starting with NTL v9.10, the current thread
      // always gets assigned index == 0. This can be convenient to know:
      // the current thread's thread pool (which is thread_local) is already
      // in use; however, the other threads' thread pools are not in use.
      // Thus, if we want, the functions Process1, Process2, Process3 (below)
      // could each independently call SetNumThreads to work with their own
      // thread pools.

    case 1:
      SetNumThreads((nThreads>1)? (nThreads-1): 1);
      bsuccess=initAllNodes(dirName, p, n, &pipe);
      pipe.end();
      // The "producer" thread, can be multi-threaded
      break;
    }
  EXEC_INDEX_END

    return bsuccess;
}

// p1 includes the transition matrix and a few extra random dimensions,
// p2 includes the same random dimensions as p1 but the identity or
// part of it instead of the transition matrix
void randomizeTransitions(GGH15ptxt& p1, GGH15ptxt& p2,
                          const NTL::Mat<long>& transition,
                          long i, long j, long nSteps, long dim)
{
    FHE_TIMER_START;
    // set the to,from indexes and the tags
    p1.setFrom((i==0)? (2*i) : (2*i -1));
    p1.setTo(2*i +1);

    p2.setFrom(2*i);
    p2.setTo((i==nSteps-1)? (2*i +1) : (2*i +2));

    std::string tag = ToString(j);
    p1.setTag(tag);
    p2.setTag(tag);

    // Set the data part of the plaintext
    long extraDims = dim - transition.NumRows();  // how many extra dimensions
    mat_l& p1Data = p1.getData();
    mat_l& p2Data = p2.getData();
    p1Data.SetDims(transition.NumRows()+extraDims,
                   transition.NumCols()+extraDims);
    p2Data.SetDims(transition.NumRows()+extraDims,
                   transition.NumCols()+extraDims);
    clear(p1Data);
    clear(p2Data);

    // Copy the transition matrix to p1
    for (long ii=0; ii<transition.NumRows(); ii++)
        for (long jj=0; jj<transition.NumCols(); jj++)
            p1Data[ii][jj] = (long) transition[ii][jj];

    // Set p2 as the identity or part of it
    long smallDim = min(transition.NumRows(), transition.NumCols());
    if (i < nSteps-1) // all but the last step
        for (long ii=0; ii<smallDim/2; ii++) p2Data[ii][ii] = 1;
    if (i > 0)      // all but the first step
        for (long ii=smallDim/2; ii<smallDim; ii++) p2Data[ii][ii] = 1;

    // Choose a small random matrix and store in the extra dimensions of p1,p2
    mat_l rand;
    setSmall(rand, extraDims, extraDims, /*sigma=*/256);
    for (long ii=0; ii<extraDims; ii++) for (long jj=0; jj<extraDims; jj++)
        {
            long row = ii + transition.NumRows();
            long col = jj + transition.NumCols();
            p1Data[row][col] = p2Data[row][col] = rand[ii][jj];
        }
}

/* Function creates a new sub-directory and moves into it. it then creates
 * a binary file "params.dat", opens it and writes the values of p into it.
 * At the end, it moves back to the initial directory.
 */
long saveParamstoNewDir(const std::string& dirName, TDMatrixParams& p)
{
    FHE_TIMER_START;
#if defined(_WIN32) // Windows mkdir takes only one parameter
    int ret = mkdir(dirName.c_str());
#else               // *nix mkdir takes two parameters
    int ret = mkdir(dirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

    char cwdOrig[1024];
    if (!getcwd(cwdOrig, sizeof(cwdOrig)))
      throw std::logic_error("Cannot get current directory");
    ret = chdir(dirName.c_str());  // move into the new directory
    if (ret != 0) {
      cout << "saveParamstoNewDir: attempted chdir to "<<dirName<<endl;
      throw std::logic_error("Cannot change directory");
    }
    FILE* handle = fopen("params.dat", "wb"); // binary file, open for writing
    long count = p.writeToFile(handle);        // write parameters to file
    fclose(handle);

    ret = chdir(cwdOrig);
    if (ret != 0) {
      cout << "saveParamstoNewDir: attempted return to "<<cwdOrig<<endl;
      throw std::logic_error("Cannot change directory");
    }
    return count;
}


//create a random transitions file - for debug purposes,
long makeRandomBP(std::string dirName, std::string name, long dim, long sig, long L)
{
    Mat< Mat<long> > trans(INIT_SIZE, L, sig);
    for (long i=0; i<L; i++) for (long iSig=0; iSig<sig; iSig++) {
	Mat<long>& M = trans[i][iSig];
	M.SetDims(dim,dim);
	if (NTL::RandomBnd(L)==0) clear(M);
	else
	  for (long row=0; row<dim; row++) for (long col=0; col<dim; col++) {
	      M[row][col] = NTL::RandomBnd(2);
	    }
      }

    // Write the transitions matrices to file

    std::string filename = dirName + "/" + name + "_dim" + ToString(dim)
                           + "_sig" + ToString(sig) + "_L" + ToString(L)
                           + ".txt";
    std::fstream fs;
    fs.open(filename.c_str(), std::ios::out); // E.g., "BPs/P_dim2_sig2_L3.txt"
    if (!fs.is_open()) {
        char cwd[1024];
        if (getcwd(cwd, sizeof(cwd)))
          throw std::logic_error("Cannot get current directory");
        std::cout << "Cannot open input file "<< filename
		  << ", working directory=" << cwd << endl;
        exit(0);
    }
    fs << trans << endl;
    fs.close();

    return 0;
}

// Read the A_n matrix from disk and randomize transition matrices
static void prepare4encoding(const NTL::Mat< NTL::Mat<long> >& trans,
			     Mat<GGH15ptxt>& mainGGH15ptxt,
			     Mat<GGH15ptxt>& dummyGGH15ptxt,
			     CRTmatrix& An, TDMatrixParams* params,
			     long nThreads)
{
    long L = trans.NumRows(); // # of steps
    long nLoops = L*trans.NumCols();
    FILE* handle;

    mainGGH15ptxt.SetDims(L, trans.NumCols()); //main branch
    dummyGGH15ptxt.SetDims(L,trans.NumCols()); //dummy branch

    // "fork" two threads, one for reading An, one for randomizing matrices
    EXEC_INDEX(2, index)
    switch (index) {

    case 0: //read An
      handle = fopen("nAn.dat", "rb");
      if (handle == 0) NTL::Error("Cannot read An");
      else {
	long nn; // nn not used for anything
        fread(&nn, sizeof(nn),1, handle);
	An.readFromFile(handle, params);
	fclose(handle);
      }
      break;

    case 1: // randomize transition matrices, multi-threaded
      SetNumThreads((nThreads>1)? (nThreads-1) : 1);
      EXEC_RANGE(nLoops, first, last);
      for (long iLoop = first; iLoop < last; iLoop++) {

        // represent each transition matrix by two randomized matrices,
        // one on the main path and the other on the dummy path

        long j = iLoop % trans.NumCols();   //i = step, j = trans.col
        long i = (iLoop - j) / trans.NumCols();

        GGH15ptxt& p1=mainGGH15ptxt[i][j];
        GGH15ptxt& p2=dummyGGH15ptxt[i][j];

	//we chose the from and to in the following manner:
        p1.setFrom((i==0)? (2*i) : (2*i -1));
        p1.setTo(2*i +1);
        p2.setFrom(2*i);
        p2.setTo((i==L-1)? (2*i +1) : (2*i +2));
        randomizeTransitions(p1, p2, trans[i][j], i, j, L, params->n);
      } //create all ptxt matrices
      EXEC_RANGE_END
    }
    EXEC_INDEX_END
}

static void encodeBranch(const std::string& dirName,
			 const Mat<GGH15ptxt>& branch,
			 CRTmatrix& An, TDMatrixParams* params,
			 SynchronizedPipe<GGH15node*>* inPipe,
			 SynchronizedPipe<encPair>* outPipe)
{
    assert(branch.NumRows()>0 && branch.NumCols()>0);

    GGH15node *toNode=NULL, *fromNode=NULL;
    TaggedCRTmatrix mat1;
    GGH15encoding* enc;

    long lastNode = branch[branch.NumRows()-1][0].getTo();
    for (long i=0; i<branch.NumRows(); i++) {
      long from = branch[i][0].getFrom();
      long to = branch[i][0].getTo();

      // Read from disk the next node that's needed for this step
      if (to >= lastNode) toNode =NULL; // no toNode for last-step encoding
      else if (inPipe != NULL)
        assert(inPipe->receive(toNode));// get node from pipe, failure forbidden
      else
        toNode = getNode(to, params);   // read it yourself

      // Encode all matrices wrt from -> to, one per alphabet symbol
      for (long j=0; j<branch.NumCols(); j++) {
	enc = new(GGH15encoding);
	encodeMatrix(dirName, *enc, branch[i][j].getData(), from, to,
                     fromNode, toNode, An, &mat1, params, lastNode+1);
	// When from==0, encodeMatrix does not use the pointer fromNode

	// Write the encoded matrices to disk
	std::string FileName = ToString(from)+"_"+ToString(to)+"_"
                               + branch[i][j].getTag()+".dat";
	FILE* handle = fopen((const char*)FileName.c_str(), "wb");
	if (handle==0) NTL::Error("Cannot write encoding to disk");

	if (outPipe!=NULL) { // Send matrix over pipe to be written
	  encPair ePair(handle,enc);
	  outPipe->send(ePair);
	} else {             // write it yourself
	  enc->writeToFile(handle);
	  fclose(handle);
	  delete enc;
	}
      }
      if (fromNode!=NULL) delete fromNode;
      fromNode = toNode;
    }
}

static void readNodesPipe(SynchronizedPipe<GGH15node*>& pipe,
			  const Mat<GGH15ptxt>& mainBranch,
			  const Mat<GGH15ptxt>& dummyBranch,
			  TDMatrixParams* params)
{
  GGH15node *node;
  // First read the nodes on the main branch in order
  long lastNode = mainBranch[mainBranch.NumRows()-1][0].getTo();
  for (long i=0; i<mainBranch.NumRows(); i++) {
    long to = mainBranch[i][0].getTo();
    if (to >= lastNode) break;
    node = getNode(to, params);
    pipe.send(node);
  }
  // Next read the nodes on the dummy branch in order
  lastNode = dummyBranch[mainBranch.NumRows()-1][0].getTo();
  for (long i=0; i<dummyBranch.NumRows(); i++) {
    long to = dummyBranch[i][0].getTo();
    if (to >= lastNode) break;
    node = getNode(to, params);
    pipe.send(node);
  }
}

static void writeEncodingPipe(SynchronizedPipe<encPair>& pipe)
{
  GGH15encoding* enc;
  FILE* handle;
  encPair ePair;

  while (true) {
    if (!pipe.receive(ePair)) break;
    handle = ePair.a;
    enc = ePair.b;
    enc->writeToFile(handle);
    delete(enc);
    fclose(handle);
  }
}

static GGH15node* getNode(long idx, TDMatrixParams* params)
{
  std::string fileName = "node"+ToString(idx)+ ".dat";
  FILE* handle = fopen(fileName.c_str(), "rb");
  if (handle==0) NTL::Error("Cannot open node file on disk");
  GGH15node* node = new GGH15node();
  node->readFromFile(handle, params);
  fclose(handle);
  return node;
}

static void evalBranch(GGH15encoding& mainPathEnc,
		       const NTL::Vec<int>& indexes,
		       const NTL::Vec<long>& symbolString,
		       TDMatrixParams* params)
{
  long nEdges = indexes.length()-1;
  for (long i = nEdges-1; i>=0; --i) {
    GGH15encoding cc;

    long from1 = indexes[i];
    long to1   = indexes[i+1];
    std::string fName1 = ToString(from1) + "_"
      + ToString(to1) + "_" + ToString(symbolString[i]) + ".dat";

    // read encodings from files and multiply into the path encodings

    FILE* handle = fopen(fName1.c_str(), "rb");
    if (handle == NULL) // error
      throw std::logic_error("Cannot open file " + fName1);

    if (i == nEdges-1) // read directly to mainPathEnc
      mainPathEnc.readFromFile(handle, *params);
    else {       // multiply into mainPathEnc
      cc.readFromFile(handle, *params);
      mainPathEnc.leftMultBy(cc);
    }
    fclose(handle);
  }
}
