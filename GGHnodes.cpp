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
 * GGHnodes.cpp - a node-by-nide implementation of the graph-based
 *                multilinear maps from [Gentry-Gorbunov-Halevi, TCC 2015]
 *****************************************************************************/
#include <unistd.h>
#include <stdexcept>
#include <sys/stat.h>

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>

#include "mat_l.h"
#include "GGH15.h"
#include "obfus.h"
#include "pipe.h"

//NTL_CLIENT

// forward decleration
static bool initOneNode(long i, long n, TDMatrixParams& p,
			SynchronizedPipe<nodePair>* pipe);

// Initializing an instance of GGH15 MMAPs, without ever keeping all of
// it in memory. Each node is generated separately and stored to disk, and
// similarly the An matrix is generated and stored to disk. This function
// also stores to disk the number of nodes, to facilitate reading them back
// from disk.

void writeNodesToFile(SynchronizedPipe<nodePair>& pipe)
{
  GGH15node* node;
  FILE* handle;
  nodePair nPair;

  while (true) {
    if (!pipe.receive(nPair)) break;
    handle = nPair.a;
    node = nPair.b;
    node->writeToActualFile(handle);
    delete(node);
    fclose(handle);
  }
}

bool initAllNodes(std::string dirName, TDMatrixParams& p,
		  long numNodes, SynchronizedPipe<nodePair>* pipe)
{
    FHE_TIMER_START;
    assert (numNodes>0);

    char cwdOrig[1024];
    getcwd(cwdOrig, sizeof(cwdOrig));

    // Create a directory and write the params in it
    saveParamstoNewDir(dirName,p);

    int ret = chdir(dirName.c_str());  // move into the new directory
    if (ret != 0) {
      cout << "initAllNodes: writing files in "<<cwdOrig<<dirName<<endl;
      throw std::logic_error("Cannot change directory");
    }

    bool success = true;
    for (long i=0; i< numNodes-1; i++) {
      success = success && initOneNode(i, numNodes, p, pipe);
    }
    ret = chdir(cwdOrig);
    return success;
}

// Initialize a single node, then write it to disk
static bool initOneNode(long i, long n, TDMatrixParams& p,
			SynchronizedPipe<nodePair>* pipe)
{
#ifdef DEBUGPRINT
    fprintf(stderr,"initonenode %s %d\n", __FILE__, __LINE__);
#endif

  if (i>0) { // init a GGH15 node
    // initialized this node, then write to file

    nodePair nPair; // an NTL pair consisting of (FILE*,node*)

    std::string fileName = "node"+ToString(i)+ ".dat";
    nPair.a = fopen(fileName.c_str(), "wb");
    if (nPair.a == 0) return false; // cannot open file

    nPair.b = new GGH15node(p, i);

    if (pipe!=NULL) // Send node over pipe to be written
      pipe->send(nPair);
    else { // write it yourself
      nPair.b->writeToFile(nPair.a);
      fclose(nPair.a);
      delete nPair.b;
    }
  }
  else { // init the right bracket
    std::string fileName = "nAn.dat";
    FILE* handle = fopen(fileName.c_str(), "wb");
    if (handle == 0) return false;

    CRTmatrix An;
    randomFill(An, /*rows=*/p.n, /*cols=*/1, p);
    fwrite(&n, sizeof(n),1, handle);
    An.writeToFile(handle);
    fclose(handle);
  }
  return true;
}
