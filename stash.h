#ifndef _SAMPLESTASH_H_
#define _SAMPLESTASH_H_
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
#include <atomic>

// SampleStash - A data structure to keep sampled points.

// For each factor f we have a stash, which is a vector of automic integers.

// Depending on the compile parameter STASH_N, for each residue class
// mod f we keep either one or two such atomic ints. If STASH_N is
// defined as zero then we don't use the stash at all.

#ifndef STASH_N
#define STASH_N 1 // default is a simple stash
#endif

/********************************************************************/
#if STASH_N == 0    // no stash at all
class SampleStash { // dummy class, just no-ops
public:
  SampleStash() {}
  void init(long s) {}
  int getPoint(long i) { return i+1; } // always empty
  void setPoint(long i, int val) {}
};

/********************************************************************/
#elif STASH_N == 1  // simple implementation, one int per residue class
class SampleStash {
  std::atomic<int>* data;
  long lData; //length of data

public:
  // Constructor, destructor
  SampleStash() { data=NULL; }
  ~SampleStash() { if (data!=NULL) delete[] data; }

  // Initialization (Note: this is NOT atomic)
  void init(long s)
  {
        if (data!=NULL || s<=0) return; // don't re-initialize
    data = new std::atomic<int>[s];
    for (long i=0; i<s; i++) data[i].store(i+1);
    lData = s;
  }

  int getPoint(long i) { return data[i].exchange(i+1); }
  // i+1 marks an empty entry

   // binary I/O, returns # of bytes read
    //long writeToFile(FILE* handle) const;
    //long readFromFile(FILE* handle);

  void setPoint(long i, int val) { data[i].store(val); }
};



/********************************************************************/
#else // STASH_N >= 2: keep 2 (or more) ints per residue class

class SampleStash {
  std::atomic<int>* data;

public:
  // Constructor, destructor
  SampleStash() { data=NULL; }
  ~SampleStash() { if (data!=NULL) delete[] data; }

  // Initialization (Note: this is NOT atomic)
  void init(long s)
  {
    if (data!=NULL || s<=0) return; // don't re-initialize
    data = new std::atomic<int>[s*STASH_N];
    for (long i=0; i<s; i++) for (long j=0; j<STASH_N; j++) {
	data[STASH_N*i +j].store(i+1);
      }
  }

  int getPoint(long i) {
    int val;
    for (long j=0; j<STASH_N; j++) {  // look for a value in j'th slots
      val = data[STASH_N*i +j].exchange(i+1);
      if (val != i+1) return val; // found something
    }
    return i+1; // nothing was found
  }

  void setPoint(long i, int val) {
    for (long j=0; j<STASH_N-1; j++) {
      val = data[STASH_N*i +j].exchange(val);
      if (val == i+1) return; // found an empty slot
    }
    data[STASH_N*(i+1)-1].store(val);// step over last value, if you got there
  }
};
#endif //STASH_N

#endif // ifndef _SAMPLESTASH_H_
