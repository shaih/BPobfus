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
#include <NTL/tools.h>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <utility>

struct DefaultMoverPolicy {
   template<class T>
   static void mover(T& dst, const T& src) { dst = std::move(src); }

};
/*
template<class T, class P=DefaultMoverPolicy>
class SynchronizedPipe {
private:

   SynchronizedPipe(const SynchronizedPipe&); // disabled
   void operator=(const SynchronizedPipe&); // disabled

public:

   SynchronizedPipe();

   bool receive(T& val);
   // called by receiver: indicates receiver is ready for data, and receiver
   // is blocked until data is sent by sender.
   // Returns true if sender has sent data, and false if sender is done sending data.

   void send(T& val);
   // called by sender: sender is blocked until receiver is ready

   void end();
   // called by sender: sender is blocked until receiver is ready
   // this indicates that the sender is done 

};


// NOTE: both receive and send *move* the data using the mover method
// of P, which by default invokes the move assignment operator for T.


*/

template<class T, class P=DefaultMoverPolicy>
class SynchronizedPipe {
private:

   class SimpleSignal {
   private:
     bool flag;
     std::mutex m;
     std::condition_variable cv;
   
     SimpleSignal(const SimpleSignal&); // disabled
     void operator=(const SimpleSignal&); // disabled
   
   public:
     SimpleSignal() : flag(false) { }
   
     void wait() 
     {
       std::unique_lock<std::mutex> lock(m);
       cv.wait(lock, [&]() { return flag; } );
       flag = false;
     }
   
     void signal()
     {
       std::lock_guard<std::mutex> lock(m);
       if (flag) NTL::TerminalError("signaling while active signal");
       flag = true;
       cv.notify_one();
     }
   };
   
   
   class CompositeSignal {
   private:
     int flag; 
     T val;
     std::mutex m;
     std::condition_variable cv;
   
     CompositeSignal(const CompositeSignal&); // disabled
     void operator=(const CompositeSignal&); // disabled
   
   public:
     CompositeSignal() : flag(0) { }
   
     bool wait(T& _val) 
     {
       bool res = false;
       std::unique_lock<std::mutex> lock(m);
       cv.wait(lock, [&]() { return flag; } );
       if (flag > 0) {
          P::mover(_val, val);
          res = true;
       }
       flag = 0;
       return res;
     }
   
     void signal(T& _val)
     {
       std::lock_guard<std::mutex> lock(m);
       if (flag) NTL::TerminalError("signaling while active signal");
       P::mover(val, _val);
       flag = 1;
       cv.notify_one();
     }
   
     void end()
     {
       std::lock_guard<std::mutex> lock(m);
       if (flag) NTL::TerminalError("signaling while active signal");
       flag = -1;
       cv.notify_one();
     }
   };



   SimpleSignal receiverSignal;
   CompositeSignal senderSignal;

   SynchronizedPipe(const SynchronizedPipe&); // disabled
   void operator=(const SynchronizedPipe&); // disabled

public:

   SynchronizedPipe() { }

   bool receive(T& val)
   {
      receiverSignal.signal();
      return senderSignal.wait(val);
   }

   void send(T& val)
   {
      receiverSignal.wait();
      senderSignal.signal(val);
   }

   void end()
   {
      receiverSignal.wait();
      senderSignal.end();
   }
   
   
};
   

// Example program follows

#include <NTL/BasicThreadPool.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>

NTL_CLIENT

void Process1(SynchronizedPipe<int>& pipe1)
{

   for (;;) {
      int x;
      cin >> x;
      if (x == 0) break;
      pipe1.send(x);
   }

   pipe1.end();
   


} 

void Process2(SynchronizedPipe<int>& pipe1, SynchronizedPipe<int>& pipe2)
{

   for (;;) {
      int x;
      if (!pipe1.receive(x)) break;
      x = 2*x;
      this_thread::sleep_for(chrono::seconds(2));
      pipe2.send(x);
   }

   pipe2.end();

} 

void Process3(SynchronizedPipe<int>& pipe2)
{

   for (;;) {
      int x;
      if (!pipe2.receive(x)) break;
      cout << x << "\n";
   }

} 


int main()
{

   SynchronizedPipe<int> pipe1, pipe2;

   SetNumThreads(4);

   NTL_EXEC_INDEX(4, index)

      switch (index) {
      case 0: break;
      // NOTE: Starting with NTL v9.10, the current thread
      // always gets assigned index == 0. This can be convenient to know:
      // the current thread's thread pool (which is thread_local) is already
      // in use; however, the other threads' thread pools are not in use. 
      // Thus, if we want, the functions Process1, Process2, Process3 (below)
      // could each independently call SetNumThreads to work with their own
      // thread pools.

      case 1: Process1(pipe1); break;

      case 2: Process2(pipe1, pipe2); break;

      case 3: Process3(pipe2); break;
      }

   NTL_EXEC_INDEX_END
}

