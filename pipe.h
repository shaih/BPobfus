#ifndef _PIPE_H_
#define _PIPE_H_
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
pipe: handles the piping using a mutex for locking. 
Piping is optional in the initialization program,but using it improves efficiancy.
****************************************************************************/
#include <NTL/tools.h>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <utility>

struct DefaultMoverPolicy {

   template<class T>
   static void mover(T& dst, T& src) { dst = std::move(src); }

};

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
      bool ret;
      receiverSignal.signal();
      //cout << "in receive  before receiversignal.wait" << endl;
       ret= senderSignal.wait(val);
       //cout << "in receive  after receiversignal.wait" << endl;
       return ret;
   }

   void send(T& val)
   {
      //cout << "in send before receiversignal.wait" << endl;
      receiverSignal.wait();
     // cout << "in send after receiversignal.wait" << endl;
      senderSignal.signal(val);
   }

   void end()
   {
  // cout << "in end before receiversignal.wait" << endl;
      receiverSignal.wait();
         //cout << "in end after receiversignal.wait" << endl;
      senderSignal.end();
   }


};

#endif // ifndef _PIPE_H_
