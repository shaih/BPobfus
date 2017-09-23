#ifndef _ARGMAP_H_
#define _ARGMAP_H_
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
/********************************************************************
 * argmap.h : Code for parsing command line arguments (from HElib)
 ********************************************************************/
#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#include <cassert>

#if (__cplusplus>199711L)
#include <unordered_map>
#include <memory>
#define TR1_NAMESPACE
#else
#include <tr1/unordered_map>
#include <tr1/memory>
#define TR1_NAMESPACE tr1::
#endif

//! @typedef
typedef std:: TR1_NAMESPACE unordered_map<std::string, const char *> argmap_t;



/**
 * Tries to parse each argument as arg=val, and returns a correspinding map.
 * It returns false if errors were detected, and true otherwise.
 **/
bool parseArgs(int argc,  char *argv[], argmap_t& argmap);

/** Easier arg parsing, example use:
 *
 *   ArgMapping amap;
 *
 *   long p = 2;
 *   amap.arg("p", p, "doc for p");
 *   long m = 0;
 *   amap.arg("m", m, "doc for m", "undefined"); // special default info
 *   long k = 0;
 *   amap.arg("k", k, "doc for k", NULL); // no default info
 *
 *   amap.parse(argc, argv); // parses and overrides initail values
 *                           // of p and p, returns false on error
 *
 *   amap.documentation(); // returns string with documentation
 *                         // for each parameter, one per line,
 **/


/* doArgProcessing: converts c-string s to value T,
 * returns upon success.  By default, we parse using
 * the istream input operator, except when T = string
 * and just convert without any parsing.
 */
template<class T>
bool doArgProcessing(T *value, const char *s)
{
  std::string ss(s);
  std::stringstream sss(ss);
  return bool(sss >> *value);
}
bool doArgProcessing(std::string *value, const char *s);

/* ArgProcessor: virtual base class */

class ArgProcessor {
public:
virtual bool process(const char *s) = 0;
};

/* ArgProcessorDerived: templated subclasses */
template<class T>
class ArgProcessorDerived : public ArgProcessor   {
public:
  T *value;

  virtual bool process(const char *s)
  {
    return doArgProcessing(value, s);
  }

  ArgProcessorDerived(T* _value) : value(_value) {}
};

class ArgMapping {
public:
  std:: TR1_NAMESPACE unordered_map< std::string, std:: TR1_NAMESPACE shared_ptr<ArgProcessor> > map;

  std::stringstream doc;

  // no documentation
  template<class T>
  void arg(const char *name, T& value)
  {
    std::TR1_NAMESPACE shared_ptr<ArgProcessor> ap =
      std::TR1_NAMESPACE shared_ptr<ArgProcessor>(new ArgProcessorDerived<T>(&value));

    assert(!map[name]);
    map[name] = ap;
  }

  // documentation + standard default info
  template<class T>
  void arg(const char *name, T& value, const char *doc1)
  {
    arg(name, value);
    doc << "\t" << name << " \t" << doc1 << "  [ default=" << value << " ]" << "\n";
  }

  // documentation + standard non-standard default info:
  // NULL => no default info
  template<class T>
  void arg(const char *name, T& value, const char *doc1, const char *info)
  {
    arg(name, value);
    doc << "\t" << name << " \t" << doc1;
    if (info)
      doc << "  [ default=" << info << " ]"  << "\n";
    else
      doc << "\n";
  }

  void note(const char *s)
  {
    doc << "\t\t   " << s << "\n";
  }

  void usage(const char *prog)
  {
    std::cerr << "Usage: " << prog << " [ name=value ]...\n";
    std::cerr << documentation();
    exit(0);
  }

  void parse(int argc, char **argv)
  {
    for (long i = 1; i < argc; i++) {
      const char *x = argv[i];
      long j = 0;
      while (x[j] != '=' && x[j] != '\0') j++;
      if (x[j] == '\0') usage(argv[0]);
      std::string name(x, j);
      const char *s = x+j+1;

      std::TR1_NAMESPACE shared_ptr<ArgProcessor> ap = map[name];
      if (!ap) return usage(argv[0]);
      if (!ap->process(s)) usage(argv[0]);
    }
  }

  std::string documentation()
  {
    return doc.str();
  }
};
#endif // _ARGMAP_H_
