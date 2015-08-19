#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED

#include <set>
#include <map>
#include <string>
#include <deque>
#include <mutex>
#include <thread>
#include <ratio>
#include <chrono>
#include <iostream>
#include "util.h"

using namespace std;

struct Logger {
  int verbosity;
  set<string> logTags;

  timed_mutex mx;
  bool mxLocked;
  thread::id mxOwner;
  map<thread::id,unsigned int> threadNum;
  
  Logger() : verbosity(0), mxLocked(false) { }
  void addTag (const char* tag);
  void addTag (const string& tag);
  void setVerbose (int v);
  bool parseLogArgs (deque<string>& argvec);

  Logger& lock();
  Logger& unlock();

  string threadName (thread::id id);
  void assignThreadName (const thread& thr);
  void clearThreadNames();
  
  typedef ostream& (*ostream_manipulator)(ostream&);
  Logger& operator<< (ostream_manipulator om) {
    clog << om;
    return unlock();
  }

  template<class T>
  Logger& operator<< (const T& t) {
    lock();
    clog << t;
    return *this;
  }
};

extern Logger logger;

#define LogAt(V)     (logger.verbosity >= (V))
#define LogWhen(TAG) (logger.logTags.find(TAG) != logger.logTags.end())
#define LogThis      (LogWhen(__FUNCTION__) || LogWhen(__FILE__))
#define LogThisAt(V) (LogAt(V) || LogThis)

#endif /* LOGGER_INCLUDED */

