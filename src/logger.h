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

// This class implements a quick-and-dirty threadsafe logger.
// EITHER explicitly call lock() and unlock(),
// OR use operator<< to stream objects to the log,
// terminated by a stream manipulator (e.g. endl or flush).
// The stream manipulator is IMPORTANT: without it,
// you will get delays, and (if you are lucky) you will notice
// "ignoring lock by..." messages appearing in the log.
class Logger {
private:
  int verbosity, lastTestedVerbosity, lastColor;
  set<string> logTags;
  bool useAnsiColor;
  vector<string> logAnsiColor;
  string threadAnsiColor, threadAnsiColorOff, logAnsiColorOff;
  
  recursive_timed_mutex mx;
  bool mxLocked;
  thread::id mxOwner;
  map<thread::id,unsigned int> threadNum;

  inline bool testLogTag (const char* tag) {
    return logTags.find(tag) != logTags.end();
  }

  inline bool testVerbosityOrLogTags (int v, const char* tag1, const char* tag2) {
    return verbosity >= v || testLogTag(tag1) || testLogTag(tag2);
  }

public:
  Logger();
  void addTag (const char* tag);
  void addTag (const string& tag);
  void setVerbose (int v);
  bool parseLogArgs (deque<string>& argvec);

  inline bool testVerbosityWithLock (int v) {
    if (verbosity >= v) {
      lock();
      lastTestedVerbosity = v;
      return true;
    }
    return false;
  }

  inline bool testVerbosityOrLogTagsWithLock (int v, const char* tag1, const char* tag2) {
    if (testVerbosityOrLogTags (v, tag1, tag2)) {
      lock();
      lastTestedVerbosity = v;
      return true;
    }
    return false;
  }

  inline bool testLogTagWithLock (const char* tag) {
    if (testLogTag(tag)) {
      lock();
      lastTestedVerbosity = -1;
      return true;
    }
    return false;
  }
  
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

#define VFUNCFILE(V) V,__FUNCTION__,__FILE__

#define LogAt(V)     (logger.testVerbosityWithLock(V))
#define LogWhen(TAG) (logger.testLogTagWithLock(TAG))
#define LogThisAt(V) (logger.testVerbosityOrLogTagsWithLock(VFUNCFILE(V)))


/* progress logging */
struct ProgressLogger {
  std::chrono::system_clock::time_point startTime;
  double lastElapsedSeconds, reportInterval;
  char* msg;
  int verbosity;
  const char *function, *file;
  ProgressLogger (int verbosity, const char* function, const char* file);
  ~ProgressLogger();
  void initProgress (const char* desc, ...);
  void logProgress (double completedFraction, const char* desc, ...);
};

#define PROGRESS_LOGGER(PLOG,V) ProgressLogger PLOG (VFUNCFILE(V))

#endif /* LOGGER_INCLUDED */

