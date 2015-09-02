#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED

#include <list>
#include <set>
#include <map>
#include <string>
#include <deque>
#include <mutex>
#include <thread>
#include <ratio>
#include <chrono>
#include <iostream>
#include <sstream>
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
  int verbosity;
  set<string> logTags;
  bool useAnsiColor;
  vector<string> logAnsiColor;
  string threadAnsiColor, ansiColorOff;
  
  recursive_timed_mutex mx;
  thread::id lastMxOwner;
  map<thread::id,string> threadName;

public:
  Logger();
  // configuration
  void addTag (const char* tag);
  void addTag (const string& tag);
  void setVerbose (int v);
  void colorOff();
  bool parseLogArgs (deque<string>& argvec);
  string args() const;
  
  inline bool testVerbosity (int v) {
    return verbosity >= v;
  }

  inline bool testLogTag (const char* tag) {
    return logTags.find(tag) != logTags.end();
  }

  inline bool testVerbosityOrLogTags (int v, const char* tag1, const char* tag2) {
    return verbosity >= v || testLogTag(tag1) || testLogTag(tag2);
  }

  string getThreadName (thread::id id);
  void setThreadName (thread::id id, const string& name);
  void nameLastThread (const list<thread>& threads, const char* prefix);
  void eraseThreadName (const thread& thr);

  void lock (int color = 0, bool banner = true);
  void unlock (bool endBanner = true);
  void lockSilently() { lock(0,false); }
  void unlockSilently() { unlock(false); }

  template<class T>
  void print (const T& t, int v) {
    lock(v,true);
    clog << t;
    unlock(true);
  }
};

extern Logger logger;

#define VFUNCFILE(V) V,__func__,__FILE__

#define LoggingAt(V)     (logger.testVerbosity(V))
#define LoggingThisAt(V) (logger.testVerbosityOrLogTags(VFUNCFILE(V)))
#define LoggingTag(T)    (logger.testLogTag(T))

#define LogStream(V,S) do { ostringstream tmpLog; tmpLog << S; logger.print(tmpLog.str(),V); } while(0)

#define LogAt(V,S)     do { if (LoggingAt(V)) LogStream(V,S); } while(0)
#define LogThisAt(V,S) do { if (LoggingThisAt(V)) LogStream(V,S); } while(0)
#define LogThisIf(X,S) do { if (X) LogStream(0,S); } while(0)


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

#define ProgressLog(PLOG,V) ProgressLogger PLOG (VFUNCFILE(V))

#endif /* LOGGER_INCLUDED */

