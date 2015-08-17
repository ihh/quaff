#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED

#include <set>
#include <string>
#include <deque>
#include "util.h"

using namespace std;

struct Logger {
  int verbosity;
  set<string> logTags;
  Logger() : verbosity(0) { }
  void addTag (const char* tag);
  void addTag (const string& tag);
  void setVerbose (int v);
  bool parseLogArgs (deque<string>& argvec);
};

extern Logger logger;

#define LogAt(V)     (logger.verbosity >= (V))
#define LogWhen(TAG) (logger.logTags.find(TAG) != logger.logTags.end())
#define LogThis      (LogWhen(__FUNCTION__) || LogWhen(__FILE__))
#define LogThisAt(V) (LogAt(V) || LogThis)

#endif /* LOGGER_INCLUDED */

