#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED

#include <set>
#include <string>
#include "util.h"

using namespace std;

struct Logger {
  int verbosity;
  set<string> logTags;
  Logger() : verbosity(0) { }
  bool parseLogArgs (int& argc, char**& argv);
};

extern Logger logger;

#define LogAt(V)     (logger.verbosity >= (V))
#define LogWhen(TAG) (logger.logTags.find(TAG) != logger.logTags.end())
#define LogThis      LogWhen(__FUNCTION__)
#define LogThisAt(V) (LogAt(V) || LogThis)

#endif /* LOGGER_INCLUDED */

