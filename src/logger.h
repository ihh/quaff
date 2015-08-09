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
  int parseLogArgs (int* argcPtr, char*** argvPtr);
};

#define LogAt(V)     (logger != NULL && logger->verbosity >= (V))
#define LogWhen(TAG) (logger != NULL && logger->logTags.find(TAG) != logger->logTags.end())
#define LogThis      LogWhen(__FUNCTION__)
#define LogThisAt(V) (LogAt(V) || LogThis)

#endif /* LOGGER_INCLUDED */

