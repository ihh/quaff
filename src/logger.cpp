#include <stdlib.h>
#include <regex>
#include "logger.h"

Logger logger;

void Logger::addTag (const char* tag) {
  addTag (string (tag));
}

void Logger::addTag (const string& tag) {
  logTags.insert (tag);
}

void Logger::setVerbose (int v) {
  verbosity = max (verbosity, v);
}

bool Logger::parseLogArgs (deque<string>& argvec) {
  regex all_v ("^-v+$");
  regex numeric_v ("^-v(\\d+)$");
  smatch sm;
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-log") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      addTag (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-verbose") {
      setVerbose (1);
      argvec.pop_front();
      return true;

    } else if (regex_match (arg, all_v)) {
      setVerbose ((unsigned int) (arg.size() - 1));
      argvec.pop_front();
      return true;

    } else if (regex_match (arg, sm, numeric_v)) {
      setVerbose (atoi (sm.str(1).c_str()));
      argvec.pop_front();
      return true;
    }
  }
  return false;
}
