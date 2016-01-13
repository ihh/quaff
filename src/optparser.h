#ifndef OPTPARSER_INCLUDED
#define OPTPARSER_INCLUDED

#include <string>
#include <deque>

struct OptParser {
  std::deque<std::string> argvec;
  std::string prog, briefText, text;
  std::deque<std::string> implicitSwitches;
  bool unlimitImplicitSwitches;
  OptParser (int argc, char** argv, const char* progName, const char* briefOptsDescription);
  std::string getCommand (const char* error = NULL);
  bool parseUnknown();
};

#endif /* OPTPARSER_INCLUDED */
