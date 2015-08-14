#include <stdlib.h>
#include "logger.h"

Logger logger;

void Logger::addTag (const char* tag) {
  logTags.insert (string (tag));
}

void Logger::setVerbose (int v) {
  verbosity = max (verbosity, v);
}

bool Logger::parseLogArgs (int& argc, char**& argv) {
  if (argc > 0) {
    const char* arg = argv[0];
    if (strcmp (arg, "-log") == 0) {
      Assert (argc > 1, "%s must have an argument", arg);
      addTag (argv[1]);
      argv += 2;
      argc -= 2;
      return true;

    } else if (strcmp (arg, "-verbose") == 0) {
      setVerbose (1);
      argv += 1;
      argc -= 1;
      return true;

    } else if (arg[0] == '-' && arg[1] == 'v') {
      int all_v = 1, v;
      for (v = 1; arg[v+1] != '\0'; ++v) {
	if (arg[v+1] != 'v') {
	  all_v = 0;
	  break;
	}
      }
      if (all_v) {
	setVerbose (v);
	argv += 1;
	argc -= 1;
	return true;

      } else {
	v = atoi (arg + 2);
	if (v >= 1) {
	  setVerbose (v);
	  argv += 1;
	  argc -= 1;
	  return true;

	}
      }
    }
  }
  return false;
}
