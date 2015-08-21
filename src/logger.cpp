#include <regex>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include "logger.h"

Logger logger;

string ansiEscape (int code) {
  return string("\x1b[") + to_string(code) + "m";
}

Logger::Logger()
  : verbosity(0), lastTestedVerbosity(-1), lastColor(-1), mxLocked(false), useAnsiColor(true)
{
  for (int col : { 7, 2, 3, 5, 6, 1, 4 })  // roughly, brighter colors then darker ones
    logAnsiColor.push_back (ansiEscape(30 + col) + ansiEscape(40));
  threadAnsiColor = ansiEscape(37) + ansiEscape(41);  // white on red
  ansiColorOff = ansiEscape(0);
}

void Logger::addTag (const char* tag) {
  addTag (string (tag));
}

void Logger::addTag (const string& tag) {
  logTags.insert (tag);
}

void Logger::setVerbose (int v) {
  verbosity = max (verbosity, v);
}

void Logger::colorOff() {
  useAnsiColor = false;
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

    } else if (arg == "-nocolor") {
      useAnsiColor = false;
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

Logger& Logger::lockAndPrint (bool showHeader) {
  thread::id myId = this_thread::get_id();
  if (!(mxLocked && mxOwner == myId)) {
    if (mx.try_lock_for (std::chrono::milliseconds(1000))) {
      if (showHeader && mxOwner != myId && threadNum.size() > 1)
	clog << (useAnsiColor ? threadAnsiColor.c_str() : "")
	     << '(' << threadName(myId) << ')'
	     << (useAnsiColor ? ansiColorOff.c_str() : "") << ' ';
      mxOwner = myId;
      mxLocked = true;
    } else if (showHeader)
      clog << (useAnsiColor ? threadAnsiColor.c_str() : "")
	   << '(' << threadName(myId) << ", ignoring lock by " << threadName(mxOwner) << ')'
	   << (useAnsiColor ? ansiColorOff.c_str() : "") << ' ';
  }
  if (showHeader && useAnsiColor && lastTestedVerbosity != lastColor) {
    lastColor = lastTestedVerbosity;
    clog << (lastTestedVerbosity < 0
	     ? logAnsiColor.front()
	     : (lastTestedVerbosity >= (int) logAnsiColor.size()
		? logAnsiColor.back()
		: logAnsiColor[lastTestedVerbosity]));
  }
  return *this;
}

Logger& Logger::unlockAndPrint() {
  thread::id myId = this_thread::get_id();
  if (mxLocked && mxOwner == myId) {
    mxLocked = false;
    mx.unlock();
  }
  if (useAnsiColor) {
    clog << ansiColorOff;
    lastColor = -1;
  }
  return *this;
}

Logger& Logger::lock() { return lockAndPrint(true); }
Logger& Logger::unlock() { return unlockAndPrint(); }

Logger& Logger::lockSilently() { return lockAndPrint(false); }
Logger& Logger::unlockSilently() { return unlockAndPrint(); }

string Logger::threadName (thread::id id) {
  string s;
  const auto& iter = threadNum.find(id);
  if (iter == threadNum.end()) {
    ostringstream o;
    o << "thread " << id;
    s = o.str();
  } else
    s = string("thread #") + to_string (iter->second);
  return s;
}

void Logger::assignThreadName (const thread& thr) {
  threadNum[thr.get_id()] = (unsigned int) (threadNum.size() + 1);
}

void Logger::clearThreadNames() {
  threadNum.clear();
}

ProgressLogger::ProgressLogger (int verbosity, const char* function, const char* file)
  : msg(NULL), verbosity(verbosity), function(function), file(file)
{ }

void ProgressLogger::initProgress (const char* desc, ...) {
  startTime = std::chrono::system_clock::now();
  lastElapsedSeconds = 0;
  reportInterval = 2;

  time_t rawtime;
  struct tm * timeinfo;

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  
  va_list argptr;
  va_start (argptr, desc);
  vasprintf (&msg, desc, argptr);
  va_end (argptr);

  if (logger.testVerbosityOrLogTagsWithLock (verbosity, function, file))
    logger << msg << ": started at " << asctime(timeinfo)
	   << flush;  // asctime has a newline, so no endl, but we need a manipulator to unlock
}

ProgressLogger::~ProgressLogger() {
  if (msg)
    free (msg);
}

void ProgressLogger::logProgress (double completedFraction, const char* desc, ...) {
  va_list argptr;
  const std::chrono::system_clock::time_point currentTime = std::chrono::system_clock::now();
  const auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds> (currentTime - startTime).count();
  const double estimatedTotalSeconds = elapsedSeconds / completedFraction;
  if (elapsedSeconds > lastElapsedSeconds + reportInterval) {
    const double estimatedSecondsLeft = estimatedTotalSeconds - elapsedSeconds;
    const double estimatedMinutesLeft = estimatedSecondsLeft / 60;
    const double estimatedHoursLeft = estimatedMinutesLeft / 60;
    const double estimatedDaysLeft = estimatedHoursLeft / 24;

    if (logger.testVerbosityOrLogTagsWithLock (verbosity, function, file)) {
      char *progMsg;
      va_start (argptr, desc);
      vasprintf (&progMsg, desc, argptr);
      va_end (argptr);

      logger << msg << ": " << progMsg << ". Estimated time left: ";
      if (estimatedDaysLeft > 2)
	logger << estimatedDaysLeft << " days";
      else if (estimatedHoursLeft > 2)
	logger << estimatedHoursLeft << " hrs";
      else if (estimatedMinutesLeft > 2)
	logger << estimatedMinutesLeft << " mins";
      else
	logger << estimatedSecondsLeft << " secs";
      logger << " (" << (100*completedFraction) << "%)" << endl;

      free(progMsg);
    }
    
    lastElapsedSeconds = elapsedSeconds;
    reportInterval = fmin (10., 2*reportInterval);
  }
}
