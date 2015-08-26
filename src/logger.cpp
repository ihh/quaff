#include <sstream>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "logger.h"
#include "regexmacros.h"

Logger logger;

// POSIX basic regular expressions
const regex all_v ("-" RE_GROUP(RE_PLUS("v")), regex::basic);
const regex numeric_v ("-v" RE_NUMERIC_GROUP, regex::basic);

// functions
string ansiEscape (int code) {
  return string("\x1b[") + to_string(code) + "m";
}

Logger::Logger()
  : verbosity(0), lastTestedVerbosity(-1), lastColor(-1), mxLocked(false), useAnsiColor(true)
{
  for (int col : { 7, 2, 3, 5, 6, 1, 2, 3, 5, 6 })  // no blue, it's invisible
    logAnsiColor.push_back (ansiEscape(30 + col) + ansiEscape(40));
  threadAnsiColor = ansiEscape(37) + ansiEscape(41);  // white on red
  ansiColorOff = ansiEscape(0);

  setThreadName (this_thread::get_id(), "main thread");
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

string Logger::args() const {
  string a;
  if (verbosity > 0)
    a += " -v" + to_string(verbosity);
  for (const auto& t : logTags)
    a += " -log " + t;
  if (!useAnsiColor)
    a += " -nocolor";
  return a;
}

Logger& Logger::lockAndPrint (bool showHeader) {
  thread::id myId = this_thread::get_id();
  if (!(mxLocked && mxOwner == myId)) {
    if (mx.try_lock_for (std::chrono::milliseconds(1000))) {
      if (showHeader && mxOwner != myId && threadName.size() > 1)
	clog << (useAnsiColor ? threadAnsiColor.c_str() : "")
	     << '(' << getThreadName(myId) << ')'
	     << (useAnsiColor ? ansiColorOff.c_str() : "") << ' ';
      mxOwner = myId;
      mxLocked = true;
    } else if (showHeader)
      clog << (useAnsiColor ? threadAnsiColor.c_str() : "")
	   << '(' << getThreadName(myId) << ", ignoring lock by " << getThreadName(mxOwner) << ')'
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

string Logger::getThreadName (thread::id id) {
  const auto& iter = threadName.find(id);
  if (iter == threadName.end()) {
    ostringstream o;
    o << "thread " << id;
    return o.str();
  }
  return iter->second;
}

void Logger::setThreadName (thread::id id, const string& name) {
  threadName[id] = name;
}

void Logger::nameLastThread (const list<thread>& threads, const char* prefix) {
  setThreadName (threads.back().get_id(),
		 string(prefix) + " thread #" + to_string(threads.size()));
}

void Logger::eraseThreadName (const thread& thr) {
  threadName.erase (thr.get_id());
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

    if (completedFraction > 0 && logger.testVerbosityOrLogTagsWithLock (verbosity, function, file)) {
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
