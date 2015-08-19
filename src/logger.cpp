#include <stdlib.h>
#include <regex>
#include <sstream>
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

Logger& Logger::lock() {
  thread::id myId = this_thread::get_id();
  if (!(mxLocked && mxOwner == myId)) {
    if (mx.try_lock_for (std::chrono::milliseconds(1000))) {
      if (mxOwner != myId && threadNum.size() > 1)
	clog << "(" << threadName(myId) << ") ";
      mxOwner = myId;
      mxLocked = true;
    } else
      clog << "(" << threadName(myId) << ", ignoring lock by " << threadName(mxOwner) << ") ";
  }
  return *this;
}

Logger& Logger::unlock() {
  thread::id myId = this_thread::get_id();
  if (mxLocked && mxOwner == myId) {
    mxLocked = false;
    mx.unlock();
  }
  return *this;
}

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
