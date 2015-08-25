#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include "util.h"
#include "stacktrace.h"
#include "logger.h"

// buffer size for popen
#define PIPE_BUF_SIZE 1024

void Warn(const char* warning, ...) {
  va_list argptr;
  fprintf(stderr,"Warning: ");
  va_start (argptr, warning);
  vfprintf(stderr,warning,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
}

void Abort(const char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  fprintf(stderr,"Abort: ");
  vfprintf(stderr,error,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
  printStackTrace();
  throw;
}

void Assert(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    fprintf(stderr,"Assertion Failed: ");
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
    va_end (argptr);
    printStackTrace();
    throw;
  }
}

void Fail(const char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  vfprintf(stderr,error,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
  exit (EXIT_FAILURE);
}

void Require(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
    va_end (argptr);
    exit (EXIT_FAILURE);
  }
}

bool Test(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
    va_end (argptr);
  }
  return assertion;
}

std::string plural (long n, const char* singular) {
  std::string s = std::to_string(n) + " " + singular;
  if (n != 1)
    s += "s";
  return s;
}

const string TempFile::dir = "/tmp";
unsigned int TempFile::count = 0;
TempFile::TempFile (const std::string& contents, const char* filenamePrefix) {
  mx.lock();
  do {
    fullPath = dir + '/' + filenamePrefix + std::to_string(getpid()) + '.' + std::to_string(++count);
  } while (access(fullPath.c_str(),F_OK ) != -1);
  mx.unlock();
  ofstream out (fullPath);
  Assert (out.is_open() && !out.fail(), "Couldn't write to temp file %s", fullPath.c_str());
  out << contents;
}

TempFile::~TempFile() {
  if (fullPath.size())
    unlink (fullPath.c_str());
}

string join (const vector<string>& s, const char* sep) {
  string j;
  if (s.size()) {
    j = s.front();
    for (size_t n = 1; n < s.size(); ++n)
      j = j + sep + s[n];
  }
  return j;
}

string pipeToString (const char* command, int* status) {
  string result;
  FILE* pipe = popen (command, "r");
  char line[PIPE_BUF_SIZE];

  while (fgets(line, PIPE_BUF_SIZE, pipe))
    result += line;

  const int s = pclose (pipe);
  if (status)
    *status = s;
  
  return result;
}

