#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
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

string pipeToString (const char* command) {
  string result;
  FILE* pipe = popen (command, "r");
  char line[PIPE_BUF_SIZE];

  while (fgets(line, PIPE_BUF_SIZE, pipe))
    result += line;

  return result;
}
