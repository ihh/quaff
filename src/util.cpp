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

std::string plural (long n, const char* singular) {
  std::string s = std::to_string(n) + " " + singular;
  if (n > 1)
    s += "s";
  return s;
}

ProgressLogger::ProgressLogger() : msg(NULL)
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

  logger << msg << ": started at " << asctime(timeinfo);
  logger.unlock();
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

    logger.lock();

    fprintf (stderr, "%s: ", msg);
    va_start (argptr, desc);
    vfprintf (stderr, desc, argptr);
    va_end (argptr);
    fprintf (stderr, ". Estimated time left: ");
    if (estimatedDaysLeft > 2)
      fprintf (stderr, "%g days", estimatedDaysLeft);
    else if (estimatedHoursLeft > 2)
      fprintf (stderr, "%g hrs", estimatedHoursLeft);
    else if (estimatedMinutesLeft > 2)
      fprintf (stderr, "%g mins", estimatedMinutesLeft);
    else
      fprintf (stderr, "%g secs", estimatedSecondsLeft);
    fprintf (stderr, " (%g%%)\n", 100*completedFraction);

    logger.unlock();
    
    lastElapsedSeconds = elapsedSeconds;
    reportInterval = fmin (10., 2*reportInterval);
  }
}
