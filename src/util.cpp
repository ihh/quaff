#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "util.h"

void Warn(const char* warning, ...) {
  va_list argptr;
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
    throw;
  }
}

clock_t progress_startTime;
double progress_lastElapsedSeconds, progress_reportInterval;
char* progress_desc = NULL;
void initProgress (const char* desc, ...) {
  progress_startTime = clock();
  progress_lastElapsedSeconds = 0;
  progress_reportInterval = 2;

  time_t rawtime;
  struct tm * timeinfo;

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  
  if (progress_desc)
    free (progress_desc);

  va_list argptr;
  va_start (argptr, desc);
  vasprintf (&progress_desc, desc, argptr);
  va_end (argptr);
  fprintf (stderr, "%s: started at %s", progress_desc, asctime(timeinfo));
}

void logProgress (double completedFraction, const char* desc, ...) {
  va_list argptr;
  const clock_t currentTime = clock();
  const double elapsedSeconds = ((double) (currentTime - progress_startTime)) / CLOCKS_PER_SEC;
  const double estimatedTotalSeconds = elapsedSeconds / completedFraction;
  if (elapsedSeconds > progress_lastElapsedSeconds + progress_reportInterval) {
    const double estimatedSecondsLeft = estimatedTotalSeconds - elapsedSeconds;
    const double estimatedMinutesLeft = estimatedSecondsLeft / 60;
    const double estimatedHoursLeft = estimatedMinutesLeft / 60;
    const double estimatedDaysLeft = estimatedHoursLeft / 24;
    fprintf (stderr, "%s: ", progress_desc);
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
    progress_lastElapsedSeconds = elapsedSeconds;
    progress_reportInterval = fmin (10., 2*progress_reportInterval);
  }
}
