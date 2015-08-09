#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/* uncomment to enable NaN checks */
#define NAN_DEBUG

/* errors, warnings, assertions */
void Abort(const char* error, ...);
void Assert(int assertion, const char* error, ...);
void Warn(const char* warning, ...);

/* progress logging */
void initProgress (const char* desc, ...);
void logProgress (double completedFraction, const char* desc, ...);

#endif /* UTIL_INCLUDED */
