#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <numeric>
#include <vector>

/* uncomment to enable NaN checks */
#define NAN_DEBUG

/* errors, warnings, assertions */
void Abort(const char* error, ...);
void Assert(int assertion, const char* error, ...);
void Warn(const char* warning, ...);

/* progress logging */
void initProgress (const char* desc, ...);
void logProgress (double completedFraction, const char* desc, ...);

/* sgn function
   http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
 */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* index sort
   http://stackoverflow.com/questions/10580982/c-sort-keeping-track-of-indices
 */
template <typename T>
std::vector<size_t> orderedIndices (std::vector<T> const& values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

#endif /* UTIL_INCLUDED */
