#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <numeric>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <cassert>

/* uncomment to enable NaN checks */
#define NAN_DEBUG

/* Errors, warnings, assertions.
   Fail(...) and Require(...) are quieter versions of Abort(...) and Assert(...)
   that do not print a stack trace or throw an exception,
   but merely call exit().
   Test(...) does not exit or throw an exception,
   just prints a warning and returns false if the assertion fails.
   Desire(...) is a macro wrapper for Test(...)
   that returns false from the calling function if the test fails.
*/
void Abort(const char* error, ...);
void Assert(int assertion, const char* error, ...);
void Warn(const char* warning, ...);
void Fail(const char* error, ...);
void Require(int assertion, const char* error, ...);
bool Test(int assertion, const char* error, ...);
#define Desire(...) do { if (!Test(__VA_ARGS__)) return false; } while (0)

/* singular or plural? */
std::string plural (long n, const char* singular);

/* stringify */
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

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

/* vector sum */
template <typename T>
std::vector<T> vector_sum(const std::vector<T>& a, const std::vector<T>& b)
{
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), 
		 std::back_inserter(result), std::plus<T>());
  return result;
}

/* vector-scalar product */
template <typename T>
std::vector<T> vector_scale(const T x, const std::vector<T>& a)
{
  std::vector<T> result = a;
  for (auto& y : result)
    y *= x;

  return result;
}
    
#endif /* UTIL_INCLUDED */
