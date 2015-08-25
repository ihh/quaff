#ifndef REGEXMACROS_INCLUDED
#define REGEXMACROS_INCLUDED

#ifdef USE_BOOST
#include <boost/regex.hpp>
using namespace boost;
#else
#include <regex>
using namespace std;
#endif

// POSIX basic regular expressions are used for maximum compatibility,
// since g++ does not stably support ECMAScript regexes yet,
// so we resort to boost.

#define RE_CHAR_CLASS(STR) "[" STR "]"
#define RE_PLUS(CLASS) CLASS CLASS "*"
#define RE_GROUP(EXPR) "\\(" EXPR "\\)"

#define RE_NUMERIC_RANGE "0123456789"
#define RE_ALPHA_RANGE "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
#define RE_ALPHANUM_RANGE RE_ALPHA_RANGE RE_NUMERIC_RANGE

#define RE_NUMERIC_CHAR_CLASS RE_CHAR_CLASS(RE_NUMERIC_RANGE)
#define RE_VARNAME_CHAR_CLASS RE_CHAR_CLASS(RE_ALPHANUM_RANGE "_")
#define RE_DNS_CHAR_CLASS RE_CHAR_CLASS("-" RE_ALPHANUM_RANGE "\\.")
#define RE_FLOAT_CHAR_CLASS RE_CHAR_CLASS("-" RE_NUMERIC_RANGE "eE+\\.")

#define RE_NUMERIC_GROUP RE_GROUP(RE_PLUS(RE_NUMERIC_CHAR_CLASS))
#define RE_VARNAME_GROUP RE_GROUP(RE_PLUS(RE_VARNAME_CHAR_CLASS))
#define RE_DNS_GROUP RE_GROUP(RE_PLUS(RE_DNS_CHAR_CLASS))
#define RE_FLOAT_GROUP RE_GROUP(RE_PLUS(RE_FLOAT_CHAR_CLASS))
#define RE_DOT_GROUP RE_GROUP(RE_PLUS("."))

#endif /* REGEXMACROS_INCLUDED */
