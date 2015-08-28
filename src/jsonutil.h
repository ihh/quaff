#ifndef JSONUTIL_INCLUDED
#define JSONUTIL_INCLUDED

#include <map>
#include <vector>
#include <string>
#include "gason.h"

using namespace std;

struct JsonUtil {
  static JsonValue* find (JsonValue& parent, const char* key);
  static JsonValue& findOrDie (JsonValue& parent, const char* key);
  static map<string,JsonValue*> objectAsMap (JsonValue& obj);
  static vector<double> numericArrayAsVector (JsonValue& arr);
  static string quoteEscaped (const string& str);
};

#endif /* JSONUTIL_INCLUDED */
