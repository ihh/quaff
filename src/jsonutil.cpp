#include "jsonutil.h"
#include "util.h"
  
JsonValue* JsonUtil::find (JsonValue& parent, const char* key) {
  Assert (parent.getTag() == JSON_OBJECT, "JSON value is not an object");
  for (auto i : parent)
    if (strcmp (i->key, key) == 0)
      return &i->value;
  return NULL;
}

JsonValue& JsonUtil::findOrDie (JsonValue& parent, const char* key) {
  JsonValue* val = find (parent, key);
  Assert (val != NULL, "Couldn't find JSON tag %s", key);
  return *val;
}

map<string,JsonValue*> JsonUtil::objectAsMap (JsonValue& obj) {
  map<string,JsonValue*> m;
  Assert (obj.getTag() == JSON_OBJECT, "JSON value is not an object");
  for (auto n : obj)
    m[string(n->key)] = &n->value;
  return m;
}

vector<double> JsonUtil::numericArrayAsVector (JsonValue& arr) {
  vector<double> v;
  Assert (arr.getTag() == JSON_ARRAY, "JSON value is not an array");
  for (auto n : arr) {
    Assert (n->value.getTag() == JSON_NUMBER, "JSON value is not a number");
    v.push_back (n->value.toNumber());
  }
  return v;
}

string JsonUtil::quoteEscaped (const string& str) {
  string esc = "\"";
  writeEscaped (str, back_inserter(esc));
  esc += "\"";
  return esc;
}
