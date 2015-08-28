#include "jsonutil.h"
#include "util.h"

const regex eofRegex (SocketTerminatorString, regex_constants::basic);

JsonMap::JsonMap (const JsonValue& value) {
  initMap (value);
}

void JsonMap::initMap (const JsonValue& value) {
  Assert (value.getTag() == JSON_OBJECT, "JSON value is not an object");
  for (auto n : value)
    m[string(n->key)] = &n->value;
}

bool JsonMap::contains (const char* key) const {
  return contains (string (key));
}

bool JsonMap::contains (const string& key) const {
  return m.find(key) != m.end();
}

JsonValue& JsonMap::operator[] (const char* key) const {
  return (*this) [string (key)];
}

JsonValue& JsonMap::operator[] (const string& key) const {
  const auto i = m.find(key);
  Require (i != m.end(), "Couldn't find %s in JSON", key.c_str());
  return *i->second;
}

ParsedJson::ParsedJson (const string& s)
  : str(s), buf(new char[s.size() + 1])
{
  strcpy (buf, str.c_str());
  status = jsonParse (buf, &endPtr, &value, allocator);
  Assert (status == JSON_OK, "JSON parsing error: %s at byte %zd\n%s\n", jsonStrError(status), endPtr - buf, endPtr);
  if (value.getTag() == JSON_OBJECT)
    initMap (value);
}

ParsedJson::~ParsedJson() {
  if (buf) delete[] buf;
}

JsonValue* JsonUtil::find (const JsonValue& parent, const char* key) {
  Assert (parent.getTag() == JSON_OBJECT, "JSON value is not an object");
  for (auto i : parent)
    if (strcmp (i->key, key) == 0)
      return &i->value;
  return NULL;
}

JsonValue& JsonUtil::findOrDie (const JsonValue& parent, const char* key) {
  JsonValue* val = find (parent, key);
  Assert (val != NULL, "Couldn't find JSON tag %s", key);
  return *val;
}

vector<double> JsonUtil::doubleVec (const JsonValue& arr) {
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

string JsonUtil::readStringFromSocket (TCPSocket* sock, const regex& terminatorRegex, int bufSize) {
  string msg;
  char* buf = new char[bufSize];
  int recvMsgSize;
  smatch sm;
  while ((recvMsgSize = sock->recv(buf, bufSize)) > 0) {
    msg.append (buf, (size_t) recvMsgSize);
    if (regex_search (msg, sm, terminatorRegex)) {
      msg = sm.prefix().str();
      break;
    }
  }
  delete[] buf;
  return msg;
}

string JsonUtil::readStringFromStream (istream& in) {
  string s;
  while (!in.eof()) {
    string line;
    getline(in,line);
    s += line;
  }
  return s;
}

ParsedJson* JsonUtil::readJson (istream& in) {
  return new ParsedJson (readStringFromStream (in));
}

ParsedJson* JsonUtil::readJson (TCPSocket* sock, const regex& terminatorRegex, int bufSize) {
  return new ParsedJson (readStringFromSocket (sock, terminatorRegex, bufSize));
}
