#include "jsonutil.h"
#include "util.h"
#include "logger.h"

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

bool JsonMap::containsType (const char* key, JsonTag type) const {
  return containsType (string(key), type);
}

bool JsonMap::containsType (const string& key, JsonTag type) const {
  return contains(key) && (*this)[key].getTag() == type;
}

JsonValue& JsonMap::operator[] (const char* key) const {
  return (*this) [string (key)];
}

JsonValue& JsonMap::operator[] (const string& key) const {
  const auto i = m.find(key);
  Require (i != m.end(), "Couldn't find %s in JSON", key.c_str());
  return *i->second;
}

ParsedJson::ParsedJson (const string& s, bool parseOrDie) {
  parse (s, parseOrDie);
}

ParsedJson::ParsedJson (istream& in, bool parseOrDie) {
  parse (JsonUtil::readStringFromStream (in), parseOrDie);
}

ParsedJson::ParsedJson (TCPSocket* sock, bool parseOrDie, const regex& terminatorRegex, int bufSize) {
  parse (JsonUtil::readStringFromSocket (sock, terminatorRegex, bufSize), parseOrDie);
}

void ParsedJson::parse (const string& s, bool parseOrDie) {
  LogThisAt(9, "Parsing string:\n" << s << endl);
  str = s;
  buf = new char[str.size() + 1];
  strcpy (buf, str.c_str());
  status = jsonParse (buf, &endPtr, &value, allocator);
  if (parsedOk()) {
    if (value.getTag() == JSON_OBJECT)
      initMap(value);
  } else {
    if (parseOrDie)
      Fail ("JSON parsing error: %s at byte %zd", jsonStrError(status), endPtr - buf);
    else
      Warn ("JSON parsing error: %s at byte %zd", jsonStrError(status), endPtr - buf);
  }
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

vector<size_t> JsonUtil::indexVec (const JsonValue& arr) {
  vector<size_t> v;
  Assert (arr.getTag() == JSON_ARRAY, "JSON value is not an array");
  for (auto n : arr) {
    Assert (n->value.getTag() == JSON_NUMBER, "JSON value is not a number");
    v.push_back ((size_t) n->value.toNumber());
  }
  return v;
}

string JsonUtil::quoteEscaped (const string& str) {
  string esc;
  write_quoted_escaped (str, back_inserter(esc));
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
  LogThisAt(9,"Read string from socket on port " << sock->getLocalPort() << ":" << endl << msg << endl);
  return msg;
}

string JsonUtil::readStringFromStream (istream& in) {
  string s;
  while (in && !in.eof()) {
    string line;
    getline(in,line);
    s += line;
  }
  return s;
}
