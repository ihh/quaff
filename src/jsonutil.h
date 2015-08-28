#ifndef JSONUTIL_INCLUDED
#define JSONUTIL_INCLUDED

#include <map>
#include <vector>
#include <string>
#include "gason.h"
#include "PracticalSocket.h"
#include "regexmacros.h"

using namespace std;

// Default size of receive buffer for sockets
#define RCVBUFSIZE 1024

// Terminator string for socket messages
#define SocketTerminatorString "# EOF"
extern const regex eofRegex;

// JSON object index
struct JsonMap {
  map<string,JsonValue*> m;
  JsonMap() { }
  JsonMap (const JsonValue& value);
  void initMap (const JsonValue& value);
  bool contains (const char* key) const;
  bool contains (const string& key) const;
  bool containsType (const char* key, JsonTag type) const;
  bool containsType (const string& key, JsonTag type) const;
  JsonValue& operator[] (const char* key) const;
  JsonValue& operator[] (const string& key) const;
};

// wrapper for JsonValue with parent string
class ParsedJson : public JsonMap {
private:
  ParsedJson (const ParsedJson&) = delete;
public:
  string str;
  char *buf, *endPtr;
  JsonValue value;
  JsonAllocator allocator;
  int status;
  ParsedJson (const string& s, bool parseOrDie = true);
  ParsedJson (istream& in, bool parseOrDie = true);
  ParsedJson (TCPSocket* sock, bool parseOrDie = true, const regex& terminatorRegex = eofRegex, int bufSize = RCVBUFSIZE);
  ~ParsedJson();
  void parse (const string& s, bool parseOrDie = true);
  bool parsedOk() const { return status == JSON_OK; }
};

// (mostly) JSON-related utility functions
struct JsonUtil {
  static JsonValue* find (const JsonValue& parent, const char* key);
  static JsonValue& findOrDie (const JsonValue& parent, const char* key);
  static vector<double> doubleVec (const JsonValue& arr);
  static vector<size_t> indexVec (const JsonValue& arr);
  static string quoteEscaped (const string& str);

  // helpers
  static string readStringFromStream (istream& in);
  static string readStringFromSocket (TCPSocket* sock, const regex& terminatorRegex = eofRegex, int bufSize = RCVBUFSIZE);
};

#endif /* JSONUTIL_INCLUDED */
