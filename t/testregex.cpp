#include <regex>
#include <iostream>
#include "../src/regexmacros.h"

using namespace std;

void match (const char* re, const char* s) {
  cerr << "Trying '" << s << "' =~ /" << re << "/g ... ";
  const regex r (re, regex_constants::basic);
  if (!regex_match (s, r)) {
    cerr << "no match (not ok)" << endl << '/' << re << "/g should match '" << s << "', but it doesn't" << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  }
  cerr << "match (ok)" << endl;
}

void match (const string& re, const char* s) {
  match (re.c_str(), s);
}

void nomatch (const char* re, const char* s) {
  cerr << "Trying '" << s << "' =~ /" << re << "/g ... ";
  const regex r (re, regex_constants::basic);
  if (regex_match (s, r)) {
    cerr << "match (not ok)" << endl << '/' << re << "/g should not match '" << s << "', but it does" << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  }
  cerr << "no match (ok)" << endl;
}

void nomatch (const string& re, const char* s) {
  nomatch (re.c_str(), s);
}

void matchgroup (const char* re, const char* sc, int n, const char* g) {
  const string s (sc);
  cerr << "Trying '" << s << "' =~ /" << re << "/g ... ";
  const regex r (re, regex_constants::basic);
  smatch sm;
  if (!regex_match (s, sm, r)) {
    cerr << "no match (not ok)" << endl << '/' << re << "/g should match '" << s << "', but it doesn't" << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  } else {
    const string mg = sm.str(n);
    if (mg != g) {
      cerr << "$" << n << " = '" << mg << "', expected '" << g << "' (not ok)" << endl;
      cout << "not ok" << endl;
      exit (EXIT_FAILURE);
    } else
      cerr << "$" << n << " = '" << mg << "' (ok)" << endl;
  }
}

void matchgroup (const string& re, const char* s, int n, const char* g) {
  matchgroup (re.c_str(), s, n, g);
}


void search (const char* re, const char* s) {
  cerr << "Trying '" << s << "' =~ /" << re << "/ ... ";
  const regex r (re, regex_constants::basic);
  if (!regex_search (s, r)) {
    cerr << "no match (not ok)" << endl << '/' << re << "/ should match '" << s << "', but it doesn't" << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  }
  cerr << "match (ok)" << endl;
}

void search (const string& re, const char* s) {
  search (re.c_str(), s);
}

void nosearch (const char* re, const char* s) {
  cerr << "Trying '" << s << "' =~ /" << re << "/ ... ";
  const regex r (re, regex_constants::basic);
  if (regex_search (s, r)) {
    cerr << "match (not ok)" << endl << '/' << re << "/ should not match '" << s << "', but it does" << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  }
  cerr << "no match (ok)" << endl;
}

void nosearch (const string& re, const char* s) {
  nosearch (re.c_str(), s);
}

void searchgroup (const char* re, const char* sc, int n, const char* g) {
  const string s (sc);
  cerr << "Trying '" << s << "' =~ /" << re << "/ ... ";
  const regex r (re, regex_constants::basic);
  smatch sm;
  if (!regex_search (s, sm, r)) {
    cerr << "no match (not ok)" << endl << '/' << re << "/ should match '" << s << "', but it doesn't" << endl;
    cout << "not ok" << endl;
    exit (EXIT_FAILURE);
  } else {
    const string mg = sm.str(n);
    if (mg != g) {
      cerr << "$" << n << " = '" << mg << "', expected '" << g << "' (not ok)" << endl;
      cout << "not ok" << endl;
      exit (EXIT_FAILURE);
    } else
      cerr << "$" << n << " = '" << mg << "' (ok)" << endl;
  }
}

void searchgroup (const string& re, const char* s, int n, const char* g) {
  searchgroup (re.c_str(), s, n, g);
}


int main (int argc, char** argv) {

  const string paramValRegex (" *" RE_VARNAME_GROUP " *: *" RE_DOT_GROUP);
  const string lineRegex (RE_DOT_GROUP);
  const string orderRegex (RE_NUMERIC_GROUP);
  const string countRegex (RE_NUMERIC_GROUP " *: *" RE_FLOAT_GROUP);
  const string remoteUserAddrRegex (RE_DOT_GROUP "@" RE_DNS_GROUP ".*");
  const string remoteAddrRegex (RE_DNS_GROUP ".*");
  const string singleRemotePortRegex (".*:" RE_NUMERIC_GROUP);
  const string multiRemotePortRegex (".*:" RE_NUMERIC_GROUP "-" RE_NUMERIC_GROUP);

  match (RE_CHAR_CLASS("abc"), "a");
  match (RE_CHAR_CLASS("abc"), "b");
  nomatch (RE_CHAR_CLASS("abc"), "d");

  matchgroup (RE_GROUP(RE_CHAR_CLASS("abc")), "a", 1, "a");
  matchgroup (RE_GROUP(RE_CHAR_CLASS("abc")), "b", 1, "b");
  nomatch (RE_GROUP(RE_CHAR_CLASS("abc")), "d");

  match ("a*", "");
  match ("a*", "aaaa");
  nomatch ("a*", "aaab");

  match (paramValRegex, " a: b  ");
  nomatch (paramValRegex, " a 3: b  ");

  matchgroup (paramValRegex, "  a : b  ", 1, "a");
  matchgroup (paramValRegex, "  a : b  ", 2, "b  ");

  match (RE_DNS_CHAR_CLASS, "a");
  match (RE_DNS_CHAR_CLASS, "1");
  match (RE_DNS_CHAR_CLASS, ".");
  match (RE_DNS_CHAR_CLASS, "-");
  nomatch (RE_DNS_CHAR_CLASS, "@");
  nomatch (RE_DNS_CHAR_CLASS, ":");

  matchgroup ("a*" RE_GROUP(RE_PLUS(RE_CHAR_CLASS("abc"))), "aaabc", 1, "bc");

  match (remoteAddrRegex, "127.0.0.1:8000");
  match (remoteAddrRegex, "localhost:8000");
  match (remoteAddrRegex, "localhost");

  cout << "ok" << endl;
  return EXIT_SUCCESS;
}
