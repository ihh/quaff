#ifndef BASE64_INCLUDED
#define BASE64_INCLUDED

#include <string>

std::string base64_encode(const char* , unsigned int len);
std::string base64_decode(std::string const& s);

#endif /* BASE64_INCLUDED */

