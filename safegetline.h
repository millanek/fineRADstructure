#ifndef SAFEGETLINE_H
#define SAFEGETLINE_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

std::istream& safeGetline(std::istream& is, std::string& t); // safely gets a line across platforms, by removing any trailing '\r' characters. If eofok=FALSE, then an error is thrown if end-of-file is encountered, on the expectation that it could be due to missing all the end-of-line characers (i.e. only '\r' are present)

#endif

