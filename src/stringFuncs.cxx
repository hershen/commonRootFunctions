#include "stringFuncs.h"
#include <sstream>

std::vector<std::string> &myFuncs::splitString(const std::string &s, char delim, std::vector<std::string> &elements) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elements.push_back(item);
  }
  return elements;
}

std::vector<std::string> myFuncs::splitString(const std::string &input, char delim) {
  std::vector<std::string> elements;
  splitString(input, delim, elements);
  return elements;
}
