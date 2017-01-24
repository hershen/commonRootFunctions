#pragma once

//std
#include <string>
#include <vector>

namespace myFuncs
{
  
  std::vector<std::string> &splitString(const std::string &s, char delim, std::vector<std::string> &elements);

  std::vector<std::string> splitString(const std::string &input, char delim);
 
  inline bool endsWith(std::string const & string, std::string const & ending)
  {
    if (ending.size() > string.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
  }
}