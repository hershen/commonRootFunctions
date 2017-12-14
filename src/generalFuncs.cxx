#include "generalFuncs.h"
#include <iostream>

namespace myFuncs {

std::vector<double> changeStringsToDoubles(const std::vector<std::string> &strings) {
  std::vector<double> doubles;
  doubles.reserve(strings.size());

  for_each(strings.begin(), strings.end(), [&doubles](const auto &string) { doubles.push_back(std::atof(string.c_str())); });

  return doubles;
}

} // namespace myFuncs
