#include <boost/lexical_cast.hpp>

namespace myFuncs {

std::string getAlpMassString(const std::string& text) {
  const std::string tmp = text.substr(4);
  return tmp.substr(0, tmp.find("_"));
}

double getAlpMassDouble(const std::string& text) { return boost::lexical_cast<double>(getAlpMassString(text)); }
} // namespace myFuncs
