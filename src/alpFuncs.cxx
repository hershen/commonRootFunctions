#include <boost/lexical_cast.hpp>

#include "TVector3.h"

namespace myFuncs {

std::string getAlpMassString(const std::string& text) {
  const std::string tmp = text.substr(4);
  return tmp.substr(0, tmp.find("_"));
}

double getAlpMassDouble(const std::string& text) { return boost::lexical_cast<double>(getAlpMassString(text)); }

bool inEMC(const double x, const double y, const double z) {
  const TVector3 vtx(x, y, z);
  const double r = vtx.Perp();

  // Barrel
  if (vtx.z() < 180 and vtx.z() > -230 and r > 92) {
    return true;
  }

  // line equation for face of forward endcap (in z-r plane) r=-2.39z + 522.2
  const double rEndcapFace = -2.39 * vtx.z() + 522.2;
  return r > rEndcapFace;
}
} // namespace myFuncs
