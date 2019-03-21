#include <limits>

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
  //Might really be r = -2.659z + 5.707 (all in meters)
  const double rEndcapFace = -2.39 * vtx.z() + 522.2;
  return r > rEndcapFace;
}

double getDistanceToEmcFace_m(const double theta_rad) {
  if(theta_rad < 0.27576202181510407315 or theta_rad > 2.47487687932795934008) { //[15.8, 141.8]
    return std::numeric_limits<double>::infinity();
  }
  if(theta_rad < 0.47249693551706384232) { //endcap
    //Front face of endcap in z-r plane: r = -2.659z+5.707 (all in meters).
    const double z = 5.707/(std::tan(theta_rad)+2.659);
    return z/std::cos(theta_rad);
  }
  if(theta_rad < 2.475) {
    return 0.92/std::sin(theta_rad);
  }
  
  throw; 
}

} // namespace myFuncs
