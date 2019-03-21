#include <string>

namespace myFuncs {

std::string getAlpMassString(const std::string& text);

double getAlpMassDouble(const std::string& text);

//x,y,z in cm
//based on https://www.slac.stanford.edu/BFROOT/www/Detector/Calorimeter/software/geometry.html
// bool inEMC(const double x, const double y, const double z); line equation for face of forward endcap is probably wrong, see getDistanceToEmcFace_m

} // namespace myFuncs
