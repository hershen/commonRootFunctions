#include <string>

#include "LinearInterp.h"

namespace myFuncs {

std::string getAlpMassString(const std::string& text);

double getAlpMassDouble(const std::string& text);

//x,y,z in cm
//based on https://www.slac.stanford.edu/BFROOT/www/Detector/Calorimeter/software/geometry.html
// bool inEMC(const double x, const double y, const double z); line equation for face of forward endcap is probably wrong, see getDistanceToEmcFace_m


//Based on https://bbr-wiki.slac.stanford.edu/bbr_wiki/index.php/File:EMC_long_cross.png
double getDistanceToEmcFace_m(const double theta_rad);

LinearInterp absDeltaThetaLab12_degMin_Run16(std::vector<double>{0.0, 0.7987624889932496}, std::vector<double>{1.5211127049830904, 0});
LinearInterp absDeltaThetaLab12_degMin_Run7(std::vector<double>{0.0, 0.7987624889932496}, std::vector<double>{1.5211127049830904, 0});

LinearInterp chi2Max_Run16(std::vector<double>{0.06051532033426172, 0.3767270194986073, 9.0, 11.0}, std::vector<double>{0.0, 100.0, 100.0, 4.8357348703171965});
LinearInterp chi2Max_Run7(std::vector<double>{0.05, 0.37, 11.0}, std::vector<double>{20.000000000000004, 100.0, 100.0});

LinearInterp minE12cmMin_Run16(std::vector<double>{0.42315909057346635, 5.95, 6.5, 11.0}, std::vector<double>{0.7, 1.5445748422438028, 3.3499999999999694, 3.577000000000015});
LinearInterp minE12cmMin_Run7(std::vector<double>{0.19, 0.25, 5.5, 11.0}, std::vector<double>{0.7, 1.0, 1.5835052732560575, 4.575059523809527});
LinearInterp minTheta_degMin_Run16(std::vector<double>{-100, 100}, std::vector<double>{22.5, 22.5});
LinearInterp minTheta_degMin_Run7(std::vector<double>{-100, 100}, std::vector<double>{22.5, 22.5});

LinearInterp minAbsAcolPhiCM_degMin_Run16(std::vector<double>{10., 10.3}, std::vector<double>{0, 0.5});
LinearInterp minAbsAcolPhiCM_degMin_Run7(std::vector<double>{10., 10.3}, std::vector<double>{0, 0.5});

const std::vector<LinearInterp> minE12cmMin_interp{minE12cmMin_Run16, minE12cmMin_Run7};
const std::vector<LinearInterp> minTheta_degMin_interp{minTheta_degMin_Run16, minTheta_degMin_Run7};
const std::vector<LinearInterp> chi2Max_interp{chi2Max_Run16, chi2Max_Run7};
const std::vector<LinearInterp> minAbsAcolPhiCM_degMin_interp{minAbsAcolPhiCM_degMin_Run16, minAbsAcolPhiCM_degMin_Run7};
const std::vector<LinearInterp> absDeltaThetaLab12_degMin_interp{absDeltaThetaLab12_degMin_Run16, absDeltaThetaLab12_degMin_Run7};



} // namespace myFuncs
