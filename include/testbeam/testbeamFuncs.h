#pragma once

// STL
#include <memory> //for unique_ptr
#include <unordered_map>

// Linux
#include <dirent.h>

// ROOT
#include "TColor.h"
#include "TCanvas.h"

// Boost
#include "boost/format.hpp"

// Mine
#include "testbeam/constants.h"
#include "histFuncs.h"

class TGraphErrors;
class TF1;

namespace myFuncs {
class PaveText;
namespace testbeam {

constexpr std::array<int, 20> representitiveRuns = {
    591, 597, 599, 600, 589,

    627, // Chinese
    611, 625, 653, 654,

    662, // Babar
    664, 666, 676, 687, 688,

    // 																								 713,
    // //Ukranian
    715,
    // 																								 716,
    717, 730, 731};

constexpr int elecColor = kBlue;
constexpr int muonColor = kRed;
constexpr int pionColor = kGreen;

static const std::unordered_map<int, int> c_pdgToColor{{11, elecColor}, {13, muonColor}, {211, pionColor}};

//-----------------------------------------------------------
// Return Nominal beam momentum string with units
//-----------------------------------------------------------
std::string getNominalMomentumString(const int runNum);

//-----------------------------------------------------------
// Return Measured beam momentum string with units
//-----------------------------------------------------------
std::string getMeasuredMomentumString(const int runNum);

//-----------------------------------------------------------
// Return crystal string
//-----------------------------------------------------------
std::string getCrystalString(const int runNum);

//-----------------------------------------------------------
// Return source distance string with units
//-----------------------------------------------------------
std::string getSourceDistanceString(const int runNum);

inline std::string getRunParamsString(const int runNum) {
  return boost::str(boost::format("Run %1%, %2%, %3%, %4%") % runNum % getCrystalString(runNum) %
                    getMeasuredMomentumString(runNum) % getSourceDistanceString(runNum));
}

//--------------------------------------------------------------------------------------------
// getRunNum
//********************************************************************************************
// Get a filename (may include directories) and return the run number string.
// This is without the partial run number.
//--------------------------------------------------------------------------------------------
inline std::string getRunNum(const std::string &filename) {
  std::string trimLeft = filename.substr(filename.find("00000", 0) + 5);
  return trimLeft.substr(0, trimLeft.find("_000"));
}

//-----------------------------------------------------------
// Get all root files in pathToFiles with runNum in their name and ending in extension.
//-----------------------------------------------------------
std::vector<std::string> getFilesRelatedToRun(const std::string pathToFiles, const int runNum, const std::string = ".root");

inline bool isCrystalChannel(const int channel) { return (channel == 1 or channel == 15); }

// Return TGraphErrors with voltage as function of time
std::unique_ptr<TGraphErrors> getWaveformGraph(const std::vector<double> &voltage,
                                               const std::vector<double> &errors = std::vector<double>());

// Get PaveText of run parameters, channel and event num
std::unique_ptr<myFuncs::PaveText> getRunChannelEventPaveText(const int runNum, const int channelNum, const int eventNum);

//===
// Check if beta is inside the range for the particle and not in range of a different particle
bool isElectron(const double beta, const int runNum);
bool isMuon(const double beta, const int runNum);
bool isPion(const double beta, const int runNum);

//====
// Check if beta is inside a 3sigma range around the beta mean assuming gaussian distribtuion
bool inElectronRange(const double beta, const int runNum, const double sigmasAway = 3.0);
bool inMuonRange(const double beta, const int runNum, const double sigmasAway = 3.0);
bool inPionRange(const double beta, const int runNum, const double sigmasAway = 3.0);
//==	=

bool isCsI(const int runNum);

int getV1730waveformLength(const int runNum);

std::unordered_map<std::string, double> getGeantFileSimParamers(const std::string &filename);

static inline int getGeantFilePdg(const std::string &filename) {
  return std::lround(getGeantFileSimParamers(filename).at("primaryPdg"));
}

static inline double getGeantFilePrimaryMeanMomentum(const std::string &filename) {
  return getGeantFileSimParamers(filename).at("primaryMeanMomentum");
}

static inline double getGeantFileMomentumResolution(const std::string &filename) {
  return getGeantFileSimParamers(filename).at("momentumResolution");
}

static inline int getRunNumAccordingToMomentum(const double momentum) {
  if (std::abs(momentum - 167) < 10)
    return 571;
  if (std::abs(momentum - 145) < 10)
    return 599;
  else if (std::abs(momentum - 124) < 10)
    return 600;
  else
    return 1;
}

// Convert ADC counts to MeV. Depends on channel
inline double adc2mev(const double adc, const int channel, const Crystal crystal = Crystal::CsI_Tl_Belle) {
  return adc * c_crystal2_adc2mev.at(crystal).at(channel);
}

// Check if seagate mounted and return it.
// If not, return local HD.
inline std::string getTestbeamDir() {
  return opendir("/home/hershen/Seagate/originalMidasFiles") ? "/home/hershen/Seagate"
                                                             : "/home/hershen/PhD/Testbeam2015/midasFiles";
}

template <class... Tdrawable>
std::unique_ptr<TCanvas> drawWaveforms(const std::string &saveFilename, const bool waitDoubleClick,
                                       Tdrawable &... drawableObjects) {
  static int iDummy = 0;
  std::unique_ptr<TCanvas> canvas(new TCanvas(("canvas_" + std::to_string(iDummy)).data(), "", 0, 0, 1200, 900));

  drawTObjects(drawableObjects...);

  if (not saveFilename.empty()) {
    myFuncs::mySaveCanvas(canvas.get(), saveFilename);
  }

  if (waitDoubleClick) {
    myFuncs::waitForDoubleClick();
  }

  return canvas;
}

//Return a function which describes the pedestal as a function of event number
TF1 getPedestalFitFunction(const int runNum, const int channel);

} // namespace testbeam
} // namespace myFuncs
