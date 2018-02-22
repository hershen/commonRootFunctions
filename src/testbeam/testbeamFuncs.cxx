#include "testbeam/testbeamFuncs.h"

// Root
#include "TChain.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPad.h"

// Mine
#include "fileFuncs.h"
#include "generalFuncs.h"
#include "histFuncs.h" //for PaveText
#include "testbeam/RunDB.h"

namespace myFuncs {
namespace testbeam {

// The rightmost x position of the run + channel + event string
const double runChannelEventXpos = 0.6;

//-----------------------------------------------------------
// getNominalMomentumString
//-----------------------------------------------------------
std::string getNominalMomentumString(const int runNum) {
  int momentum = std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()); // Convert to int to get rid of decimal place
  return std::to_string(momentum) + " MeV/c";
}

//-----------------------------------------------------------
// getMeasuredMomentumString
//-----------------------------------------------------------
std::string getMeasuredMomentumString(const int runNum) {
  int momentum = std::lround(RunDB::instance()[runNum].getMeasuredBeamMomentum()); // Convert to int to get rid of decimal place
  return std::to_string(momentum) + " MeV/c";
}
//-----------------------------------------------------------
// getCrystalString
//-----------------------------------------------------------
std::string getCrystalString(const int runNum) {

  Crystal crystal = RunDB::instance()[runNum].getCrystal();
  if (crystal == Crystal::CsI_Tl_Belle)
    return "CsI(Tl) PIN";
  else if (crystal == Crystal::CsI_Tl_Babar)
    return "CsI(Tl) PP";
  else if (crystal == Crystal::CsI_Ukrainian)
    return "CsI AMCRYS PP";
  else if (crystal == Crystal::CsI_Chinese)
    return "CsI SICCAS PP";
  else
    return "Unknown crystal";
}

//-----------------------------------------------------------
// getSourceDistanceString
//-----------------------------------------------------------
std::string getSourceDistanceString(const int runNum) {
  int distance = std::lround(RunDB::instance()[runNum].getSourceDistance()); // Convert to int to get rid of decimal place
  if (distance < 0)
    return "no ^{60}Co";
  return "^{60}Co " + std::to_string(distance) + " cm away";
}

std::vector<std::string> getFilesRelatedToRun(const std::string pathToFiles, const int runNum, const std::string extension) {
  const std::vector<std::string> allFilenames = myFuncs::getFilesEndingWith(pathToFiles, extension);
  std::vector<std::string> filenames;
  // Find relevant files
  const std::string runString = std::to_string(runNum);
  for (const auto filename : allFilenames)
    if (filename.find(runString) != std::string::npos)
      filenames.push_back(filename);

  return filenames;
}

std::unique_ptr<TGraphErrors> getWaveformGraph(const std::vector<double>& voltage, const std::vector<double>& errors) {

  // Create time vector
  std::vector<double> times;
  times.reserve(voltage.size());
  for (size_t i = 0; i < voltage.size(); ++i)
    times.push_back(i * 2.0);

  std::unique_ptr<TGraphErrors> graph;

  if (errors.size() == 0) // no errors
    graph = std::unique_ptr<TGraphErrors>(new TGraphErrors(voltage.size(), times.data(), voltage.data(), nullptr, nullptr));
  else // with errors
    graph = std::unique_ptr<TGraphErrors>(new TGraphErrors(voltage.size(), times.data(), voltage.data(), nullptr, errors.data()));

  graph->SetTitle(";Time (ns); Voltage (ADC counts)");
  graph->SetMarkerSize(0.5);

  return graph;
}

std::unique_ptr<myFuncs::PaveText> getRunChannelEventPaveText(const int runNum, const int channelNum, const int eventNum) {
  std::unique_ptr<myFuncs::PaveText> pt(new myFuncs::PaveText(runChannelEventXpos));
  pt->AddText(
      (getRunParamsString(runNum) + ", event " + std::to_string(eventNum) + ", channel " + std::to_string(channelNum)).data());

  return pt;
}

bool inElectronRange(const double beta, const int runNum, const double sigmasAway) {
  const double sigma = electronBetaSigma.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
  return std::abs(1.0 - beta) < sigmasAway * sigma;
}

bool inMuonRange(const double beta, const int runNum, const double sigmasAway) {
  const double mean = muonBetaMean.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
  const double sigma = muonBetaSigma.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
  return std::abs(mean - beta) < sigmasAway * sigma;
}

bool inPionRange(const double beta, const int runNum, const double sigmasAway) {
  const double mean = pionBetaMean.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
  const double sigma = pionBetaSigma.at(std::lround(RunDB::instance()[runNum].getNominalBeamMomentum()));
  return std::abs(mean - beta) < sigmasAway * sigma;
}

bool isElectron(const double beta, const int runNum) {
  return inElectronRange(beta, runNum, 5.0) and not inMuonRange(beta, runNum, 4.0);
}
bool isMuon(const double beta, const int runNum) {
  return inMuonRange(beta, runNum, 4.0) and not inPionRange(beta, runNum, 4.0) and not inElectronRange(beta, runNum, 5.0);
}
bool isPion(const double beta, const int runNum) { return inPionRange(beta, runNum, 4.0) and not inMuonRange(beta, runNum, 4.0); }

int getV1730waveformLength(const int runNum) { return isCsI(runNum) ? 5000 : 17500; }

bool isCsI(const int runNum) {
  const auto crystal = RunDB::instance()[runNum].getCrystal();
  if (crystal == Crystal::CsI_Ukrainian or crystal == Crystal::CsI_Chinese)
    return true;
  else
    return false;
}

std::unordered_map<std::string, double> getGeantFileSimParamers(const std::string& filename) {

  const std::vector<std::string> branches{"beamParticle_pdgID", "nominalParticleMomentum_Mev_c", "particleMomentumResolution"};
  int pdg = 0;                    // Initialized so clang doesn't complain
  double primaryMeanMomentum = 0; // Initialized so clang doesn't complain
  double momentumResolution = 0;  // Initialized so clang doesn't complain
  const std::vector<void*> pointers{&pdg, &primaryMeanMomentum, &momentumResolution};

  std::unique_ptr<TChain> chain(myFuncs::openChain_setBranch(filename, "simParameters", branches, pointers));

  chain->GetEntry(0);

  std::unordered_map<std::string, double> map;
  map["primaryPdg"] = pdg;
  map["primaryMeanMomentum"] = primaryMeanMomentum;
  map["momentumResolution"] = momentumResolution;

  return map;
}

TF1 getPedestalFitFunction(const int runNum, const int channel) {
  const pedestalFitParams fitParams = RunDB::instance()[runNum].getPedestalFitParams(channel);
  TF1 func(("pedesalFunc_" + std::to_string(runNum) + "_" + std::to_string(channel)).c_str(), "pol2", 0, 1e6);
  func.FixParameter(0, fitParams.p0);
  func.FixParameter(1, fitParams.p1);
  func.FixParameter(2, fitParams.p2);
  return func;
}

std::vector<double> getTimes(const size_t size, const double dt) {
  std::vector<double> times;
  times.reserve(size);

  for (uint iSample = 0; iSample < size; ++iSample)
    times.push_back(iSample * dt);

  return times;
}

// Load filtering parameters
FilterInfo loadFilteringParams(const std::string& filterParamsFilename) {
  FilterInfo filterInfo;
  const auto params(myFuncs::getParams(filterParamsFilename));

  // CR-RC4
  if (params.find("CR_RC4_timeConstant_ns") != params.end()) {
    filterInfo.CR_RC4_timeConstant_ns = std::atof(params.at("CR_RC4_timeConstant_ns").c_str());
    filterInfo.figureText += (boost::format("CR-RC^{4} (%.0f ns)") % filterInfo.CR_RC4_timeConstant_ns).str();
    filterInfo.fileText += (boost::format("CR-RC4_%.0fns") % filterInfo.CR_RC4_timeConstant_ns).str();
  } else {
    filterInfo.CR_RC4_timeConstant_ns = 0.0;
  }

  // HPF
  if (params.find("HPF_coefficientFile") != params.end()) {
    filterInfo.HPF_CoeficientFile = params.at("HPF_coefficientFile");
    filterInfo.HPF_filterOrder = std::atof(params.at("HPF_filterOrder").c_str());
    filterInfo.HPFcutoffFreq_GHz = std::atof(params.at("HPF_cutOff").c_str());
    filterInfo.fir1Coeficients = myFuncs::changeStringsToDoubles(myFuncs::readFile(filterInfo.HPF_CoeficientFile));
  } else {
    filterInfo.HPF_CoeficientFile = "";
    filterInfo.HPF_filterOrder = 0.0;
    filterInfo.HPFcutoffFreq_GHz = 0.0;
  }

  // reductionFactor
  if (params.find("reductionFactor") != params.end()) {
    filterInfo.reductionFactor = std::atoi(params.at("reductionFactor").c_str());
  } else {
    filterInfo.reductionFactor = 1;
  }
  if (filterInfo.reductionFactor > 1) {
    filterInfo.fileText += "_reduction" + std::to_string(filterInfo.reductionFactor);
  }

  return filterInfo;
}
FilterInfo loadFilteringParamsFromTTree(const std::string& filename) {
  TFile file(filename.c_str(),"READ");
  TTree* paramsTree = (TTree*)file.Get("params");
  char filterCharStar[100];
  FilterInfo filterInfo;
  paramsTree->SetBranchAddress("filter", &filterCharStar);
  paramsTree->SetBranchAddress("reductionFactor", &filterInfo.reductionFactor);
  paramsTree->SetBranchAddress("CR_RC4_timeConstant_ns", &filterInfo.CR_RC4_timeConstant_ns);
  paramsTree->SetBranchAddress("HPForder", &filterInfo.HPF_filterOrder);
  paramsTree->SetBranchAddress("HPFcutoffFreq_GHz", &filterInfo.HPFcutoffFreq_GHz);
  paramsTree->GetEntry(0);

  filterInfo.fileText = filterCharStar;
  return filterInfo;
}
} // namespace testbeam

} // namespace myFuncs
