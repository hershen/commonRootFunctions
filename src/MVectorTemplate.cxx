#include "MVectorTemplate.h"

#include <iostream>
#include <memory>

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"
#include "TVirtualFitter.h"

#include "fileFuncs.h"
#include "mathFuncs.h"

#define cline std::cout << "line = " << __LINE__ << std::endl;

// temperary
#include "TGraph.h"

// For TMiuit2
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"

namespace myFuncs {

MVectorTemplate::MVectorTemplate()
    : m_tF1(TF1()), m_treeName("vectorTemplate"), m_templateValuesBranchName("m_templateValues"), m_dxBranchName("rDx"),
      m_numAveragedFuncsBranchName("rNumAveragedFuncs"), m_peakIdxBranchName("m_peakIdx"), m_numAveragedFuncs(0), m_peakIdx(0),
      m_dx(0.), m_minAmplitudeLimit(0), m_maxAmplitudeLimit(-1), m_minChi2_NdfLimit(-1), m_maxChi2_NdfLimit(3),
      m_minXshiftLimit(0), m_maxXshiftLimit(0), m_minPedestalLimit(0), m_maxPedestalLimit(0), m_xValueOfFirstTemplateEntry(0),
      m_doubleNumbersEqualThershold(1e-20), m_debugLevel(0), m_averageMode(keepWeight), m_amplitudeFitEnabled(true),
      m_pedestalFitEnabled(true), m_xShiftFitEnabled(true), m_useAmplitudeLimits(false), m_useXshiftLimits(false),
      m_usePedestalLimits(false) {}

MVectorTemplate::MVectorTemplate(const std::vector<double> &newVec, const double newDx) : MVectorTemplate() {
  setDx(newDx);
  addVector(newVec, 0, m_tF1);
}

void MVectorTemplate::setTF1ParNames() {
  m_tF1.SetParName(0, "Amplitude");
  m_tF1.SetParName(1, "Pedestal");
  m_tF1.SetParName(2, "xShift");
}

void MVectorTemplate::resetTemplateValues() {
  m_templateValues.clear();
  m_numAveragedFuncs = 0;
  m_tF1 = TF1();
}

void MVectorTemplate::resetTemplateRange() {
  if (m_tF1.IsValid())
    m_tF1.SetRange(getXvalueOfFirstTemplateEntry(), getXvalueOfFirstTemplateEntry() + m_dx * (m_templateValues.size() - 1));
}

void MVectorTemplate::setDx(const double newDx) {
  m_dx = newDx;

  // Reset tempalte range only if there are template values
  if (m_templateValues.size())
    resetTemplateRange();
}

void MVectorTemplate::setXvalueOfFirstTemplateEntry(const double xValueOfFirstTemplateEntry) {
  m_xValueOfFirstTemplateEntry = xValueOfFirstTemplateEntry;
  resetTemplateRange();
}

double MVectorTemplate::calcSimplePedestal(const std::vector<double> &vector, const double percentage) const {

  const long numElementsToAverage = std::distance(vector.begin(), vector.begin() + (vector.size() * percentage));

  if (numElementsToAverage <= 1)
    return vector[0];

  return std::accumulate(vector.begin(), vector.begin() + numElementsToAverage, 0.0) / numElementsToAverage;
}

double MVectorTemplate::getAmplitudeGuess(const std::vector<double> &vector, const double pedestal) const {
  // Search for maximum around m_peakIdx
  const int elementsBeforeAfter = vector.size() * 0.1;

  auto peakIdx = getEffPeakIdx();
  if (peakIdx > vector.size() - 1)
    peakIdx = vector.size() - 1;

  long firstElement = peakIdx - elementsBeforeAfter;
  if (firstElement < 0)
    firstElement = 0;

  long lastElement = peakIdx + elementsBeforeAfter;
  if (lastElement > static_cast<long>(vector.size()) - 1)
    lastElement = vector.size() - 1;

  const double maxElement = *std::max_element(vector.begin() + firstElement, vector.begin() + lastElement);

  if (getDebugLevel() > 15)
    std::cout << "MVectorTemplate::getAmplitudeGuess: maxElement = " << maxElement << std::endl;

  // Amplitude initial guess is the value of vector where the template maximum is. This will not work for fast varying functions,
  // or if the xShift is large
  double amplitudeGuess = (maxElement - pedestal);

  // Make sure amplitude guess is not out of limits
  if (m_useAmplitudeLimits) {
    amplitudeGuess = std::min(amplitudeGuess, m_maxAmplitudeLimit);
    amplitudeGuess = std::max(amplitudeGuess, m_minAmplitudeLimit);
  }

  return amplitudeGuess;
}

size_t MVectorTemplate::getEffPeakIdx() const {
  long effPeakIdx = m_peakIdx + getXvalueOfFirstTemplateEntry() / m_dx;
  if (effPeakIdx < 0)
    return 0;
  return effPeakIdx;
}

void MVectorTemplate::addFirstVector(const std::vector<double> &newVector) {
  // ------------------------------------------------------------------
  // Take pedestal as first 5% of entries
  // ------------------------------------------------------------------
  const double pedestal = calcSimplePedestal(newVector);

  // ------------------------------------------------------------------
  // Search for maximum (and what idx it is at - for later guessing the time shift) of newVector so peak can be normalized
  // ------------------------------------------------------------------
  const auto maxElementIt = std::max_element(newVector.begin(), newVector.end());
  const double maximum = *maxElementIt;
  m_peakIdx = std::distance(newVector.begin(), maxElementIt);

  const double maxMinusPedestal = maximum - pedestal;
  // ------------------------------------------------------------------
  // Populate m_templateValues, normalizing and subtracting pedestal
  // ------------------------------------------------------------------
  m_templateValues.clear();
  m_templateValues.reserve(newVector.size());

  const double ov_maxMinusPedestal = 1.0 / maxMinusPedestal;

  for (size_t idx = 0; idx < newVector.size(); ++idx)
    m_templateValues.push_back((newVector[idx] - pedestal) * ov_maxMinusPedestal);

  if (m_debugLevel > 20) {
    std::cout << "MVectorTemplate::addFirstVector : pedestal = " << pedestal << ", maximum = " << maximum
              << ", maxMinusPedestal = " << maxMinusPedestal << std::endl;
    for (size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
      std::cout << "MVectorTemplate::addFirstVector : temp[" << iTemplateIdx << "] = " << m_templateValues[iTemplateIdx]
                << std::endl;
  }

  // ------------------------------------------------------------------
  // Create TF1
  // ------------------------------------------------------------------
  m_tF1 = TF1("templateTF1", this, &MVectorTemplate::TF1Eval, 0, m_dx * (m_templateValues.size() - 1), 3, "MVectorTemplate",
              "TF1Eval");
  setTF1ParNames();
  m_tF1.SetParameters(maxMinusPedestal, pedestal, 0.0);

  std::cout << "MVectorTemplate::addFirstVector : remember to setXvalueOfFirstTemplateEntry!" << std::endl;
}

bool MVectorTemplate::goodFit(const TFitResultPtr &fitResult) const {
  const double fittedAmplitude = fitResult->Parameter(0);
  const double fittedPedestal = fitResult->Parameter(1);
  const double fittedXshift = fitResult->Parameter(2);

  const int fitStatus = fitResult;
  if (fitStatus) {
    if (fitStatus != 4000) {
      if (m_debugLevel > 0)
        std::cerr << "MVectorTemplate::goodFit : fit failed, fit status = " << fitStatus << ". Not adding vector!" << std::endl;
      return false;
    } else if (m_debugLevel > 0)
      std::cout << "MVectorTemplate::goodFit : MINOS falied, fit status = " << fitStatus << ". Vector MIGHT still be added!!!"
                << std::endl;
  }

  const double chi2_ndf = fitResult->Chi2() / fitResult->Ndf();

  if (chi2_ndf < m_minChi2_NdfLimit || chi2_ndf > m_maxChi2_NdfLimit) {
    if (m_debugLevel > 0)
      std::cerr << "MVectorTemplate::goodFit : fit failed, Chi2/NDF of fit = " << chi2_ndf << " is out of limits. limits are ("
                << m_minChi2_NdfLimit << "," << m_maxChi2_NdfLimit << "). Not adding vector!" << std::endl;
    return false;
  }

  if (m_debugLevel > 10)
    std::cout << "MVectorTemplate::goodFit : fittedAmplitude = " << fittedAmplitude << ", fittedPedestal = " << fittedPedestal
              << ", fittedXshift = " << fittedXshift << ", chi2 / NDF = " << chi2_ndf << std::endl;

  if (m_useAmplitudeLimits and (std::abs(fittedAmplitude - m_minAmplitudeLimit) < m_doubleNumbersEqualThershold ||
                                std::abs(fittedAmplitude - m_maxAmplitudeLimit) < m_doubleNumbersEqualThershold)) {
    if (m_debugLevel > 10)
      std::cerr << "MVectorTemplate::goodFit : fit failed, reached amplitude limit. fittedAmplitude = " << fittedAmplitude
                << ". limits are (" << m_minAmplitudeLimit << "," << m_maxAmplitudeLimit << "). Not adding vector!" << std::endl;
    return false;
  }

  if (m_useXshiftLimits and (std::abs(fittedXshift - m_minXshiftLimit) < m_doubleNumbersEqualThershold ||
                             std::abs(fittedXshift - m_maxXshiftLimit) < m_doubleNumbersEqualThershold)) {
    if (m_debugLevel > 10)
      std::cerr << "MVectorTemplate::goodFit : fit failed, reached xShift limit. fittedXshift = " << fittedXshift
                << ". limits are (" << m_minXshiftLimit << "," << m_maxXshiftLimit << "). Not adding vector!" << std::endl;
    return false;
  }

  if (m_usePedestalLimits and (std::abs(fittedPedestal - m_minPedestalLimit) < m_doubleNumbersEqualThershold ||
                               std::abs(fittedPedestal - m_maxPedestalLimit) < m_doubleNumbersEqualThershold)) {
    if (m_debugLevel > 10)
      std::cerr << "MVectorTemplate::goodFit : fit failed, reached pedestal limit. fittedPedestal = " << fittedPedestal
                << ". limits are (" << m_minPedestalLimit << "," << m_maxPedestalLimit << "). Not adding vector!" << std::endl;
    return false;
  }

  return true;
}

TFitResult MVectorTemplate::fitFunctionToVector(const std::vector<double> &vector, const double std, TF1 &function) {

  // Histrogram with entries of vector - used to fit the template
  // Center of each bin corresponds to correct x value (that's why there's a -m_dx/2. term)
  TH1D fitHist("fitHist", "fitHist", vector.size(), -m_dx / 2.0, vector.size() * m_dx - m_dx / 2.0);
  for (size_t idx = 1; idx <= vector.size(); ++idx) {
    fitHist.SetBinContent(idx, vector[idx - 1]);
    fitHist.SetBinError(idx, std);
  }

  // ------------------------------------------------------------------
  // Initial guesses for fit
  // ------------------------------------------------------------------
  const double pedestalGuess = isPedestalFitEnabled() ? calcSimplePedestal(vector) : function.GetParameter(1);

  const double amplitudeGuess = isAmplitudeFitEnabled() ? getAmplitudeGuess(vector, pedestalGuess) : function.GetParameter(0);

  const double xShiftGuess = isXshiftFitEnabled() ? 0.0 : function.GetParameter(2);

  if (m_debugLevel > 15)
    std::cout << "MVectorTemplate::fitFunctionToVector : amplitudeGuess = " << amplitudeGuess
              << ", pedestalGuess = " << pedestalGuess << ", xShiftGuess = " << xShiftGuess << std::endl;

  setTF1Parameters(function, amplitudeGuess, pedestalGuess, xShiftGuess); // Amplitude, pedestal, xShift

  // ------------------------------------------------------------------
  // Fit!
  // ------------------------------------------------------------------

  // First fit is used to find the xShift
  TFitResultPtr fitResult = fitTemplate(fitHist, function);

  if (!goodFit(fitResult))
    return TFitResult();

  // We then Change the function range so it takes into account the xShift. Otherwise, the ends of the function are not correctly
  // represented and the averaging is not good.
  double minRange, maxRange;
  function.GetRange(minRange, maxRange);
  if (getDebugLevel() > 15)
    std::cout << "MVectorTemplate::fitFunctionToVector : About to fit in the range = (" << minRange + fitResult->Parameter(2)
              << "," << maxRange + fitResult->Parameter(2) << ")" << std::endl;

  function.SetRange(minRange + fitResult->Parameter(2), maxRange + fitResult->Parameter(2));

  // Redo the fit
  fitResult = fitTemplate(fitHist, function);

  // Set back the correct function range
  function.SetRange(minRange, maxRange);

  if (!goodFit(fitResult))
    return TFitResult();

  return *fitResult;
}

std::vector<double> MVectorTemplate::getTimeAlignedVector(const std::vector<double> &newVector, const double xOfNewVector_0) {

  std::vector<double> alignedVector;
  alignedVector.reserve(getTemplateSize());

  for (size_t iTemplateIdx = 0; iTemplateIdx < getTemplateSize(); ++iTemplateIdx) {

    // x value corresponding to iTemplateIdx
    const double x = getXvalueOfFirstTemplateEntry() + (iTemplateIdx * m_dx);

    // idx in newVector which is the first one before x
    const size_t idx = (x - xOfNewVector_0) / m_dx;

    assert(idx < newVector.size());

    // Linear interpolation
    const double x1 = xOfNewVector_0 + static_cast<double>(idx) * m_dx;
    const double x2 = xOfNewVector_0 + static_cast<double>(idx + 1) * m_dx;
    const double y1 = newVector[idx];
    const double y2 = newVector[idx + 1];

    double timeAlignedVectorValue = linearInterpolate(x1, x2, y1, y2, x);

    alignedVector.push_back(timeAlignedVectorValue);
  }

  return alignedVector;
}

void MVectorTemplate::clipTemplateEnds(const double xOfNewVector_0, const size_t newVectorSize) {

  // Find first and last relevant template indices

  // First Idx

  // Has to be max(double, double) because if both are size_t, and xOfNewVector_0 - getXvalueOfFirstTemplateEntry() < 0, it
  // becomes ~infinity
  const size_t firstTemplateIdx = std::max(0.0, (xOfNewVector_0 - getXvalueOfFirstTemplateEntry()) / m_dx + 1);

  // Last idx
  const double timeOfLastTemplateEntry = getXvalueOfFirstTemplateEntry() + (m_templateValues.size() - 1) * m_dx;
  const double timeOfLastNewVectorValue = xOfNewVector_0 + (newVectorSize - 1) * m_dx;

  if (getDebugLevel() > 20)
    std::cout << "MVectorTemplate::clipTemplateEnds : timeOfLastTemplateEntry = " << timeOfLastTemplateEntry
              << ", timeOfLastNewVectorValue = " << timeOfLastNewVectorValue << std::endl;

  const size_t lastTemplateIdx =
      std::min(m_templateValues.size() - 1,
               static_cast<size_t>(m_templateValues.size() - 1 - (timeOfLastTemplateEntry - timeOfLastNewVectorValue) / m_dx));

  if (getDebugLevel() > 10)
    std::cout << "MVectorTemplate::clipTemplateEnds : firstTemplateIdx = " << firstTemplateIdx
              << ", lastTemplateIdx = " << lastTemplateIdx << std::endl;

  // clip End of template
  if (lastTemplateIdx + 1 < m_templateValues.size())
    m_templateValues.resize(lastTemplateIdx + 1); // lastTemplateIdx + 1 == size of new template vector;

  // clip Start of template
  if (firstTemplateIdx > 0) {
    m_templateValues.erase(m_templateValues.begin(),
                           m_templateValues.begin() +
                               firstTemplateIdx); // erase first int(firstTemplateIdx) elements in m_templateValues.
    m_xValueOfFirstTemplateEntry += firstTemplateIdx * m_dx;
  }

  resetTemplateRange();
}
TFitResult MVectorTemplate::addVector(const std::vector<double> &newVector, const double std, TF1 &function) {
  // TODO implement with TMinuit
  // TODO average using a weighted average
  // ------------------------------------------------------------------
  //   //Sanity checks
  // ------------------------------------------------------------------
  if (newVector.empty()) {
    std::cerr << "MVectorTemplate::addVector : Error - trying to add empty vector. Not adding. " << std::endl;
    return -1;
  }

  if (std < 0) {
    std::cerr << "MVectorTemplate::addVector : Error - negative standard deviation (" << std << "). Not adding. " << std::endl;
    return -2;
  }

  // ------------------------------------------------------------------
  // Set first vector as template, normalize, etc.
  // ------------------------------------------------------------------
  if (m_numAveragedFuncs == 0) {

    addFirstVector(newVector);
    ++m_numAveragedFuncs;
    return TFitResult();
  } // First vector added

  else { // from second vector

    if (getDebugLevel() > 50) {
      std::cout << "MVectorTemplate::addVector : template values before adding vector:" << std::endl;

      for (size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
        std::cout << "template[" << iTemplateIdx << "] = " << m_templateValues[iTemplateIdx] << "\n";

      std::cout << std::endl;
    }

    // Fit function
    TFitResult fitResult = fitFunctionToVector(newVector, std, function);

    if (!fitResult.IsValid())
      return fitResult;

    // ------------------------------------------------------------------
    // Positive fittedXshift means newVector[0] is forward in time relative to m_templateValues[0], i.e. :
    // In order to align the newVector, newVector[0] should correspond to time (-fittedXshift)
    // not sure if true:
    // newVector[idx] = m_templateValues[idx - fittedXshift/m_dx]
    // ------------------------------------------------------------------

    const double fittedXshift = fitResult.Parameter(2);
    const double xOfNewVector_0 = /*getXvalueOfFirstTemplateEntry()*/ -fittedXshift;

    if (getDebugLevel() > 20)
      std::cout << "MVectorTemplate::addVector : xOfNewVector_0 = " << xOfNewVector_0 << std::endl;

    // 		std::vector<double> times;
    // 		times.reserve(newVector.size());
    // 		for(int idx = 0; idx < newVector.size(); ++idx) {
    // 			times.push_back(idx * 2.0);
    // 		}
    // 		TCanvas c("c","",0,0,1200,900);
    // 		c.Divide(1,2);
    // 		c.cd(1);
    // 		TGraph newVectorGraph(times.size(), times.data(), newVector.data());
    // 		newVectorGraph.Draw("AP");
    // 		function.Draw("Same");
    // 		c.cd(2);
    //
    // 		std::vector<double> templateTimes;
    // 		templateTimes.reserve(getTemplateValues().size());
    // 		for(int idx = 0; idx < getTemplateValues().size(); ++idx)
    // 			templateTimes.push_back(effXofFirstTemplateEntry + idx*2.0);
    //
    // 		TGraph templateGraph(templateTimes.size(), templateTimes.data(), getTemplateValues().data());
    // 		templateGraph.Draw("AP");
    // 		gPad->WaitPrimitive();

    // fittedXshift rescales the time of newVector. newVector[0] corresponds to time xOfNewVector_0, newVector[1] corresponds to
    // time (xOfNewVector_0 + m_dx)...

    // Clip first and end of template, if they don't overlap in time with the newVector. We want to keep only the parts that can
    // be averaged.
    clipTemplateEnds(xOfNewVector_0, newVector.size());

    assert(getXvalueOfFirstTemplateEntry() >= xOfNewVector_0);

    // newAlignedVector has it's values time aligned to template values, time of 0'th element of newTimeAlignedVector has same
    // time as templateVlaues[0]
    std::vector<double> newTimeAlignedVector = getTimeAlignedVector(newVector, xOfNewVector_0);

    /*
if(getXvalueOfFirstTemplateEntry >= xOfNewVector_0 )
{
template:                      Idx0---------Idx1---------Idx2
Idx(size-4)-------------Idx(size-3)-------------Idx(size-2)-------------Idx(size-1) newVector:
Idx0---------Idx1---------Idx2---------Idx3     Idx(size-3)-------------Idx(size-2)-------------Idx(size-1) In this example both
vectors are same length and xOfNewVector_0 = getXvalueOfFirstTemplateEntry - 1.5m_dx firstTemplateIdx = 0, lastTemplateIdx =
newVector.size()-3
    }
else
{
template:  Idx0--------Idx1--------Idx2          Idx(size-2)-------------Idx(size-1)
newVector:                   Idx0--------Idx1                Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
In this example both vectors are same length and xOfNewVector_0 = getXvalueOfFirstTemplateEntry + 1.5m_dx
firstTemplateIdx = 2, lastTemplateIdx = m_templateValues.size()-1
}*/

    if (getDebugLevel() > 10)
      std::cout << "MVectorTemplate::addVector : About to average in the range (" << getXvalueOfFirstTemplateEntry() << ","
                << getXvalueOfFirstTemplateEntry() + (getTemplateSize() - 1) * m_dx << ")" << std::endl;

    // ------------------------------------------------------------------
    // Average newVector with template
    // ------------------------------------------------------------------

    const double ov_numAveragedFuncsP1 = 1.0 / static_cast<double>(m_numAveragedFuncs + 1);

    // Average
    double peakVal = -DBL_MAX;

    const double fittedAmplitude = fitResult.Parameter(0);
    const double fittedPedestal = fitResult.Parameter(1);

    for (size_t iTemplateIdx = 0; iTemplateIdx < getTemplateSize(); ++iTemplateIdx) {
      double newVectorValue = newTimeAlignedVector[iTemplateIdx];

      if (getAverageMode() == noWeight) {

        // Normalize newVectorValue so it can be added to template values.
        newVectorValue = (newVectorValue - fittedPedestal) / fittedAmplitude;

        // average vectors
        m_templateValues[iTemplateIdx] += newVectorValue;
      } else if (getAverageMode() == keepWeight) { // This assumes the result will later be normalized so peak is at 1.
        const double templateValue = fittedPedestal + fittedAmplitude * m_templateValues[iTemplateIdx];
        m_templateValues[iTemplateIdx] =
            (templateValue * static_cast<double>(m_numAveragedFuncs) + newVectorValue) * ov_numAveragedFuncsP1;
      }

      if (getDebugLevel() > 30) {
        std::cout << "MVectorTemplate::addVector : Averaging template iTemplateIdx " << iTemplateIdx
                  << ", time = " << getXvalueOfFirstTemplateEntry() + iTemplateIdx * m_dx
                  << ". newVectorValue = " << newVectorValue
                  << ", m_templateValues[iTemplateIdx] after averaging = " << m_templateValues[iTemplateIdx] << std::endl;
      }

      // Store maximum value (used to normalize later) and its index
      if (m_templateValues[iTemplateIdx] > peakVal) {
        peakVal = m_templateValues[iTemplateIdx];
        m_peakIdx = iTemplateIdx;
      }
    }

    // ------------------------------------------------------------------
    // re-Normalize vector - keep peak == 1
    // Make sure pedestal is zero
    // ------------------------------------------------------------------
    if (getAverageMode() == keepWeight)
      normalizeAndZeroTemplateValues();

    if (getDebugLevel() > 20) {
      std::cout << "template values after adding vector:" << std::endl;
      for (size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
        std::cout << "template[" << iTemplateIdx << "] = " << m_templateValues[iTemplateIdx] << "\n";

      std::cout << std::endl;
    }

    ++m_numAveragedFuncs;

    return fitResult;
  } // add, from second vector

  return TFitResult();
}

TFitResultPtr MVectorTemplate::fitTemplate(TH1D &fitHist, TF1 &function) {

  // Attempt to fit at current x.

  TFitResultPtr fitResult = fitHist.Fit(&function, "QSNR");
  //"Q" Quiet mode (minimum printing)
  //"M" More. Improve fit results. It uses the IMPROVE command of TMinuit (see TMinuit::mnimpr). This algorithm attempts to
  // improve the found local minimum by searching for a better one. "N" Do not store the graphics function, do not draw "S" The
  // result of the fit is returned in the TFitResultPtr (see below Access to the Fit Result)

  if (fitGood(fitResult) or !isXshiftFitEnabled())
    return fitResult;

  // Fit failed. Move xShift and retry
  double xShiftGuess = m_dx * m_templateValues.size() * 0.05;
  if (m_debugLevel >= 20)
    std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
  function.SetParameter(2, xShiftGuess);
  fitResult = fitHist.Fit(&function, "QSNR");
  if (fitGood(fitResult))
    return fitResult;

  // move xShift and retry
  xShiftGuess = -m_dx * m_templateValues.size() * 0.05;
  if (m_debugLevel >= 20)
    std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
  function.SetParameter(2, xShiftGuess);
  fitResult = fitHist.Fit(&function, "QSNR");
  if (fitGood(fitResult))
    return fitResult;

  xShiftGuess = m_dx * m_templateValues.size() * 0.1;
  if (m_debugLevel >= 20)
    std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
  function.SetParameter(2, xShiftGuess);
  fitResult = fitHist.Fit(&function, "QSNR");
  if (fitGood(fitResult))
    return fitResult;

  // move xShift and retry
  xShiftGuess = -m_dx * m_templateValues.size() * 0.1;
  if (m_debugLevel >= 20)
    std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
  function.SetParameter(2, xShiftGuess);
  fitResult = fitHist.Fit(&function, "QSNR");

  // Return fitResult in any case
  return fitResult;
}

bool MVectorTemplate::fitGood(const TFitResultPtr &fitResult) const {

  const int fitStatus = fitResult;

  // Fit status
  if (fitStatus != 0 and fitStatus != 4000) // 4000 is MINOS problem
    return false;

  // chi2/ndf
  const double chi2_ndf = fitResult->Chi2() / fitResult->Ndf();
  if (chi2_ndf < m_minChi2_NdfLimit or chi2_ndf > m_maxChi2_NdfLimit)
    return false;

  return true;
}

void MVectorTemplate::enableAmplitudeFit(const bool enableAmplitudeFit) {
  m_amplitudeFitEnabled = enableAmplitudeFit;
  setTF1Parameters();
}

void MVectorTemplate::enablePedestalFit(const bool enablePedestalFit) {
  m_pedestalFitEnabled = enablePedestalFit;
  setTF1Parameters();
}

void MVectorTemplate::enableXshiftFit(const bool enableXshiftFit) {
  m_xShiftFitEnabled = enableXshiftFit;
  setTF1Parameters();
}

double MVectorTemplate::TF1Eval(double *var, double *params) {
  //   double amplitude = params[0];
  //   double pedestal = params[1];
  //   double xShift = params[2];

  const double effectiveX = var[0] - params[2] - getXvalueOfFirstTemplateEntry();

  if (effectiveX <= 0)
    return params[1]; // + params[0]*m_templateValues[0];
  if (effectiveX >= m_dx * static_cast<double>(m_templateValues.size() - 1))
    return params[1] + params[0] * m_templateValues.back();

  const int idx = static_cast<int>(effectiveX / m_dx);

  // Linear interpolation
  const double x1 = static_cast<double>(idx) * m_dx;
  const double x2 = static_cast<double>(idx + 1) * m_dx;
  const double y1 = m_templateValues[idx];
  const double y2 = m_templateValues[idx + 1];

  // 	std::cout << "MVectorTemplate::TF1Eval :var[0] = " << var[0] << ",  x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 <<
  // ", y2 = " << y2 << ", effectiveX = " << effectiveX << ", TF1Eval = " << params[1] + params[0]*linearInterpolate(x1, x2, y1,
  // y2, effectiveX) <<  std::endl;

  return params[1] + params[0] * linearInterpolate(x1, x2, y1, y2, effectiveX);
}

void MVectorTemplate::setTF1Parameters(TF1 &function, const double amplitude, const double pedestal, const double xShift) {
  if (!function.IsValid())
    return;

  if (m_amplitudeFitEnabled) {
    function.SetParameter(0, amplitude);
    function.SetParLimits(0, m_minAmplitudeLimit, m_maxAmplitudeLimit);
  } else
    function.FixParameter(0, amplitude);

  if (m_pedestalFitEnabled) {
    function.SetParameter(1, pedestal);
    function.SetParLimits(1, m_minPedestalLimit, m_maxPedestalLimit);
  } else
    function.FixParameter(1, pedestal);

  if (m_xShiftFitEnabled) {
    function.SetParameter(2, xShift);
    function.SetParLimits(2, m_minXshiftLimit, m_maxXshiftLimit);
  } else
    function.FixParameter(2, xShift);
}

void MVectorTemplate::setTF1Parameters() {
  if (!m_tF1.IsValid())
    return;
  setTF1Parameters(m_tF1, m_tF1.GetParameter(0), m_tF1.GetParameter(1), m_tF1.GetParameter(2));
}

void MVectorTemplate::setAmplitudeLimits(const double newMinAmplitudeLimit, const double newMaxAmplitudeLimit) {
  if (std::abs(newMinAmplitudeLimit - newMaxAmplitudeLimit) > m_doubleNumbersEqualThershold)
    m_useAmplitudeLimits = true;
  else
    m_useAmplitudeLimits = false;
  m_minAmplitudeLimit = newMinAmplitudeLimit;
  m_maxAmplitudeLimit = newMaxAmplitudeLimit;
  setTF1Parameters();
}

void MVectorTemplate::setXshiftLimits(const double newMinXshiftLimit, const double newMaxXshiftLimit) {
  if (std::abs(newMinXshiftLimit - newMaxXshiftLimit) > m_doubleNumbersEqualThershold)
    m_useXshiftLimits = true;
  else
    m_useXshiftLimits = false;
  m_minXshiftLimit = newMinXshiftLimit;
  m_maxXshiftLimit = newMaxXshiftLimit;
  setTF1Parameters();
}

void MVectorTemplate::setPedestalLimits(const double newMinPedestalLimit, const double newMaxPedestalLimit) {
  if (std::abs(newMinPedestalLimit - newMaxPedestalLimit) > m_doubleNumbersEqualThershold)
    m_usePedestalLimits = true;
  else
    m_usePedestalLimits = false;
  m_minPedestalLimit = newMinPedestalLimit;
  m_maxPedestalLimit = newMaxPedestalLimit;
  setTF1Parameters();
}

int MVectorTemplate::saveTemplateToTFile(const std::string &fullFileName, const std::string &treeDescription) {

  std::cout << "MVectorTemplate::saveTemplateToTFile: Saving template to " << fullFileName << std::endl;

  // ------------------------------------------------------------------
  // Create TFile
  // ------------------------------------------------------------------
  TFile outputFile(fullFileName.data(), "RECREATE");

  if (outputFile.IsZombie()) {
    std::cerr << "MVectorTemplate::saveTemplateToTfile : could not open file " << fullFileName.data() << ". Aborting"
              << std::endl;
    return -1;
  }

  // ------------------------------------------------------------------
  // Create TTree
  // ------------------------------------------------------------------
  TTree tree(m_treeName.data(), treeDescription.data());

  // ------------------------------------------------------------------
  // Create pointers to variables that will be saved
  // ------------------------------------------------------------------
  std::vector<double> *templateValuesPointer = &m_templateValues;

  // ------------------------------------------------------------------
  // Create branches
  // ------------------------------------------------------------------
  tree.Branch(m_templateValuesBranchName.data(), "std::vector<double>", &templateValuesPointer);
  tree.Branch(m_dxBranchName.data(), &m_dx);
  tree.Branch(m_peakIdxBranchName.data(), &m_peakIdx, "m_peakIdx/l");
  tree.Branch(m_numAveragedFuncsBranchName.data(), &m_numAveragedFuncs, "m_numAveragedFuncs/l");
  tree.Branch("xValueOfFirstTemplateEntry", &m_xValueOfFirstTemplateEntry, "m_xValueOfFirstTemplateEntry/D");

  // ------------------------------------------------------------------
  // Save
  // ------------------------------------------------------------------
  tree.Fill();
  outputFile.Write();
  outputFile.Close();
  return 0;
}

int MVectorTemplate::loadTemplateFromTFile(const std::string &fullFileName) {

  std::vector<std::string> branchNamesV = {m_templateValuesBranchName.data(), m_dxBranchName.data(), m_peakIdxBranchName.data(),
                                           m_numAveragedFuncsBranchName.data()};
  std::vector<double> *m_templateValuesPointer = 0;

  std::vector<void *> pointerV = {&m_templateValuesPointer, &m_dx, &m_peakIdx, &m_numAveragedFuncs};

  std::unique_ptr<TChain> inputChain(
      myFuncs::openChain_setBranch(fullFileName.data(), m_treeName.data(), branchNamesV, pointerV));

  if (!inputChain) {
    std::cout << "MVectorTemplate::loadTemplateFromTFile : Error - Could not load file " << fullFileName.data() << ". Aborting"
              << std::endl;
    return -1;
  }

  if (inputChain->GetEntries() < 1) {
    std::cout << "MVectorTemplate::loadTemplateFromTFile : Error - No entries in file " << fullFileName.data() << ". Aborting"
              << std::endl;
    return -2;
  }

  // Add xValueOfFirstTemplateEntry branch. Earlier versions of MVectorTemplate didn't have it so it's addition is not automatic.
  myFuncs::setChainBranch(inputChain.get(), "xValueOfFirstTemplateEntry", &m_xValueOfFirstTemplateEntry);

  inputChain->GetEntry(0);

// Prevent clang static analyzer from thinking this derefrences null pointer - makes this line invisible to the analyzer
#ifndef __clang_analyzer__
  m_templateValues = *m_templateValuesPointer;
#endif

  delete m_templateValuesPointer;

  m_tF1 = TF1("templateTF1", this, &MVectorTemplate::TF1Eval, 0, m_dx * (m_templateValues.size() - 1), 3, "MVectorTemplate",
              "TF1Eval");
  resetTemplateRange();
  setTF1ParNames();
  return 0;
}

void MVectorTemplate::normalizeAndZeroTemplateValues() {
  const double pedestal = calcSimplePedestal(m_templateValues);
  const double peakVal = m_templateValues[m_peakIdx];
  for (size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
    m_templateValues[iTemplateIdx] = (m_templateValues[iTemplateIdx] - pedestal) / (peakVal - pedestal);
}

} // namespace myFuncs
