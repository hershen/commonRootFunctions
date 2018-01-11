#include "catch.hpp"

#include "VectorTemplate.h"

TEST_CASE("Test empty template", "[empty]") {
  myFuncs::VectorTemplate vectorTemplate;
  SECTION("Test get/setDx") {
    vectorTemplate.setDx(0.0);
    CHECK(vectorTemplate.getDx() == Approx(0));

    CHECK_THROWS_AS(vectorTemplate.setDx(-1.0), std::invalid_argument);

    vectorTemplate.setDx(5);
    CHECK(vectorTemplate.getDx() == Approx(5));

    vectorTemplate.setDx(5.5);
    CHECK(vectorTemplate.getDx() == Approx(5.5));
  }
}
// MVectorTemplate();
//
// template <class T>
// MVectorTemplate(const std::vector<T> &newVec, const double newDx);
//
// inline unsigned int getNumAveragedFuncs() const { return m_numAveragedFuncs; }
// // Should be used with care
// inline void setNumAveragedFuncs(const unsigned int numAveraged) { m_numAveragedFuncs = numAveraged; }
//

//
// // inline void setVectorValues(const std::vector<double> &newValues) {m_templateValues = newValues; m_numAveragedFuncs = 0;}
// inline std::vector<double> getTemplateValues() const { return m_templateValues; }
//
// void resetTemplateValues();
//
// // ------------------------------------------------------------------
// // enableAmplitudeFit
// // Currently no way of actually setting amplitude of fit. It is guessed when a vector is added.
// // ------------------------------------------------------------------
// void enableAmplitudeFit(const bool enableAmplitudeFit);
// inline bool isAmplitudeFitEnabled() const { return m_amplitudeFitEnabled; }
//
// // ------------------------------------------------------------------
// // enablePedestalFit
// // Currently no way of actually setting pedestal of fit. It is guessed when a vector is added.
// // ------------------------------------------------------------------
// void enablePedestalFit(const bool enablePedestalFit);
// inline bool isPedestalFitEnabled() const { return m_pedestalFitEnabled; }
//
// void enableXshiftFit(const bool enableXshiftFit);
// inline bool isXshiftFitEnabled() const { return m_xShiftFitEnabled; }
//
// inline unsigned int getTemplateSize() const { return m_templateValues.size(); }
//
// // Adds a new vector to the template. The values are averaged (taking into account the relative wight of the new vector
// // according to the number of previous vectors added.  m_tF1 parameters reset at the end (otherwise they keep values of last
// // vector added.
// template <class T>
// TFitResult addVector(const std::vector<T> &newVector, const double std, TF1 &function);
//
// inline TF1 *getTF1() { return &m_tF1; }
//
// inline int getDebugLevel() const { return m_debugLevel; }
// inline void setDebugLevel(const int newDebugLevel) { m_debugLevel = newDebugLevel; }
//
// inline void getAmplitudeLimits(double &minAmplitude, double &maxAmplitude) const {
//   minAmplitude = m_minAmplitudeLimit;
//   maxAmplitude = m_maxAmplitudeLimit;
// }
// void setAmplitudeLimits(const double newMinAmplitudeLimit, const double newMaxAmplitudeLimit);
//
// inline void getChi2_NdfLimits(double &minChi2_Ndf, double &maxChi2_Ndf) const {
//   minChi2_Ndf = m_minChi2_NdfLimit;
//   maxChi2_Ndf = m_maxChi2_NdfLimit;
// }
// void setChi2_NdfLimits(const double newMinChi2_NdfLimit, const double newMaxChi2_NdfLimit) {
//   m_minChi2_NdfLimit = newMinChi2_NdfLimit;
//   m_maxChi2_NdfLimit = newMaxChi2_NdfLimit;
// }
//
// inline void getXshiftLimits(double &minXshift, double &maxXshift) const {
//   minXshift = m_minXshiftLimit;
//   maxXshift = m_maxXshiftLimit;
// }
// void setXshiftLimits(const double newMinXshiftLimit, const double newMaxXshiftLimit);
//
// inline void getPedestalLimits(double &minPedestal, double &maxPedestal) const {
//   minPedestal = m_minPedestalLimit;
//   maxPedestal = m_maxPedestalLimit;
// }
// void setPedestalLimits(const double newMinPedestalLimit, const double newMaxPedestalLimit);
//
// // Shouldn't really be used.
// inline void setPeakIdx(const size_t peakIdx) { m_peakIdx = peakIdx; }
//
// inline size_t getPeakIdx() const { return m_peakIdx; }
//
// int saveTemplateToTFile(const std::string &fullFileName, const std::string &treeDescription);
//
// int loadTemplateFromTFile(const std::string &fullFileName);
//
// void setTF1Parameters(TF1 &function, const double amplitude, const double pedestal, const double xShift);
//
// // This is used to keep track where the x axis 0 is.
// // It can change if items are removed from the beggining of the template, for example.
// inline double getXvalueOfFirstTemplateEntry() const { return m_xValueOfFirstTemplateEntry; }
// void setXvalueOfFirstTemplateEntry(const double xValueOfFirstTemplateEntry);
//
// TFitResult fitFunctionToVector(const std::vector<double> &vector, const double std, TF1 &function);
//
// AverageMode getAverageMode() const { return m_averageMode; }
// void setAverageMode(const AverageMode &averageMode) { m_averageMode = averageMode; }
//
// // Normalize peak m_templateValues to 1. And make sure pedestal is 0.
// void normalizeAndZeroTemplateValues();
