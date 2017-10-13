#include "filterFuncs.h"

// STL
#include <cmath>

//Mine
#include "mathFuncs.h"

namespace myFuncs {
namespace DSP {
// Implementations taken from "Recursive algorithms for real time digital CR-(RC)^n pulse shaping", M. Nakhostin, IEEE
// transactions on nuclear science, Vol 58, no 5, october 2011

//----------------------------------------------------------
// filterCR_RC
//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC(const std::vector<Type> &xs, const double tau, const double T) {

  //---------------------------
  // xs.size() == 0
  //---------------------------
  if (xs.size() == 0) {
    return std::vector<double>();
  }

  //---------------------------
  // xs.size() >= 1
  //---------------------------
  const double normFactor = 1.0; //T / tau; // Peak amplitude = tau^4 / T.
  // return vector
  std::vector<double> ys;
  ys.reserve(xs.size());
  ys.push_back(xs[0]); // ys[0]
  if (xs.size() == 1) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 2
  //---------------------------

  // Define constants

  const double a = 1. / tau;
  const double alpha = std::exp(-T * a);
  const double xPreviousMultiplier = -alpha * (1.0 + a * T);
  const double yPreviousMultiplier = 2 * alpha;

  ys.push_back(yPreviousMultiplier * ys[0] + xs[1] + xPreviousMultiplier * xs[0]); // ys[1]
  if (xs.size() == 2) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 3
  //---------------------------
  const double yPrevious2Multiplier = -alpha * alpha;

  for (unsigned int idx = 2; idx < xs.size(); ++idx) {
    ys.push_back(yPreviousMultiplier * ys[idx - 1] + yPrevious2Multiplier * ys[idx - 2] + xs[idx] +
                 xPreviousMultiplier * xs[idx - 1]); // ys[>=2]
  }

  return myFuncs::scaleVector(ys, normFactor);
}

// Declare all usefule implementations in order to avoid linker problems!
template std::vector<double> filterCR_RC<double>(const std::vector<double> &xs, const double tau, const double T);
template std::vector<double> filterCR_RC<unsigned int>(const std::vector<unsigned int> &xs, const double tau, const double T);
template std::vector<double> filterCR_RC<int>(const std::vector<int> &xs, const double tau, const double T);

//----------------------------------------------------------
// filterCR_RC2
//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC2(const std::vector<Type> &xs, const double tau, const double T) {
  //---------------------------
  // xs.size() == 0
  //---------------------------
  if (xs.size() == 0) {
    return std::vector<double>();
  }

  //---------------------------
  // xs.size() >= 1
  //---------------------------

  // return vector
  std::vector<double> ys;
  ys.reserve(xs.size());
  ys.push_back(0.0); // y[0]
  if (xs.size() == 1) {
    return ys;
  }
  //---------------------------
  // xs.size() >= 2
  //---------------------------

  // Define constants
  const double normFactor = 1.0; //T / tau / tau; // Peak amplitude = tau^2 / T.
  const double a = 1. / tau;
  const double alpha = std::exp(-T * a);
  const double xPreviousMultiplier = T * alpha * (1.0 - 0.5 * a * T);

  ys.push_back(xPreviousMultiplier * xs[0]); // y[1]
  if (xs.size() == 2) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 3
  //---------------------------
  const double xPrevious2Multiplier = -T * alpha * alpha * (1.0 + 0.5 * a * T);
  const double yPreviousMultiplier = 3 * alpha;
  ys.push_back(yPreviousMultiplier * ys[1] + xPreviousMultiplier * xs[1] + xPrevious2Multiplier * xs[0]); // y[2]
  if (xs.size() == 3) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 4
  //---------------------------
  const double yPrevious2Multiplier = -3 * alpha * alpha;
  ys.push_back(yPreviousMultiplier * ys[2] + yPrevious2Multiplier * ys[1] + xPreviousMultiplier * xs[2] +
               xPrevious2Multiplier * xs[1]); // y[3]
  if (xs.size() == 4) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() > 4
  //---------------------------
  const double yPrevious3Multiplier = alpha * alpha * alpha;

  for (unsigned int idx = 4; idx < xs.size(); ++idx)
    ys.push_back(yPreviousMultiplier * ys[idx - 1] + yPrevious2Multiplier * ys[idx - 2] + yPrevious3Multiplier * ys[idx - 3] +
                  xPreviousMultiplier * xs[idx - 1] + xPrevious2Multiplier * xs[idx - 2]);

  return myFuncs::scaleVector(ys, normFactor);
} // namespace DSP

// Declare all usefule implementations in order to avoid linker problems!
template std::vector<double> filterCR_RC2<double>(const std::vector<double> &xs, const double tau, const double T);
template std::vector<double> filterCR_RC2<unsigned int>(const std::vector<unsigned int> &xs, const double tau, const double T);
template std::vector<double> filterCR_RC2<int>(const std::vector<int> &xs, const double tau, const double T);

//----------------------------------------------------------
// filterCR_RC4
//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC4(const std::vector<Type> &xs, const double tau, const double T) {
  //---------------------------
  // xs.size() == 0
  //---------------------------
  if (xs.size() == 0) {
    return std::vector<double>();
  }

  //---------------------------
  // xs.size() >= 1
  //---------------------------

  // return vector
  std::vector<double> ys;
  ys.reserve(xs.size());
  ys.push_back(0.0); // y[0]
  if (xs.size() == 1) {
    return ys;
  }

  //---------------------------
  // xs.size() >= 2
  //---------------------------
  // Define constants
  const double normFactor = 1.0; //T / tau / tau / tau / tau; // Peak amplitude = tau^4 / T.
  const double a = 1. / tau;
  const double alpha = std::exp(-T * a);
  constexpr double oneO24 = 1.0 / 24.0;
  const double T3 = T * T * T;
  const double T4 = T3 * T;
  const double xPreviousMultiplier = oneO24 * alpha * (-a * T4 + 4.0 * T3);

  ys.push_back(xPreviousMultiplier * xs[0]); // y[1]
  if (xs.size() == 2) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 3
  //---------------------------
  const double xPrevious2Multiplier = oneO24 * alpha * alpha * (-11.0 * a * T4 + 12.0 * T3);
  const double yPreviousMultiplier = 5.0 * alpha;
  ys.push_back(yPreviousMultiplier * ys[1] + xPreviousMultiplier * xs[1] + xPrevious2Multiplier * xs[0]); // y[2]
  if (xs.size() == 3) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 4
  //---------------------------
  const double xPrevious3Multiplier = oneO24 * alpha * alpha * alpha * (-11.0 * a * T4 - 12.0 * T3);
  const double yPrevious2Multiplier = -10.0 * alpha * alpha;
  ys.push_back(yPreviousMultiplier * ys[2] + yPrevious2Multiplier * ys[1] + xPreviousMultiplier * xs[2] +
                xPrevious2Multiplier * xs[1] + xPrevious3Multiplier * xs[0]); // y[3]
  if (xs.size() == 4) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 5
  //---------------------------
  const double xPrevious4Multiplier = oneO24 * alpha * alpha * alpha * alpha * (-a * T4 - 4.0 * T3);
  const double yPrevious3Multiplier = 10.0 * alpha * alpha * alpha;
  ys.push_back(yPreviousMultiplier * ys[3] + yPrevious2Multiplier * ys[2] + yPrevious3Multiplier * ys[1] +
                xPreviousMultiplier * xs[3] + xPrevious2Multiplier * xs[2] + xPrevious3Multiplier * xs[1] +
                xPrevious4Multiplier * xs[0]); // y[4]
  if (xs.size() == 5) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() >= 6
  //---------------------------
  const double yPrevious4Multiplier = -5.0 * alpha * alpha * alpha * alpha;
  ys.push_back(yPreviousMultiplier * ys[4] + yPrevious2Multiplier * ys[3] + yPrevious3Multiplier * ys[2] +
                yPrevious4Multiplier * ys[1] + xPreviousMultiplier * xs[4] + xPrevious2Multiplier * xs[3] +
                xPrevious3Multiplier * xs[2] + xPrevious4Multiplier * xs[1]); // y[5]
  if (xs.size() == 6) {
    return myFuncs::scaleVector(ys, normFactor);
  }

  //---------------------------
  // xs.size() > 6
  //---------------------------
  const double yPrevious5Multiplier = alpha * alpha * alpha * alpha * alpha;

  for (unsigned int idx = 6; idx < xs.size(); ++idx)
    ys.push_back(yPreviousMultiplier * ys[idx - 1] + yPrevious2Multiplier * ys[idx - 2] + yPrevious3Multiplier * ys[idx - 3] +
                  yPrevious4Multiplier * ys[idx - 4] + yPrevious5Multiplier * ys[idx - 5] + xPreviousMultiplier * xs[idx - 1] +
                  xPrevious2Multiplier * xs[idx - 2] + xPrevious3Multiplier * xs[idx - 3] + xPrevious4Multiplier * xs[idx - 4]);

  return myFuncs::scaleVector(ys, normFactor);
} // namespace myFuncs

// Declare all usefule implementations in order to avoid linker problems!
template std::vector<double> filterCR_RC4<double>(const std::vector<double> &xs, const double tau, const double T);
template std::vector<double> filterCR_RC4<unsigned int>(const std::vector<unsigned int> &xs, const double tau, const double T);
template std::vector<double> filterCR_RC4<int>(const std::vector<int> &xs, const double tau, const double T);

} // namespace DSP
} // namespace myFuncs
