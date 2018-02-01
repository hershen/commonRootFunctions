#include "filterFuncs.h"

// STL
#include <cmath>

// Mine
#include "fftFuncs.h"
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
  const double normFactor = T / tau; // Peak amplitude = tau / T.
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
  const double normFactor = T / tau / tau; // Peak amplitude = tau^2 / T.
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
  const double normFactor = T / tau / tau / tau / tau; // Peak amplitude = tau^4 / T.
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

template <typename Type>
std::vector<double> filterBrickwall(const std::vector<Type> &xs, const std::size_t numToBrick, const Type valueToBrick) {
  auto fft = myFuncs::fftR2C(xs);
  for (size_t i = 0; i < numToBrick; ++i) {
    fft.first[i] = valueToBrick;
    fft.second[i] = valueToBrick;
  }

  return myFuncs::fftC2R(xs.size(), fft.first, fft.second);
}
template std::vector<double> filterBrickwall<double>(const std::vector<double> &xs, const std::size_t numToBrick,
                                                     const double valueToBrick);
// template std::vector<double> filterBrickwall<unsigned int>(const std::vector<unsigned int> &xs, const std::size_t numToBrick,
// const unsigned int valueToBrick); template std::vector<double> filterBrickwall<int>(const std::vector<int> &xs, const
// std::size_t numToBrick, const int valueToBrick);

std::pair<std::vector<double>, std::vector<double>> getCR_RCnCoefficients(const int n, const double tau,
                                                                          const double samplingFrequency) {

  const double a = 1.0 / tau;
  const double T = 1.0 / samplingFrequency;
  const double alpha = std::exp(-T / tau);

  const double norm = T; // I don't know why this is needed. It doesn't come out of the z transform of the time response
  std::vector<double> nominators;
  std::vector<double> denominators;
  if (n == 1) {
    nominators = {a * norm, -a * alpha * (1 + a * T) * norm};
    denominators = {1.0, -2 * alpha, alpha * alpha};
  } else if (n == 2) {
    nominators = {0, a * a * T * alpha * (2 - a * T) * norm, -a * a * T * alpha * alpha * (2 + a * T) * norm};
    denominators = {2, -6 * alpha, 6 * alpha * alpha, -2 * alpha * alpha * alpha};
  } else if (n == 3) {
    nominators = {0, a * a * a * T * T * alpha * (3 - a * T) * norm, a * a * a * T * T * alpha * alpha * (-4 * a * T) * norm,
                  -a * a * a * T * T * alpha * alpha * alpha * (3 + a * T) * norm};
    denominators = {6, -24 * alpha, 36 * alpha * alpha, -24 * alpha * alpha * alpha, 6 * alpha * alpha * alpha * alpha};
  } else if (n == 4) {
    nominators = {0, a * a * a * a * alpha * T * T * T * (4 - a * T) * norm,
                  a * a * a * a * alpha * alpha * T * T * T * (12 - 11 * a * T) * norm,
                  a * a * a * a * alpha * alpha * alpha * T * T * T * (-12 - 11 * a * T) * norm,
                  a * a * a * a * alpha * alpha * alpha * alpha * T * T * T * (-4 - a * T) * norm};
    denominators = {24,
                    -120 * alpha,
                    240 * alpha * alpha,
                    -240 * alpha * alpha * alpha,
                    120 * alpha * alpha * alpha * alpha,
                    -24 * alpha * alpha * alpha * alpha * alpha};
  } else {
    std::cout << "filterFuncs::getCR_RCnCoefficients: n = " << n << " not supported. \n";
  }

  return std::make_pair(nominators, denominators);
}

} // namespace DSP
} // namespace myFuncs
