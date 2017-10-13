#include "filterFuncs.h"

#include "TMath.h"
#include <iostream>

namespace myFuncs {
namespace DSP {

template <typename Type>
std::vector<double> filterCR_RC2(std::vector<Type> xV, double tau, double T) {
  // return vector
  std::vector<double> yV = std::vector<double>(xV.size(), 0);

  // ------------------------------------------------------------------
  // Sanity checks
  // ------------------------------------------------------------------
  if (xV.size() < 1) {
    std::cout << "filterCR_RC2: empty input vector xV. Aborting." << std::endl;
    return yV;
  }
  double effT = T;
  if (std::abs(effT) < 1e-9)
    effT = 2. / 500. * tau;

  double a = 1. / tau;
  double alpha = TMath::Exp(-effT / tau);

  double xPrevious = xV[0], xPrevious2 = xV[0];
  double yPrevious = 0., yPrevious2 = 0., yPrevious3 = 0.;

  for (uint idx = 0; idx < xV.size(); ++idx) {
    if (idx == 0) {
    } // All previous are zeros
    else if (idx == 1) {
      xPrevious = xV[idx - 1];
      yPrevious = yV[idx - 1];
    } else if (idx == 2) {
      xPrevious = xV[idx - 1];
      xPrevious2 = xV[idx - 2];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
    } else // idx > 2
    {
      xPrevious = xV[idx - 1];
      xPrevious2 = xV[idx - 2];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
      yPrevious3 = yV[idx - 3];
    }

    yV[idx] = 3. * alpha * yPrevious - 3. * pow(alpha, 2) * yPrevious2 + pow(alpha, 3) * yPrevious3 +
              effT * alpha * (1. - 0.5 * a * effT) * xPrevious - effT * std::pow(alpha, 2) * (1. + 0.5 * a * effT) * xPrevious2;
  }

  // Normalize
  for (uint idx = 0; idx < yV.size(); ++idx) {
    yV[idx] = yV[idx] * effT / pow(tau, 2);
  }
  return yV;
}

// Declare all usefule implementations in order to avoid linker problems!
template std::vector<double> filterCR_RC2<double>(std::vector<double> xV, double tau, double T);
template std::vector<double> filterCR_RC2<unsigned int>(std::vector<unsigned int> xV, double tau, double T);
template std::vector<double> filterCR_RC2<int>(std::vector<int> xV, double tau, double T);

template <typename Type>
std::vector<double> filterCR_RC4(std::vector<Type> xV, double tau, double T) {
  // return vector
  std::vector<double> yV = std::vector<double>(xV.size(), 0);

  // ------------------------------------------------------------------
  // Sanity checks
  // ------------------------------------------------------------------
  if (xV.size() < 1) {
    std::cout << "filterCR_RC4: empty input vector xV. Aborting." << std::endl;
    return yV;
  }
  if (std::abs(T) < 1e-9)
    T = 2. / 500. * tau;
  double a = 1. / tau;
  double alpha = TMath::Exp(-T / tau);

  double xPrevious = xV[0], xPrevious2 = xV[0], xPrevious3 = xV[0], xPrevious4 = xV[0];
  double yPrevious = 0., yPrevious2 = 0., yPrevious3 = 0., yPrevious4 = 0., yPrevious5 = 0.;

  for (uint idx = 0; idx < xV.size(); ++idx) {
    if (idx == 0) {
    } // All previous are zeros
    else if (idx == 1) {
      xPrevious = xV[idx - 1];
      yPrevious = yV[idx - 1];
    } else if (idx == 2) {
      xPrevious = xV[idx - 1];
      xPrevious2 = xV[idx - 2];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
    } else if (idx == 3) {
      xPrevious = xV[idx - 1];
      xPrevious2 = xV[idx - 2];
      xPrevious3 = xV[idx - 3];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
      yPrevious3 = yV[idx - 3];
    } else if (idx == 4) {
      xPrevious = xV[idx - 1];
      xPrevious2 = xV[idx - 2];
      xPrevious3 = xV[idx - 3];
      xPrevious4 = xV[idx - 4];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
      yPrevious3 = yV[idx - 3];
      yPrevious4 = yV[idx - 4];
    } else // idx > 4
    {
      xPrevious = xV[idx - 1];
      xPrevious2 = xV[idx - 2];
      xPrevious3 = xV[idx - 3];
      xPrevious4 = xV[idx - 4];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
      yPrevious3 = yV[idx - 3];
      yPrevious4 = yV[idx - 4];
      yPrevious5 = yV[idx - 5];
    }

    yV[idx] = 5. * alpha * yPrevious - 10. * pow(alpha, 2) * yPrevious2 + 10. * pow(alpha, 3) * yPrevious3 -
              5. * pow(alpha, 4) * yPrevious4 + pow(alpha, 5) * yPrevious5 +
              1. / 24. * (-a * pow(T, 4) * alpha + 4. * pow(T, 3) * alpha) * xPrevious +
              1. / 24. * (-11. * a * pow(T, 4) * pow(alpha, 2) + 12. * pow(T, 3) * pow(alpha, 2)) * xPrevious2 +
              1. / 24. * (-12. * pow(T, 3) * pow(alpha, 3) - 11. * a * pow(T, 4) * pow(alpha, 3)) * xPrevious3 +
              1. / 24. * (-a * pow(T, 4) * pow(alpha, 4) - 4. * pow(T, 3) * pow(alpha, 4)) * xPrevious4;
  }
  for (uint idx = 0; idx < yV.size(); ++idx) {
    yV[idx] = yV[idx] * T / pow(tau, 4);
  }
  return yV;
}

// Declare all usefule implementations in order to avoid linker problems!
template std::vector<double> filterCR_RC4<double>(std::vector<double> xV, double tau, double T);
template std::vector<double> filterCR_RC4<unsigned int>(std::vector<unsigned int> xV, double tau, double T);
template std::vector<double> filterCR_RC4<int>(std::vector<int> xV, double tau, double T);

std::vector<double> filterCR_RC(std::vector<double> xV, double tau, double T) {
  // return vector
  std::vector<double> yV = std::vector<double>(xV.size(), 0);

  // ------------------------------------------------------------------
  // Sanity checks
  // ------------------------------------------------------------------
  if (xV.size() < 1) {
    std::cout << "filterCR_RC: empty input vector xV. Aborting." << std::endl;
    return yV;
  }
  if (std::abs(T) < 1e-9)
    T = 2. / 500. * tau;
  double a = 1. / tau;
  double alpha = TMath::Exp(-T / tau);

  double xPrevious = 0., yPrevious = 0., yPrevious2 = 0.;
  for (uint idx = 0; idx < xV.size(); ++idx) {
    if (idx == 0)
      xPrevious = xV[0];
    else if (idx == 1) {
      xPrevious = xV[idx - 1];
      yPrevious = yV[idx - 1];
    } else // idx > 1
    {
      xPrevious = xV[idx - 1];
      yPrevious = yV[idx - 1];
      yPrevious2 = yV[idx - 2];
    }

    yV[idx] = 2. * alpha * yPrevious - pow(alpha, 2) * yPrevious2 + xV[idx] - alpha * (1. + a * T) * xPrevious;
  }

  // Normalize
  for (size_t idx = 0; idx < yV.size(); ++idx)
    yV[idx] = yV[idx] * T / tau;

  return yV;
}

} // namespace DSP
} // namespace myFuncs
