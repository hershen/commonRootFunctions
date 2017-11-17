#include <vector>

// Mine
#include "fftFuncs.h"
#include "generalFuncs.h"

namespace myFuncs {
namespace DSP {

//----------------------------------------------------------
// filterCR_RC

// Do a digital CR-RC filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC(const std::vector<Type> &xs, const double tau = 1., const double T = 0.);

//----------------------------------------------------------
// filterCR_RC2

// Do a digital CR-RC^2 filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC2(const std::vector<Type> &xs, const double tau = 1., const double T = 0.);

//----------------------------------------------------------
// filterCR_RC4

// Do a digital CR-RC^4 filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC4(const std::vector<Type> &xs, const double tau = 1., const double T = 0.);

template <typename Type>
std::vector<double> filterBrickwall(const std::vector<Type> &xs, const std::size_t numToBrick, const Type valueToBrick = 0);

//----------------------------------------------------------
// filterByTransferFunction
// Multiply the fft of timeDomainValues by transferFunction and return the transform back to time domain.
//----------------------------------------------------------
template <class timeType, class transferType>
std::vector<double> filterByTransferFunction(const std::vector<timeType> &timeDomainValues,
                                             const std::vector<transferType> &transferFunction) {
  auto fftPair = myFuncs::fftR2C(timeDomainValues);

  std::vector<std::complex<double>> fftComplex;
  for (size_t i = 0; i < fftPair.first.size(); ++i) {
    fftComplex.push_back(std::complex<double>(fftPair.first[i], fftPair.second[i]));
  }

  return myFuncs::fftC2R(timeDomainValues.size(), myFuncs::multiplyElementByElement(fftComplex, transferFunction));
}
template <class timeType, class... Ts>
std::vector<double> filterByTransferFunction(const std::vector<timeType> &timeDomainValues, Ts &... others) {
  return filterByTransferFunction(timeDomainValues, myFuncs::multiplyElementByElement(others...));
}

} // namespace DSP

} // namespace myFuncs
