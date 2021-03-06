#pragma once

#include <algorithm> //for for_each, reverse
#include <iostream>
#include <iterator>
#include <numeric> //for inner_product
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
std::vector<double> filterCR_RC(const std::vector<Type>& xs, const double tau = 1., const double T = 0.);

//----------------------------------------------------------
// filterCR_RC2

// Do a digital CR-RC^2 filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC2(const std::vector<Type>& xs, const double tau = 1., const double T = 0.);

//----------------------------------------------------------
// filterCR_RC4

// Do a digital CR-RC^4 filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC4(const std::vector<Type>& xs, const double tau = 1., const double T = 0.);

template <typename Type>
std::vector<double> filterBrickwall(const std::vector<Type>& xs, const std::size_t numToBrick, const Type valueToBrick = 0);

//----------------------------------------------------------
// filterByTransferFunction
// Multiply the fft of timeDomainValues by transferFunction and return the transform back to time domain.
//----------------------------------------------------------
template <class timeType, class transferType>
std::vector<double> filterByTransferFunction(const std::vector<timeType>& timeDomainValues,
                                             const std::vector<transferType>& transferFunction) {
  auto fftPair = myFuncs::fftR2C(timeDomainValues);

  std::vector<std::complex<double>> fftComplex;
  for (size_t i = 0; i < fftPair.first.size(); ++i) {
    fftComplex.push_back(std::complex<double>(fftPair.first[i], fftPair.second[i]));
  }

  return myFuncs::fftC2R(timeDomainValues.size(), myFuncs::multiplyElementByElement(fftComplex, transferFunction));
}
template <class timeType, class... Ts>
std::vector<double> filterByTransferFunction(const std::vector<timeType>& timeDomainValues, Ts&... others) {
  return filterByTransferFunction(timeDomainValues, myFuncs::multiplyElementByElement(others...));
}

// Based on Discrete time signal processing,Oppenheim, Schafer, 3rd edition, eq. 3.66 on page 133.
// y[n] = -sum_{k=1}^N(a_k/a_0)y[n-k] + sum_{k=0}^M(b_k/a_0)x[n-k]
// Nominators - b_k, denominators a_k
template <class coefficientType, class inputType>
std::vector<decltype(coefficientType() * inputType())>
filter(std::vector<coefficientType> nominators, std::vector<coefficientType> denominators, const std::vector<inputType>& inputs) {
  using returnType = decltype(coefficientType() * inputType());
  std::vector<returnType> outputs;
  outputs.reserve(inputs.size());
  // Sanity checks
  if (denominators.empty() or denominators.front() == 0) {
    std::cerr << "filterFuncs::filter: Error: denominators must have non zero first element element \n";
    return outputs;
  }

  // Scale nominators by denominators[0]
  const auto ov_denom0 = 1.0 / denominators.front();
  std::transform(nominators.begin(), nominators.end(), nominators.begin(),
                 [&ov_denom0](const auto& nom) { return nom * ov_denom0; });

  // Scale denominators by denominators[0]
  // scaledDenominators has denominators.size() - 1 elements!
  std::transform(std::next(denominators.begin()), denominators.end(), denominators.begin(),
                 [&ov_denom0](const auto& denom) { return denom * ov_denom0; });
  denominators.resize(denominators.size() - 1);

  const auto term1 = [&nominators, &denominators, &outputs]() {
    return outputs.size() >= denominators.size()
               ? std::inner_product(denominators.begin(), denominators.end(), outputs.rbegin(), returnType(0))
               : std::inner_product(outputs.rbegin(), outputs.rend(), denominators.begin(), returnType(0));
  };

  const auto term2 = [&nominators, &denominators, &inputs](const uint n) {
    const size_t numTerms = std::min({static_cast<size_t>(n + 1), nominators.size(), inputs.size()});

    // Had problems when trying to implement with std::inner_product because of either reverse iterator problems, or when some
    // containers are empty
    double term2 = 0.0;
    for (size_t idx = 0; idx < numTerms; ++idx) {
      term2 += nominators[idx] * inputs[n - idx];
    }
    return term2;
  };

  // Take care of y[0]
  if (!inputs.empty()) {
    outputs.push_back(term2(0));
  }

  for (uint n = 1; n < inputs.size(); ++n) {
    // std::cout << "term1(" << n << ") = " << term1() << ", term2(" << n << ") = " << term2(n) << "\n";
    outputs.push_back(-term1() + term2(n));
  }
  return outputs;
}

template <class coefficientType, class inputType>
std::vector<decltype(coefficientType() * inputType())> filter(const std::vector<coefficientType>& originalNominators,
                                                              const coefficientType originalDenominator,
                                                              const std::vector<inputType>& inputs) {
  return myFuncs::DSP::filter(originalNominators, std::vector<coefficientType>{originalDenominator}, inputs);
}

template <class coefficientType, class inputType>
std::vector<decltype(coefficientType() * inputType())>
filter(const std::pair<std::vector<coefficientType>, std::vector<coefficientType>>& coeficientPair,
       const std::vector<inputType>& inputs) {
  return myFuncs::DSP::filter(coeficientPair.first, coeficientPair.second, inputs);
}

std::pair<std::vector<double>, std::vector<double>> getCR_RCnCoefficients(const int n, const double tau,
                                                                          const double samplingFrequency);

// void print(const std::vector<double> &y) {
//   for (const auto &item : y)
//     std::cout << item << " ";
//   std::cout << "\n";
// }
// zero-phase forward and backward filtering.
// Inspiration form MATLASB's filtfilt (R2017a)
template <class coefficientType, class inputType>
std::vector<decltype(coefficientType() * inputType())> filtfilt(const std::vector<coefficientType>& nominators,
                                                                const std::vector<coefficientType>& denominators,
                                                                const std::vector<inputType>& inputs) {

  const int filterOrder = std::max(nominators.size(), denominators.size());
  const uint nfact = std::max(1, 3 * (filterOrder - 1));

  using returnType = decltype(coefficientType() * inputType());
  std::vector<returnType> outputs;

  if (inputs.size() <= nfact) {
    std::cerr << "filterFuncs::filter: Error: originalDenominators must have non zero first element element \n";
    return outputs;
  }

  if (nominators.size() == 0 or denominators.size() == 0) {
    return outputs;
  }

  outputs.reserve(inputs.size());

  // Use outputs vector for intermediate steps as well
  std::vector<returnType> adjustedInputs;
  adjustedInputs.reserve(inputs.size() + nfact * 2);

  //--------------------------------------------------------------------
  // Prepair inputs:
  // Append things at begginng and end of inputs to reduce transients
  //--------------------------------------------------------------------
  const inputType firstInputTimes2 = 2 * inputs.front();
  const inputType lastInputTimes2 = 2 * inputs.back();

  std::for_each(
      std::make_reverse_iterator(inputs.begin() + nfact + 1), inputs.rend() - 1,
      [&firstInputTimes2, &adjustedInputs](const inputType& input) { adjustedInputs.push_back(firstInputTimes2 - input); });
  adjustedInputs.insert(adjustedInputs.end(), inputs.begin(), inputs.end());
  std::for_each(inputs.rbegin() + 1, inputs.rbegin() + nfact + 1, [&lastInputTimes2, &adjustedInputs](const inputType& input) {
    adjustedInputs.push_back(lastInputTimes2 - input);
  });

  // for_each(adjustedInputs.begin(), adjustedInputs.end(), [](auto n) { std::cout << n << " "; });

  //--------------------------------------------------------------------
  // Now filter
  //--------------------------------------------------------------------
  adjustedInputs = myFuncs::DSP::filter(nominators, denominators, adjustedInputs);
  std::reverse(adjustedInputs.begin(), adjustedInputs.end());
  adjustedInputs = myFuncs::DSP::filter(nominators, denominators, adjustedInputs);

  outputs.insert(outputs.end(), adjustedInputs.rbegin() + nfact, adjustedInputs.rend() - nfact);
  return outputs;
}

template <class coefficientType, class inputType>
std::vector<decltype(coefficientType() * inputType())> filtfilt(const std::vector<coefficientType>& nominators,
                                                                const coefficientType denominator,
                                                                const std::vector<inputType>& inputs) {
  return myFuncs::DSP::filtfilt(nominators, std::vector<coefficientType>{denominator}, inputs);
}

} // namespace DSP

} // namespace myFuncs
