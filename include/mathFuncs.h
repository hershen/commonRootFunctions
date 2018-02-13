#pragma once

// std
#include <iostream>
#include <iterator>
#include <numeric> //For std::accumulate

// Root
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TRandom3.h"

// For TMiuit2
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"

class TGraphErrors;

// Todo
// Think what happens if sampleMean, sampleStd have unsigned numbers (what happens with minuses?)

namespace myFuncs {

// trying to compile with g++
//   int vecSize = 0;
//   TF1 globalFunc;

//--------------------------------------------------------------------------------------------
// analytical_RC_CRn
//********************************************************************************************
// Analytical response of an RC-(CR)^n filter to a step function starting at time t = 0.
// Function is 1/n! * (t/tau)^n * exp(-t/tau)
// tau is the shaping time.
// Function is shifted in time so that it starts at time startTime.
// amplitude controls the amplitude.
// max time is an approximate time for the function to decrease below exp(-expCoefficient)
//--------------------------------------------------------------------------------------------
TF1 analytical_RC_CRn(int n, double tau = 500., double amplitude = 1., double startTime = 0., double maxTime = 0.);

//--------------------------------------------------------------------------------------------
// calcResiduals
//********************************************************************************************
// Function to calculate residuals between a a model function (TF1) and y(x) yValues taken at points xValuesV.
// Residual defined as f(x) - yValue
// x values are in vector xValues of type xValueType
// y values are in vector yValues of type yValueType
// If there are more yValues than xValues, the residuals are calculated only for the xValues given.
// More xValues than yValues are not allowed.
// Does not check for zero values in stds (can devide by 0)
//--------------------------------------------------------------------------------------------
// Not implemented with boost zip itirators because ROOT doesn't compile Boost libraries well.
template <typename xValType, typename yValType>
std::vector<double> calcResiduals(const std::vector<xValType>& xValues, const std::vector<yValType>& yValues,
                                  const std::vector<double>& stds, const TF1& modelFunc) {
  // ------------------------------------------------------------------
  // Sanity checks
  // ------------------------------------------------------------------
  // Allow cases where there are more y values than x values
  if (yValues.size() < xValues.size())
    throw std::invalid_argument("yValuesV.size() < xValuesV.size()");

  if (stds.size() < xValues.size())
    throw std::invalid_argument("stds.size() < xValuesV.size()");

  // Output vector
  std::vector<double> residuals;
  residuals.reserve(xValues.size());

  auto xIt = xValues.begin();
  auto yIt = yValues.begin();
  auto stdIt = stds.begin();
  // Loop on num elements = xValues.size()
  while (xIt != xValues.end()) {
    // 			std::cout << "xVal = " << *xIt << ", yVal = " << *yIt << std::endl;
    residuals.push_back((modelFunc.Eval(*xIt) - static_cast<double>(*yIt)) / *stdIt);
    ++xIt;
    ++yIt;
    ++stdIt;
  }

  return residuals;
}

// Overloaded
template <typename xValType, typename yValType>
std::vector<double> calcResiduals(const std::vector<xValType>& xValues, const std::vector<yValType>& yValues, const double std,
                                  const TF1& modelFunc) {
  return calcResiduals(xValues, yValues, std::vector<double>(xValues.size(), std), modelFunc);
}

// Not tested
std::vector<double> calcResiduals(const TGraphErrors& graphErrors, const TF1& modelFunc);

//--------------------------------------------------------------------------------------------
// convertArray2TF1Internal
//********************************************************************************************
// Internal function for convertVector2TF1.
double convertArray2TF1Internal(double* var, double* params);

//--------------------------------------------------------------------------------------------
// convertVector2TF1
//********************************************************************************************
// Converts a vector of double values into a TF1.
// Only parameter 0 should regularly change! Others change the function itself.

// It is assumed that the values in vecValues are evenly spaced with spacing dT. I.e. vecValues[0] = f(0), vecValues[1] = f(dT),
// vecValues[2] = f(2dT)...  timeShift can shift the x axis.  Values between the ones in vecValues are evaluated using a linear
// interpolation between the points.  Evals for x axis values smaller than the first entry in vecValues return vecValues[0]. Evals
// for values greater than the last entry in vecValues return vecValues[last].  The params array given to convertArray2TF1Internal
// is made up of {timeShift, dT, vecValues}  All parameters are fixed, except timeshift which floats.  There is a file scope
// "global" variable vecSize which holds the number of parameters in vecValues (because the internal function doesn't know this.
// The  All values of vecValues are given to the internal function
//--------------------------------------------------------------------------------------------
TF1 convertVector2TF1(double dT, std::vector<double> vecValues, double timeShift);

double CFDfuncInternal(double* var, double* params);

//--------------------------------------------------------------------------------------------
// CFD
//********************************************************************************************
// Implementing a CFD for TF1
//--------------------------------------------------------------------------------------------
TF1 CFD(TF1 func, double DLY, double fraction);

// Assumes that it is a positive signal
double CFtime(TF1 func, double DLY, double fraction /*, double initialGuess = 0.*/);

//--------------------------------------------------------------------------------------------
// cfdCR_RCnZeroCrossingTime
//********************************************************************************************
// Returns the zero crossing time of a CFD function from an analytical CR-(RC)^n output of a step function starting at t=0

double cfdCR_RCnZeroCrossingTime(int n, double fraction, double tau, double DLY);

//--------------------------------------------------------------------------------------------
// addGaussianNoise
//********************************************************************************************
// Add gaussian noise to an input vector of type Type. Output vector is alway double.
// If seed is not given the default constructor is used. NOTE - the number sequence generated is always the same in this case.
// If a seed is provided, it is used.
//--------------------------------------------------------------------------------------------
template <typename Type>
std::vector<double> addGaussianNoise(std::vector<Type> inputValues, double mean, double sigma, int seed = -1);

//--------------------------------------------------------------------------------------------
// getGaussianFit
//********************************************************************************************
// Fit a gaussian to hist.
// If xMinInitial or xMaxInitial are given, the fit is performed only in the range (xMinInitial, xMaxInitial)
// If xMinFracOfSigma or xMaxFracOfSigma are given, the fit is re-performed in the range (mean - sigma*xMinFracOfSigma, mean +
// sigma*xMaxFracOfSigma) where mean and sigma are taken from the first fit. I.e., the fit can be run again in the range (mean -
// 2std, mean + 2std).  The function returns the fit result in fitResult
//--------------------------------------------------------------------------------------------
TF1 getGaussianFit(TH1D hist, TFitResultPtr& fitResult, double xMinInitial, double xMaxInitial, double xMinFracOfSigma,
                   double xMaxFracOfSigma);

//--------------------------------------------------------------------------------------------
// getGaussianFitResult - Overloaded
//********************************************************************************************
// Same as before, only without the fitResult parameter.
//--------------------------------------------------------------------------------------------
inline TF1 getGaussianFit(TH1D hist, double xMinInitial, double xMaxInitial, double xMinFracOfSigma, double xMaxFracOfSigma);

//--------------------------------------------------------------------------------------------
// getCorrelationMatrix
//********************************************************************************************
// Calculate the correlation matrix from the covariance matrix.
// I.e. M_ij = M_ij / sqrt(M_ii) / sqrt(M_jj)
// If one of the diagonal elements M_ii is zero, return original matrix
//--------------------------------------------------------------------------------------------
TMatrixD getCorrelationMatrix(const TMatrixD& covarianceMatrix);

//--------------------------------------------------------------------------------------------
// printMatrix
//********************************************************************************************
// Print a matrix on screen
// headings - the heading of each row / column
// width - each element will be padded to fill width spaces
// aligh - align to left, right, etc.
// precision - the precision syntax to be used (as in printf )
//--------------------------------------------------------------------------------------------
void printMatrix(const TMatrixD matrix, const std::vector<std::string>& headings, int width = 10,
                 std::ios_base& align(std::ios_base& str) = std::left, const std::string& precision = ".3g");

//--------------------------------------------------------------------------------------------
// drawMatrix
//********************************************************************************************
// Print a matrix on screen
// headings - the heading of each row / column
// precision - the precision syntax to be used (as in printf )
// zMin, zMax - used to set z axis limits
//--------------------------------------------------------------------------------------------
TCanvas* drawMatrix(const TMatrixD matrix, std::string title, const std::vector<std::string>& xAxisHeadings,
                    const std::vector<std::string>& yAxisHeadings, const double zMin = 0.0, const double zMax = 0.0,
                    const std::string& precision = ".3g");

//--------------------------------------------------------------------------------------------
// sumVector
//********************************************************************************************
// Returns sum(x[i])
// wrapper around std::accumulate with inital value of 0.0 so not to confuse it with 0 in which case everything is converted to an
// int.
//--------------------------------------------------------------------------------------------
template <typename T>
T sumVector(const std::vector<T>& vector, const T& initialValue = 0.0) {
  return std::accumulate(vector.begin(), vector.end(), initialValue);
}

//--------------------------------------------------------------------------------------------
// sumVectorSquared
//********************************************************************************************
// Returns sum(x[i]*x[i])
//--------------------------------------------------------------------------------------------
template <typename T>
T sumVectorSquared(const std::vector<T>& vector, const T& initialValue = 0.0) {
  auto sum = initialValue;

  for (auto& element : vector)
    sum += element * element;

  return sum;
}

template <class T1, class T2, class outputT>
std::vector<outputT> sumVectors(const std::vector<T1>& v1, const std::vector<T2>& v2) {

  const size_t outputSize = std::min(v1.size(), v2.size());

  std::vector<outputT> output;
  output.reserve(outputSize);

  for (size_t i = 0; i < outputSize; ++i)
    output.push_back(v1[i] + v2[i]);

  return output;
}

template <class T>
std::vector<T> sumVectors(const std::vector<T>& v1, const std::vector<T>& v2) {
  return sumVectors<T, T, T>(v1, v2);
}

template <class T>
std::vector<T> scaleVector(const std::vector<T>& input, const T factor) {

  // output vector
  std::vector<T> output;
  output.reserve(input.size());

  // Scale vector
  for (auto it = input.begin(); it != input.end(); ++it)
    output.push_back(*it * factor);

  return output;
}

// Maybe can use generate algorithm
template <class T>
std::vector<T> addToVector(const std::vector<T>& vector, const T val) {
  // Create output vector
  std::vector<T> output;
  output.reserve(vector.size());

  std::transform(vector.begin(), vector.end(), std::back_inserter(output), [&](const T element) { return element + val; });
  return output;
}

// Linearly interploate the y value at x given (x0,y0), (x1,y1)
template <class Tx, class Ty>
double linearInterpolate(const Tx x0, const Tx x1, const Ty y0, const Ty y1, const double x) {
  if (x1 - x0 == 0.0) {
    throw std::invalid_argument("mathFuncs::linearInterpolate: zero denominator");
  }
  return static_cast<double>(y1 - y0) / static_cast<double>(x1 - x0) * (x - x0) + y0;
}

// Novosibirsk function
// See H. Ikeda et al. / Nuclear Instruments and Methods in Physics Research A 441 (2000) 401-426
double novosibirsk(const double x, const double norm, const double peak, const double width, const double eta);

// Returns a novosibirsk(see above) TF1 with range between minValue and maxValue
TF1 getNovosibirskTF1(const double minValue, const double maxValue);

double getNovosibirskAmplitude(const double normalization, const double eta);
inline double getNovosibirskAmplitude(const TF1& novo) {
  return getNovosibirskAmplitude(novo.GetParameter("Normalization"), novo.GetParameter("#eta"));
}

// Solve parabola f = p0 + p1*x + p2*x^2 = 0
// pair = (-b+sqrt(delta) )/2a, (-b-sqrt(delta) )/2a
std::pair<double, double> solveParabola(const double p0, const double p1, const double p2);

// f = x/y
// Return sigma_f = sqrt(errX^2/y^2 + x^2/y^4*errY^2 + )
inline double xDivYerror(const double x, const double errX, const double y, const double errY, const double covXY = 0.0) {
  return std::sqrt(errX * errX / y / y + x * x / y / y / y / y * errY * errY - 2.0 * x * covXY / y / y / y);
}

// return a number that has at max 2 significant digits.
// If the 2 significant digits are xy, then if xy > 35, return a number with only x as significant digit
// If xy <= 35, return a number with xy as significant digits.
// y is always rounded using the smaller digits.
// If error < 0 returns error
double round_35rule(double error);

// Round x and keep the first digitsToKeep digits.
// If x == 0.0 or digitsToKeep < 0, return x
// I.e. if x = 5.67890
// roundKeepDigits(5.67890, 0) = 5.67890
// roundKeepDigits(5.67890, 1) = 6
// roundKeepDigits(5.67890, 2) = 5.7
// roundKeepDigits(5.67890, 3) = 5.68
// roundKeepDigits(5.67890, 4) = 5.679
//...
double roundKeepDigits(const double x, const int digitsToKeep);

// If x = w*10^n where 0 < |w| < 10, returns n
inline int exponent10(const double x) { return std::floor(std::log10(std::abs(x))); }

// Rounds x to the same precision as error.
// I.e. number of significant digits in error determines number of significant digits in returned value.
// Assumes error has a maximum of 2 significant digits
double roundAccordingToError(const double x, const double error);

// Return sample mean of samples in values
template <class InputIt>
double sampleMean(InputIt first, InputIt last) {
  const long numElements = std::distance(first, last);
  if (numElements == 0) {
    return 0.0;
  }
  // std::cout << "sum " << std::accumulate(first, last, 0.0) <<
  // " numElements " << numElements <<
  // " returning " << std::accumulate(first, last, 0) / static_cast<double>(numElements) << "\n";
  return std::accumulate(first, last, 0.0) / static_cast<double>(numElements);
}

template <class T>
double sampleMean(const std::vector<T>& values) {
  return sampleMean(values.begin(), values.end());
}

// Sample standard deviation
// Devides by number of elements (not #elements - 1)
template <class InputIt>
double sampleStd(InputIt first, InputIt last) {
  double numElements = 0;
  const double sample_mean = sampleMean(first, last);

  double sum = 0;
  for_each(first, last, [&sum, &numElements, sample_mean](const auto element) {
    sum += (element - sample_mean) * (element - sample_mean);
    ++numElements;
  });
  return numElements > 1 ? std::sqrt(sum / static_cast<double>(numElements)) : 0.0;
}

template <class T>
double sampleStd(const std::vector<T>& values) {
  return sampleStd(values.begin(), values.end());
}

// return pair of weighted sample mean and it's standard deviation.
// If size of two input vectors is different, throw.
// Taken from https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
template <class T>
std::pair<double, double> weightedAverageAndStd(const std::vector<T>& values, const std::vector<double>& stds) {
  if (values.size() != stds.size()) {
    std::cout << "mathFuncs::weightedAverageAndStd : Error : size values (" << values.size() << ") != size stds(" << stds.size()
              << ")" << std::endl;
    throw;
  }

  double nom = 0.0;
  double denom = 0.0;

  for (ulong i = 0; i < values.size(); ++i) {
    const double one_sigma2 = 1.0 / stds[i] / stds[i]; // = 1/sigma_i^2
    nom += values[i] * one_sigma2;                     // = sum(x_i/sigma_i^2)
    denom += one_sigma2;                               // = sum(1/sigma_i^2)
  }

  return std::make_pair(nom / denom, std::sqrt(1.0 / denom));
}

// For unbiased estimator of the population std of a normallly distributed sample (it's different from the unbiased esitmator of
// the variance!!)  https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation

template <class T>
std::vector<double> vectorToDouble(const std::vector<T>& input) {
  std::vector<double> output;
  output.reserve(input.size());

  for (const auto sample : input) {
    output.push_back(static_cast<double>(sample));
  }

  return output;
}

// Shift vector to left or right, interpolating between points.
// If dBins = +1, the shift will be by 1 sample distance to the right.
// If dBins = -1, the shift will be by 1 sample distance to the left.
// This is similar to : dBins = +2, f(x) = x^2 => shiftVector(inputs, +2) will return f(x-2)
// If dt isn't provided, shift is in units of bins. If dt is provided it defines the time between bins and shift is in units of
// time.
template <class Tvector>
std::vector<double> shiftVector(const std::vector<Tvector> inputs, double shift, const double dt = 0.0) {

  std::vector<double> outputs;
  outputs.reserve(inputs.size());

  // Sanity checks
  if (inputs.empty()) {
    return outputs;
  }

  // Convert shift to units of bins
  if (dt != 0.0) {
    shift = shift / dt;
  }

  // Pad beginning
  const int pads = shift <= inputs.size() ? std::ceil(shift) : inputs.size();
  for (int i = 0; i < pads; ++i) {
    outputs.push_back(inputs[0]);
  }

  const double fractionalShift = shift - std::floor(shift);
  const int lastIdx = shift < 0 ? inputs.size() - 1 : inputs.size() - static_cast<int>(std::abs(std::ceil(shift)));
  // std::cout << "fractionalShift = " << fractionalShift << ", shift = " << shift << ", lastIdx = " << lastIdx << "\n";
  for (int idx = shift >= 0 ? 0 : std::abs(static_cast<int>(shift)); idx < lastIdx; ++idx) {
    const double x = idx + fractionalShift;
    const Tvector y1 = inputs[idx];
    const Tvector y2 = inputs[idx + 1];

    outputs.push_back(linearInterpolate(idx, idx + 1, y1, y2, x));
    // std::cout << "idx = " << idx << ", linearInterpolate(x1, x2, y1, y2, x) = " << linearInterpolate(idx, idx + 1, y1, y2, x)
    // << "\n";
  }

  // Pad end
  for (size_t i = outputs.size(); i < inputs.size(); ++i) {
    outputs.push_back(inputs.back());
  }

  return outputs;
}

//----------------------------------------------
// parabola = p0 + p1*x + p2*x^2
//----------------------------------------------
template <class T>
inline double parabola_xMax(const T p1, const T p2) {
  if (p2 != 0) {
    return -p1 / 2.0 / p2;
  }
  return 0.0;
}

template <class T>
inline double parabola_maxValue(const T p0, const T p1, const T p2) {
  const auto xMax = parabola_xMax(p1, p2);
  return p0 + p1 * xMax + p2 * xMax * xMax;
}

inline double parabola_xMax(const TFitResultPtr& fitResult) { return parabola_xMax(fitResult->Value(1), fitResult->Value(2)); }

//----------------------------------------------
// parabola using root TFitResult
//----------------------------------------------
inline double parabola_maxValue(const TFitResultPtr& fitResult, const double xMax) {
  return fitResult->Value(0) + fitResult->Value(1) * xMax + fitResult->Value(2) * xMax * xMax;
}

inline double parabola_maxValue(const TFitResultPtr& fitResult) { return parabola_maxValue(fitResult, parabola_xMax(fitResult)); }

//----------------------------------------------
// linear = p0 + p1*x
//----------------------------------------------
// Return when linear line crosses value
template <class Tfunction, class Tvalue>
inline double linear_crossValue(const Tfunction p0, const Tfunction p1, const Tvalue value = 0) {
  if (p1 != 0) {
    return (value - p0) / static_cast<double>(p1);
  }
  return 0.0;
}

// Return when linear line crosses value
// Assumes fitResult is valid
template <class Tvalue>
inline double linear_crossValue(const TFitResultPtr& fitResult, const Tvalue value) {
  return linear_crossValue(fitResult->Value(0), fitResult->Value(1), value);
}

//----------------------------------------------
// Element by element sum of 2 vectors
//----------------------------------------------
template <class T1, class T2>
std::vector<decltype(T1() * T2())> sumElementByElement(const std::vector<T1>& inputs1, const std::vector<T2>& inputs2) {
  std::vector<decltype(T1() * T2())> outputs;
  outputs.reserve(inputs1.size());
  std::transform(inputs1.begin(), inputs1.end(), inputs2.begin(), std::back_inserter(outputs),
                 [](const auto element1, const auto element2) { return element1 + element2; });
  return outputs;
}

//----------------------------------------------
// Find maximum element checking every n'th element
// I.e. element first, element first + n, first + 2n, etc.
//----------------------------------------------
template <class ForwardIt>
ForwardIt findMaxEvery_n(ForwardIt first, ForwardIt last, const size_t n) {
  if (first == last or n == 0) {
    return first;
  }

  const size_t numElements = std::distance(first, last);
  ForwardIt largest = first;
  std::advance(first, n);
  size_t currentElement = n;
  for (; currentElement < numElements; std::advance(first, n), currentElement += n) {
    if (*largest < *first) {
      largest = first;
    }
  }
  return largest;
}

template <class T>
typename T::const_iterator findMaxEvery_n(const T& container, const size_t n) {
  return findMaxEvery_n(container.begin(), container.end(), n);
}

//----------------------------------------------
// Find maximum element in 2 passes
// First pass looks at elements first, first + n1, first + 2n1, etc.
// If first returned element N, second pass evaluates all elements [N-n+1, N+n-1].
// If these are outside first or last, then they are used as end points.
//----------------------------------------------
template <class ForwardIt>
ForwardIt findMaxEvery_n_ThenBetween(ForwardIt first, ForwardIt last, const size_t n) {
  if (first == last or n == 0) {
    return first;
  }

  const auto firstPassItr = findMaxEvery_n(first, last, n);
  auto newFirst = std::prev(firstPassItr, n - 1);
  if (std::distance(first, newFirst) < 0) {
    newFirst = first;
  }
  auto newLast = std::next(firstPassItr, n - 1);
  if (std::distance(newLast, last) < 0) {
    newLast = last;
  }
  return findMaxEvery_n(newFirst, newLast, 1);
}
template <class T>
typename T::const_iterator findMaxEvery_n_ThenBetween(const T& container, const size_t n) {
  return findMaxEvery_n_ThenBetween(container.begin(), container.end(), n);
}

template <class T>
std::vector<double> averageEach_n(const std::vector<T>& vector, const size_t n) {
  std::vector<double> output;

  if (n == 0) {
    output = vectorToDouble(vector);
    return output;
  }

  const double one_n = 1.0 / static_cast<double>(n);
  auto it = vector.begin();
  for (size_t idx = 0; idx + n - 1 < vector.size(); idx += n) {
    auto nextIt = std::next(it, n);
    const double average = std::accumulate(it, nextIt, 0.0) * one_n;
    output.push_back(average);
    it = nextIt;
  }
  return output;
}
} // namespace myFuncs
