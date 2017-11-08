#pragma once

#include <string>
#include <vector>

// ROOT
#include "TComplex.h"

namespace myFuncs {

//--------------------------------------------
// Does the real to complex FFT.
// Normalizes output by 1/sqrt(input.size) - output has units of (input units) ~= sqrt(power).
// The 1/sqrt(input.size)  makes fft and inverse fft more symetric.
//
// Return pair of vectors - first is real part, second is imaginary part.
// Done like this because fftC2R doesn't use TComplex, but needs 2 input vectors - real and imaginary part.

// Insperation taken from CDMS's PulseTools::RealToComplexFFT
// Which in turn was taken from MATLAB's "Power Spectral Density Estimates Using FFT" example (I think)
//--------------------------------------------
std::pair<std::vector<double>, std::vector<double>> fftR2C(const std::vector<double> &input, const std::string &options = "ES");

//--------------------------------------------
// Does the complex to real FFT.
// Normalizes output by 1/sqrt(input.size)

// Requires the number of points in the time domain (they are unrecoverable because of the rounding down of N/2 +1 freq domain
// points)
//--------------------------------------------
std::vector<double> fftC2R(const size_t numTimePoints, const std::vector<double> &realParts, const std::vector<double> &imagParts,
                           const std::string &options = "ES");

//--------------------------------------------
// Takes in a real time domain sequence and produces the psd in units of |input units|^2/Hz.
// Output units are (input units)^2/Hz ~= power/Hz.

// Properly multiplies by 2 all freuqncies except the DC component and (Nyquist component for even lengthed inputs).

// 	The kth output corresponds to frequency k/input.size() in units of cycles/sample

// Insperation taken from CDMS's PulseTools::RealToComplexFFT
// Which in turn was taken from MATLAB's "Power Spectral Density Estimates Using FFT" example (I think)
//--------------------------------------------
std::vector<double> realSequence2psd(const std::vector<double> &input, const double sampleingFrequency = 1.0);

//--------------------------------------------
// Get the frequencies (x-axis) corresponding to a real fft
// N is the size (number of samples) of the original signal (not the FFT length!!! - because of the rounding down of N/2, we don't
// know exactly how many original samples there were)  The output is in units of [1 / samplingInterval]. I.e. if samplingFrequency
// is in ns, output is in GHz.  The output range is (0,samplingFrequency/2).

//--------------------------------------------
std::vector<double> getRealFftfrequencies(const size_t N, const double samplingFrequency);

//--------------------------------------------
// Pad a vector with zeros at the end, until it's size is newSize.
// Return result as new vector
// If newSize < vector.size(), function won't do anything
//--------------------------------------------
template <typename T>
std::vector<T> padWithZeros(const std::vector<T> &vector, const size_t newSize) {
  std::vector<T> padded(vector);
  padded.reserve(newSize);
  while (padded.size() < padded.capacity())
    padded.push_back(0);

  return padded;
}

} // namespace myFuncs
