#include "fftFuncs.h"

// STD
#include <iostream>
#include <memory>

// ROOT
#include "TVirtualFFT.h"

namespace myFuncs {

std::pair<std::vector<double>, std::vector<double>> fftR2C(const std::vector<double> &input, const std::string &options) {

  // Should be const, but can't because TVirtualFFT::FFT needs to accept non-const pointer
  // This has units of cycles per sample
  int inputSize(input.size());

  if (inputSize == 0) {
    std::cerr << "fftFuncs::fftR2C - ERROR! empty input passed into this function" << std::endl;
    return std::pair<std::vector<double>, std::vector<double>>();
  }

  // Do the FFT
  std::unique_ptr<TVirtualFFT> fftr2c(TVirtualFFT::FFT(1, &inputSize, ("R2C " + options).data()));
  fftr2c->SetPoints(input.data());
  fftr2c->Transform();

  // Output temperary vectors

  // vectors have to be right size, because we're going to overwrite their internal arrays.
  // There are inputSize/2 + 1 complex numbers in the output (fraction rounded down).
  std::vector<double> real(inputSize / 2 + 1); // fraction rounded down
  std::vector<double> imag(inputSize / 2 + 1); // fraction rounded down

  // Read transformed points
  fftr2c->GetPointsComplex(real.data(), imag.data());

  // Create output vector
  std::vector<TComplex> output;
  output.reserve(real.size());

  // Take care to normalize by 1/sqrt(inputSize) because root's fft doesn't do this.
  const double ov_Sqrt_inputSize = 1.0 / std::sqrt(inputSize);

  for (unsigned int iPoint = 0; iPoint < real.size(); ++iPoint) {
    real[iPoint] = real[iPoint] * ov_Sqrt_inputSize;
    imag[iPoint] = imag[iPoint] * ov_Sqrt_inputSize;
  }

  return std::make_pair(real, imag);
}

std::vector<double> realSequence2psd(const std::vector<double> &input, const double sampleingFrequency) {

  if (input.size() == 0) {
    std::cerr << "fftFuncs::realFft2psd - ERROR! empty input passed into this function" << std::endl;
    return std::vector<double>();
  }

  // Perform fft
  // Holds 2 vectors - first is real parts, second is imaginary parts.
  const auto fft(fftR2C(input));

  std::vector<double> psd;
  psd.reserve(fft.first.size());

  // We assume input is fft of real values. Therefore, Y_0 = Y_N, and only (N/2)+1 complex numbers were calculated. N is the size
  // of input.  This means that when calculating the PSD, we need to take into account all the frequencies that were "dropped".  We
  // do this by multiplying all their corresponding frequencies by 2.  The DC component and the Nyquist component for N odd are
  // never "dropped", so we don't multiply them by 2.

  const double ov_samplingFrequency = 1.0 / sampleingFrequency;

  // Lambda for |c|^2
  auto magSquare = [](const double real, const double imag) -> double { return real * real + imag * imag; };

  // DC component doesn't need to be multiplied by 2
  psd.push_back(magSquare(fft.first[0], fft.second[0]) * ov_samplingFrequency);

  // Loop on all except first (DC) and last (Nyquist) components.
  // They are taken care of seperately
  for (size_t iEntry = 1; iEntry < fft.first.size() - 1; ++iEntry)
    psd.push_back(magSquare(fft.first[iEntry], fft.second[iEntry]) * 2.0 * ov_samplingFrequency);

  if (input.size() % 2 == 0) // input size is even - don't need to multiply Nyquist frequency by 2
    psd.push_back(magSquare(fft.first.back(), fft.second.back()) * ov_samplingFrequency);
  else // input size is odd - need to multiply Nyquist frequency by 2
    psd.push_back(magSquare(fft.first.back(), fft.second.back()) * 2.0 * ov_samplingFrequency);

  return psd;
}

std::vector<double> getRealFftfrequencies(const size_t N, const double samplingFrequency) {
  std::vector<double> xValues;
  xValues.reserve(N);

  const auto numFftFrequencies(N / 2 + 1); // Rounds the division down

  for (size_t i = 0; i < numFftFrequencies; ++i)
    xValues.push_back(static_cast<double>(i) / N * samplingFrequency); // Normalized frequency of ith number is i/N.

  return xValues;
}

std::vector<double> fftC2R(const size_t numTimePoints, const std::vector<double> &realParts, const std::vector<double> &imagParts,
                           const std::string &options) {

  // Non suitable variable because that's what root forces us to use
  int numPoints = numTimePoints;

  // Prepare FFT function
  std::unique_ptr<TVirtualFFT> fftc2r(TVirtualFFT::FFT(1, &numPoints, ("C2R " + options).data()));

  fftc2r->SetPointsComplex(realParts.data(), imagParts.data());
  fftc2r->Transform();

  // Set output vector size because we're going to read into the interal array in next step
  std::vector<double> output(numTimePoints);

  // Read transformed points
  fftc2r->GetPoints(output.data());

  // Normalize by 1/sqrt(numTimePoints)
  for (size_t iPoint = 0; iPoint < output.size(); ++iPoint)
    output[iPoint] /= std::sqrt(numTimePoints);

  return output;
}
} // namespace myFuncs
