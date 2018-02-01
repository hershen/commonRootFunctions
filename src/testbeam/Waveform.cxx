#include "testbeam/Waveform.h"

// STL
#include <iostream>

namespace myFuncs {
namespace testbeam {

Waveform::Waveform(const std::vector<double> &samples, const double dt) : m_samples(samples), m_dt(dt) {
  if (m_dt <= 0.0) {
    throw std::invalid_argument("testbeam::Waveform::Waveform: dt = " + std::to_string(m_dt) + " <= 0");
  }
}

double Waveform::getStd(const size_t firstIdx, const size_t lastIdx) const {
  if (lastIdx - firstIdx == 0) {
    return 0.0;
  }

  return sampleStd(std::next(m_samples.begin(), firstIdx), std::next(m_samples.begin(), lastIdx + 1));
}

double Waveform::getMean(const size_t firstIdx, const size_t lastIdx) const {
  if (lastIdx - firstIdx == 0) {
    return 0.0;
  }
  return myFuncs::sampleMean(std::next(m_samples.begin(), firstIdx), std::next(m_samples.begin(), lastIdx + 1));
}

std::pair<size_t, double> Waveform::getMaximumIdx_value(const size_t every) const {
  const auto maxItr = myFuncs::findMaxEvery_n_ThenBetween(m_samples.begin(), m_samples.end(), every);
  if (maxItr == m_samples.end()) {
    return std::make_pair(0, 0.0);
  }
  return std::make_pair(std::distance(m_samples.begin(), maxItr), *maxItr);
}

Waveform Waveform::operator+(const Waveform &rWaveform) const {
  if (getDt() != rWaveform.getDt()) {
    throw std::invalid_argument("Waveform::operator+: own dt = " + std::to_string(getDt()) +
                                ", received waveform dt = " + std::to_string(rWaveform.getDt()));
  }
  if (m_samples.size() != rWaveform.getSamples().size()) {
    throw std::invalid_argument("Waveform::operator+: own number of elements = " + std::to_string(m_samples.size()) +
                                ", received waveform number of elements = " + std::to_string(rWaveform.getSamples().size()));
  }

  std::vector<double> output;
  std::transform(m_samples.begin(), m_samples.end(), rWaveform.getSamples().begin(), std::back_inserter(output),
                 std::plus<double>());
  return Waveform(output, getDt());
}

} // namespace testbeam
} // namespace myFuncs
