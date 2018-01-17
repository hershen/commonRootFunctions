#pragma once

// STL
#include <cstdint>
#include <iterator>
#include <memory>
#include <vector>
// ROOT
// #include "TFitResult.h"
// #include "TGraphErrors.h"
// #include "TH1D.h"
#include "TMath.h"

// Mine
#include "mathFuncs.h"

class TGraph;
class TGraphErrors;
class TH1D;

namespace myFuncs {
namespace testbeam {

template <class T>
class Waveform {
public:
  Waveform(const std::vector<T> &samples, const double dt) : m_samples(samples), m_dt(dt) {
    if (m_dt <= 0.0) {
      throw std::invalid_argument("testbeam::Waveform::Waveform: dt = " + std::to_string(m_dt) + " <= 0");
    }
  }

  // Get (sample!) standard deviation between first and last
  // No bounds checking
  inline double getStd(const size_t firstIdx, const size_t lastIdx) const {
    if (lastIdx - firstIdx == 0) {
      return 0.0;
    }

    return sampleStd(std::next(m_samples.begin(), firstIdx), std::next(m_samples.begin(), lastIdx + 1));
  }

  inline double getStd() const { return myFuncs::sampleStd(m_samples.begin(), m_samples.end()); }

  // Get (sample!) mean between first and last
  double getMean(const size_t firstIdx, const size_t lastIdx) const {
    if (lastIdx - firstIdx == 0) {
      return 0.0;
    }
    return myFuncs::sampleMean(std::next(m_samples.begin(), firstIdx), std::next(m_samples.begin(), lastIdx + 1));
  }

  inline double getMean() const { return myFuncs::sampleMean(m_samples.begin(), m_samples.end()); }

  inline const std::vector<T> &getSamples() const { return m_samples; }

  inline double getDt() const { return m_dt; }

  Waveform<double> transformToDouble() const {
    std::vector<double> samplesDouble;
    samplesDouble.reserve(m_samples.size());
    std::for_each(m_samples.begin(), m_samples.end(),
                  [&](const T element) { samplesDouble.push_back(static_cast<double>(element)); });
    return Waveform<double>(samplesDouble, getDt());
  }

  // Get vector of times.
  const std::vector<double> &getTimes() const {
    // If already calculated, return vector
    if (m_times.empty()) {
      fillTimes();
    }
    return m_times;
  }

  // Return pair of maximum idx and maximum value.
  // Searches every every'th element, then again from (maximum found - every + 1, maximum found + every -1)
  std::pair<size_t, T> getMaximumIdx_value(const size_t every = 20) const {
    const auto maxItr = myFuncs::findMaxEvery_n_ThenBetween(m_samples.begin(), m_samples.end(), every);
    if (maxItr == m_samples.end()) {
      return std::make_pair(0, 0.0);
    }
    return std::make_pair(std::distance(m_samples.begin(), maxItr), *maxItr);
  }

  // Access element idx.
  // No bounds checking!
  T &operator[](size_t idx) { return m_samples[idx]; }
  const T operator[](size_t idx) const { return m_samples[idx]; }

  // Average each n samples (n=0,1 returns identity)
  // Returns new waveform
  Waveform<double> averageEach_n(const size_t n) const {
    return Waveform<double>(myFuncs::averageEach_n(m_samples, n), getDt() * static_cast<double>(n));
  }

  Waveform<double> timeShift(const double shift) const {
    return Waveform<double>(myFuncs::shiftVector(m_samples, shift, getDt()), getDt());
  }

  template <class newWaveformType>
  Waveform<decltype(T() + newWaveformType())> operator+(const Waveform<newWaveformType> &rWaveform) const {
    if (getDt() != rWaveform.getDt()) {
      throw std::invalid_argument("Waveform::operator+: own dt = " + std::to_string(getDt()) +
                                  ", received waveform dt = " + std::to_string(rWaveform.getDt()));
    }
    if (m_samples.size() != rWaveform.getSamples().size()) {
      throw std::invalid_argument("Waveform::operator+: own number of elements = " + std::to_string(m_samples.size()) +
                                  ", received waveform number of elements = " + std::to_string(rWaveform.getSamples().size()));
    }

    using returnType = decltype(T() + newWaveformType());

    std::vector<returnType> output;
    std::transform(m_samples.begin(), m_samples.end(), rWaveform.getSamples().begin(), std::back_inserter(output),
                   [](const T own, const newWaveformType theirs) { return own + theirs; });
    return Waveform<returnType>(output, getDt());
  }

private:
  std::vector<T> m_samples;
  mutable std::vector<double> m_times;
  const double m_dt;

  void fillTimes() const {
    // Prepare times vector
    m_times.reserve(m_samples.size());
    for (uint iSample = 0; iSample < m_samples.size(); ++iSample)
      m_times.push_back(iSample * m_dt);
  }

  template <typename Iterator>
  inline double calcStd(Iterator first, Iterator last) const {
    return TMath::StdDev(first, last);
  }
};

} // namespace testbeam
} // namespace myFuncs
