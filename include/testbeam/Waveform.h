#pragma once

// STL
#include <cstdint>
#include <memory>
#include <vector>

// Mine
#include "filterFuncs.h"
#include "mathFuncs.h"

class TGraph;
class TGraphErrors;
class TH1D;

namespace myFuncs {
namespace testbeam {

class Waveform {
public:
  Waveform() = default;

  Waveform(const std::vector<double> &samples, const double dt);

  template <class T>
  Waveform(const std::vector<T> &samples, const double dt) : Waveform(std::vector<double>(), dt) {
    m_samples.reserve(m_samples.size());
    std::for_each(samples.begin(), samples.end(), [&](const T element) { m_samples.push_back(static_cast<double>(element)); });
  }

  // Get (sample!) standard deviation between first and last
  // No bounds checking
  double getStd(const size_t firstIdx, const size_t lastIdx) const;

  inline double getStd() const { return myFuncs::sampleStd(m_samples.begin(), m_samples.end()); }

  // Get (sample!) mean between first and last
  // No bounds checking on idices
  double getMean(const size_t firstIdx, const size_t lastIdx) const;

  inline double getMean() const { return myFuncs::sampleMean(m_samples.begin(), m_samples.end()); }

  inline const std::vector<double> &getSamples() const { return m_samples; }

  inline void setSamples(const std::vector<double> newSamples) { m_samples = newSamples; }

  inline double getDt() const { return m_dt; }

  // Return pair of maximum idx and maximum value.
  // Searches every every'th element, then again from (maximum found - every + 1, maximum found + every -1)
  std::pair<size_t, double> getMaximumIdx_value(const size_t every = 20) const;

  // Access element idx.
  // No bounds checking!
  double &operator[](size_t idx) { return m_samples[idx]; }
  double operator[](size_t idx) const { return m_samples[idx]; }

  // Average each n samples (n=0,1 does nothing)
  inline void averageEach_n(const size_t n) {
    m_samples = myFuncs::averageEach_n(m_samples, n);
    m_dt = getDt() * static_cast<double>(n);
  }

  inline void timeShift(const double shift) { m_samples = myFuncs::shiftVector(m_samples, shift, getDt()); }

  Waveform operator+(const Waveform &rWaveform) const;

  // Remove given pedestal from waveform
  void removePedestal(const double pedestal) { m_samples = myFuncs::addToVector(m_samples, -pedestal); }

  // Get maxelement - average of first elementsForPedestal elements
  // No bounds checking on elementsForPedestal
  inline double getMaxMinusPedestal(const size_t elementsForPedestal) const {
    return getMaximumIdx_value().second - getMean(0, elementsForPedestal - 1);
  }

  inline size_t size() const { return m_samples.size(); }

private:
  std::vector<double> m_samples;
  double m_dt;
};

} // namespace testbeam

namespace DSP {
// Untested
template <class coefficientType>
void filter(const std::vector<coefficientType> &originalNominators, const std::vector<coefficientType> &originalDenominators,
            myFuncs::testbeam::Waveform &waveform) {
  waveform.setSamples(myFuncs::DSP::filter(originalNominators, originalDenominators, waveform.getSamples()));
}

// Untested
template <class coefficientType>
void filter(const std::pair<std::vector<coefficientType>, std::vector<coefficientType>> &coeficientPair,
            myFuncs::testbeam::Waveform &waveform) {
  waveform.setSamples(myFuncs::DSP::filter(coeficientPair, waveform.getSamples()));
}
} // namespace DSP

} // namespace myFuncs
