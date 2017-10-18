#pragma once

// STL
#include <cstdint>
#include <memory>
#include <vector>

class TGraphErrors;
class TH1D;

namespace myFuncs {
namespace testbeam {

class Waveform {
public:
  Waveform(const std::vector<uint32_t> samples, const double dt);

  // Get standard deviation between first and last
  double getStd(const unsigned int first, const unsigned int last) const;

  // Overloaded - Get standard deviation between m_samples[0] and m_samples[size * 12%]
  inline double getStd() const { return getStd(0, m_samples.size() * 0.12); }

  // Get mean between first and last
  double getMean(const unsigned int first, const unsigned int last) const;

  // Overloaded - Get mean deviation between m_samples[0] and m_samples[size * 12%]
  double getMean() const { return getMean(0, m_samples.size() * 0.12); }

  inline const std::vector<uint32_t> &getSamples() const { return m_samples; }
  std::vector<double> &getSamplesDouble() const;

  // Get vector of times.
  std::vector<double> &getTimes() const;

  // Return time at maximum and maximum.
  // Calculated with 2nd degree polynomial between first and last
  std::pair<double, double> getMaxPoly2(const unsigned int first, const unsigned int last) const;

  // Overloaded - First and last calculated 60-80%
  inline std::pair<double, double> getMaxPoly2() const { return getMaxPoly2(m_samples.size() * 0.6, m_samples.size() * 0.8); }

  // Get simple amplitude = maximum sample - pedestal
  // Pedestal is taken to be getMean()
  // Because of this, it's not the most efficient because it loops on values again.
  double getSimpleAmplitude() const;

  // Produce a TGraphErros from the wavefrom.
  // x axis errors are 0.
  // y axis errors are a constant getStd().
  TGraphErrors getGraphErrors();

  // Produce a TH1D from the wavefrom.
  // All bin errors aer a constant getStd().
  // Seems to be about 30% - 40% slower than using getGraphErrors().
  TH1D getHistWithErrors();

private:
  std::vector<uint32_t> m_samples;
  mutable std::vector<double> m_samples_double;
  mutable std::vector<double> m_times;
  const double m_dt;
};

} // namespace testbeam
} // namespace myFuncs
