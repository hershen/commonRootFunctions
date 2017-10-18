#include "testbeam/Waveform.h"
#include <iostream>

// ROOT
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMath.h"

using namespace myFuncs::testbeam;

Waveform::Waveform(const std::vector<uint32_t> samples, const double dt) : m_samples(samples), m_dt(dt) {}

double Waveform::getStd(const unsigned first, const unsigned last) const {
  return TMath::RMS(m_samples.begin() + first, m_samples.begin() + last);
}

double Waveform::getMean(const unsigned int first, const unsigned int last) const {
  return TMath::Mean(m_samples.begin() + first, m_samples.begin() + last);
}

std::vector<double> &Waveform::getSamplesDouble() const {

  // If already calculated, return vector
  if (m_samples_double.size() == 0) {
    // Transform original samples to double
    m_samples_double.reserve(m_samples.size());
    for (const auto sample : m_samples)
      m_samples_double.push_back(static_cast<double>(sample));
  }

  return m_samples_double;
}

std::vector<double> &Waveform::getTimes() const {

  // If already calculated, return vector
  if (m_times.size() == 0) {
    // Prepare times vector
    m_times.reserve(m_samples.size());
    for (uint iSample = 0; iSample < m_samples.size(); ++iSample)
      m_times.push_back(iSample * m_dt);
  }

  return m_times;
}

std::pair<double, double> Waveform::getMaxPoly2(const unsigned int first, const unsigned int last) const {

  // Sanity
  if (first >= last)
    return {0.0, 0.0};

  const double std = getStd();
  const long numPoints = last - first;

  std::vector<double> errors(numPoints, std);

  TGraphErrors graph(numPoints, getTimes().data() + first, getSamplesDouble().data() + first, 0, errors.data());
  auto fitResult = graph.Fit("pol2", "MSQ");

  const double timeAtMax = -fitResult->Value(1) / fitResult->Value(2) / 2.0;
  const double maxVal = fitResult->Value(0) + fitResult->Value(1) * timeAtMax + fitResult->Value(2) * timeAtMax * timeAtMax;

  // 	std::cout << "Chi2 / ndf = " << fitResult->Chi2() << "/" << fitResult->Ndf() << ", prob = " << fitResult->Prob() <<
  // std::endl; 	std::cout << "time - (" << 	getTimes()[first] << " , " << getTimes()[first + numPoints] << std::endl;
  //

  return {timeAtMax, maxVal};
}

double Waveform::getSimpleAmplitude() const {
  double amp = -1.0;

  const double pedestal = getMean();
  for (unsigned int idx = m_samples.size() * 0.4; idx < m_samples.size(); ++idx) {
    if ((m_samples[idx] - pedestal) > amp) {
      amp = m_samples[idx] - pedestal;
		}
  }

  return amp;
}

TGraphErrors Waveform::getGraphErrors() {

  std::vector<double> stds(m_samples.size(), getStd());
  TGraphErrors graph(m_samples.size(), getTimes().data(), getSamplesDouble().data(), 0, stds.data());

  graph.SetTitle(";Time (ns); Voltage (ADC counts)");
  graph.SetMarkerSize(0.5);

  return graph;
}

TH1D Waveform::getHistWithErrors() {

  const double std = getStd();

  TH1D fitHist("fitHist", "fitHist", m_samples.size(), -m_dt / 2.0, m_samples.size() * m_dt - m_dt / 2.0);
  for (size_t idx = 1; idx <= m_samples.size(); ++idx) {
    fitHist.SetBinContent(idx, m_samples[idx - 1]);
    fitHist.SetBinError(idx, std);
  }

  fitHist.SetTitle(";Time (ns); Voltage (ADC counts)");

  return fitHist;
}
