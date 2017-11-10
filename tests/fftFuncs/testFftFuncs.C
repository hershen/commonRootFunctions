#include "fftFuncs.h"
#include "mathFuncs.h"
#include "myRootStyle.h"

bool testRealSequence2psd(const size_t N) {

  std::vector<double> x;
  x.reserve(N);

  TF1 cosFunc(
      "cosFunc",
      "1 + TMath::Cos(2.0*TMath::Pi()/2.0 * x) +  TMath::Cos(2.0*TMath::Pi()/5.0 * x) + TMath::Cos(2.0*TMath::Pi()/10.0 * x)", 0,
      N);

  // Create cosine samples
  for (auto i = 0; i < N; ++i) {
    x.push_back(cosFunc.Eval(i));
  }

  const auto psd = myFuncs::realSequence2psd(x);

  if (psd.size() != (N / 2 + 1)) {
    std::cerr << "testRealSequence2psd: Failed - size of psd is " << psd.size() << ", expected " << N / 2 + 1 << std::endl;
    return false;
  }

  const double xSquared = myFuncs::sumVectorSquared(x);
  const double psdSum = myFuncs::sumVector(psd);

  if (std::abs(xSquared - psdSum) > 1e-9) {
    std::cerr << "testRealSequence2psd: Perseval's theorem failed. sum(x[i]*x[i]) =  " << xSquared
              << ", sum(psd[i]) (already squared) =  " << psdSum << std::endl;
    return false;
  }

  // Check individual entries
  if (N == 100) {
    for (size_t i = 0; i < psd.size(); ++i) {

      if (i == 0) {
        if (std::abs(psd[i] - 100) > 1e-9) {
          std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 100 " << std::endl;
          return false;
        }
      } else if (i == 10) {
        if (std::abs(psd[i] - 50) > 1e-9) {
          std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 50 " << std::endl;
          return false;
        }
      } else if (i == 20) {
        if (std::abs(psd[i] - 50) > 1e-9) {
          std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 50 " << std::endl;
          return false;
        }
      } else if (i == 50) {
        if (std::abs(psd[i] - 100) > 1e-9) {
          std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 100 " << std::endl;
          return false;
        }
      } else if (std::abs(psd[i]) > 1e-9) {
        std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 0 " << std::endl;
        return false;
      }
    }
  }

  return true;
}

void draw() {
  /* Produces a waveform containing 2 frequencies.
  Filters out 1 frequency by zeroing the fft values.
  Reconstructs the time domain waveform of the filtered signal */

  myFuncs::setMyRootStyle();
  const double f1 = 5e-2; // 50 MHz
  const double f2 = 1e-2; // 10 MHz
  const double amp1 = 1;
  const double amp2 = 2;

  const size_t numPoints = 1e2;
  const double dt = 2.0; // ns;
  TF1 *func = new TF1("func", "[0]*TMath::Sin(2*TMath::Pi()*[1]*x) + [2]*TMath::Sin(2*TMath::Pi()*[3]*x)", 0, numPoints * dt);
  func->SetParameters(amp1, f1, amp2, f2);
  func->SetNpx(1000);
  func->SetLineColor(kGreen);

  std::vector<double> measurements;
  measurements.reserve(numPoints);
  std::vector<double> times;
  times.reserve(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    measurements.push_back(func->Eval(i * dt));
    times.push_back(i * dt);
  }

  TGraph *graph = new TGraph(numPoints, times.data(), measurements.data());
  TCanvas *timeCanvas = new TCanvas("timeCanvas", "", 0, 0, 1200, 900);
  graph->SetTitle(";Time (Arbitrary);Voltage (Arbitrary)");
  graph->Draw("AP");
  func->Draw("Same");

  const auto frequencies = myFuncs::getRealFftfrequencies(measurements.size(), 1. / dt);

  auto fftFiltered = myFuncs::fftR2C(measurements);
  for (size_t i = 0; i < frequencies.size(); ++i) {
    if (std::abs(frequencies[i] - f2) < 1e-9) {
      fftFiltered.first[i] = 0.0;
      fftFiltered.second[i] = 0.0;
    }
  }
  const auto filtered = myFuncs::fftC2R(measurements.size(), fftFiltered.first, fftFiltered.second);
  TGraph *graphFiltered = new TGraph(numPoints, times.data(), filtered.data());
  graphFiltered->SetMarkerColor(kRed);
  graphFiltered->Draw("SAMEP");
  TF1 *sineFunc = new TF1("sineFunc", "[0]*TMath::Sin(2*TMath::Pi()*[1]*x)", 0, numPoints * dt);
  sineFunc->SetParameters(amp1, f1);
  sineFunc->SetNpx(1000);
  sineFunc->Draw("SAME");

  // PSD
  const auto psdUnfiltered = myFuncs::realSequence2psd(measurements, 1. / dt);
  const auto psdFiltered = myFuncs::realSequence2psd(filtered, 1. / dt);

  std::cout << "filtered.size() = " << filtered.size() << "\n";
  std::cout << "measurements.size() = " << measurements.size() << "\n";

  TGraph *psdUnfilteredGraph = new TGraph(frequencies.size(), frequencies.data(), psdUnfiltered.data());
  TGraph *psdFilteredGraph = new TGraph(frequencies.size(), frequencies.data(), psdFiltered.data());
  psdFilteredGraph->SetMarkerColor(kRed);
  TCanvas *psdCanvas = new TCanvas("psdCanvas", "", 0, 0, 1200, 900);
  TMultiGraph *multi = new TMultiGraph;
  multi->Add(psdUnfilteredGraph);
  multi->Add(psdFilteredGraph);
  multi->SetTitle(";Frequency (Arbitrary);PSD (Power/Frequency)");
  multi->Draw("AP");
}

void testFftFuncs() {

  // The 100 is important because we check individual results of the psd
  if (!testRealSequence2psd(100)) {
    std::cerr << "testRealSequence2psd with input 100 Failed!!!" << std::endl;
  }

  if (!testRealSequence2psd(101)) {
    std::cerr << "testRealSequence2psd with input 101 Failed!!!" << std::endl;
  }

  draw();
}
