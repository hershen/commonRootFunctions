#pragma once

// std
#include <iostream>
#include <numeric> //for iota

// Boost
#include "boost/format.hpp"

// Root
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"

namespace myFuncs {

// Draw histograms so that the one with the largest maximum is drawn first and the rest afterwards on the same canvas
void drawHistograms_highestFirst(const std::vector<TH1 *> &histVector);

TCanvas *newCanvas(std::string canvasName = "");

// Draw histogram on new canvas. Canvas name is <histname>Canvas
TCanvas *drawNewCanvas(TH1 *hist);

// Draw graph on new canvas. Canvas name is <graph>Canvas. Draw option is "AP"
TCanvas *drawNewCanvas(TGraph *graph);

// Add a label above bin number bin
// Takes care of logy, but not of logx pads
void putLabelAboveBin(const TH1 *hist, const size_t bin, const std::string &text, const double textHeight = 0.0,
                      const double yOffset = 0.0);

// Compute the integral of a TGraph. Linear interpolation is assumed between the points.
// First sorts the points in ascending order, so they don't have to be sorted.
double integrateTGraph(const TGraph &graph);

// Normalize the TGraph so that the integral = area.
// If area == 0.0 do nothing
// If integral of the tgraph == 0.0 do nothing
void normalize(TGraph &graph, const double area = 1.0);

// Normalize histogram to area
// If area == 0.0 do nothing
// If integral of histogram with options == 0.0 do nothing
void normalize(TH1D &hist, const double area = 1.0, const std::string &options = "width");

// Get bin width with the specified precision".
inline const std::string getBinWidth(const TH1 &hist, const std::string &precision = "g") {
  return boost::str(boost::format("%1$" + precision) % hist.GetBinWidth(1));
}

/**
 * Extract a root object rootObjName from a root file file. The file is assumed to be not zombie. If rootObjName doesn't exist in
 * file, print error. rootClass should be a pointer (unless you figure out how to use it otherwise).
 */
template <class rootClass>
rootClass getRootObjectFromFile(TFile &file, const std::string &rootObjName) {
  rootClass rootObj = static_cast<rootClass>(file.Get(rootObjName.data()));
  if (!rootObj) {
    std::string filename = file.GetName();
    std::cout << "myFuncs::getRootObjectFromFile: Could not find " << rootObjName << " in " << filename << std::endl;
    ;
  }
  return rootObj;
}

// Bin a TGraph into a hist. Hist should already have the required bins.
// The content of each bin will be the average y values of the points corresponding to that bin.
// If errorOnMean is false, the y error for each point will be the population STD in that bin.
// If errorOnMean is true, the y error for each point will be the (population STD in that bin) / sqrt(number of entries in the
// bin).
void binGraph(const TGraph &graph, TH1 &hist, const bool errorOnMean = false);

// Convert a histogram to a TGraphErrors.
// Error on x is the bin width.
// Error on y is the bin error.
TGraphErrors histToGraph(const TH1 &hist);

// template <class histClass>
// histClass scaleXaxis(const histClass& hist, const double scale)
// {
//   histClass newHist = *(static_cast<histClass*> (hist.Clone()) );
// 	newHist.Clear();
//   for(iBin = 1; iBin <= hist.GetNbinsX(); ++iBin)
// 		newHist.
//
//   return newHist;
//
// }

// Save canvas as fileType, inside figureDump folder.
// Overwrites files if they exist
inline void mySaveCanvas(const TCanvas *canvas, const std::string &filename, const std::string &fileType) {
  // Create figure dump dir. If already exists, or there's a problem, it returns -1.
  gSystem->mkdir("figureDump");
  canvas->SaveAs(("figureDump/" + filename + "." + fileType).data());
}

// Save canvas as png, pdf, C, inside figureDump folder.
// Overwrites files if they exist
inline void mySaveCanvas(const TCanvas *canvas, const std::string &filename) {
  mySaveCanvas(canvas, filename, "png");
  mySaveCanvas(canvas, filename, "pdf");
  mySaveCanvas(canvas, filename, "C");
}

inline void mySaveCanvas(const TCanvas *canvas) { mySaveCanvas(canvas, canvas->GetName()); }

class PaveText : public TPaveText {
public:
  PaveText(const Double_t x1, const Double_t y1, const Double_t x2, const Double_t y2, Option_t *option = "NDCNB")
      : TPaveText(x1, y1, x2, y2, option) {
    this->SetTextFont(22);
    //     this->SetFillColor(kWhite);
    this->SetFillStyle(0);  // Make fill color transparent
    this->SetBorderSize(0); // No border
  }

  // Constructor which puts text at the top of a histogram. Just need to set x2
  // Numbers for x1,y1,y2 assume bottom and left pad margin = 0.16, top pad margin = 0.05.
  // TODO - think if might be accessed with gPad
  PaveText(const Double_t x2, Option_t *option = "NDCNB")
      : PaveText((gPad ? gPad->GetLeftMargin() : gStyle->GetPadLeftMargin()), 0.95, x2, 1, option) {
    this->SetTextAlign(12);
  }
};

// Nice legend
class Legend : public TLegend {
public:
  Legend(const Double_t x1, const Double_t y1, const Double_t x2, const Double_t y2) : TLegend(x1, y1, x2, y2) {
    this->SetTextFont(22);
    this->SetFillStyle(0);  // Transpartent fill color
    this->SetBorderSize(1); // No shadow
  }
};

template <class T>
TGraph getResidualsGraph(const std::vector<T> &xValues, const std::vector<double> residuals) {
  if (xValues.size() != residuals.size())
    throw;

  TGraph residualsGraph(xValues.size(), xValues.data(), residuals.data());

  const auto maxResidual =
      std::abs(*std::max_element(residuals.begin(), residuals.end(),
                                 [](const T residual1, const T residual2) { return std::abs(residual1) < std::abs(residual2); }));
  residualsGraph.SetMaximum(std::ceil(maxResidual));
  residualsGraph.SetMinimum(-std::ceil(maxResidual));
  residualsGraph.GetYaxis()->SetNdivisions(std::ceil(maxResidual), false); // false - no optimization - forces current value

  return residualsGraph;
}

double FWHM(const TH1 &hist);

inline double FWHM_2355(const TH1 &hist) { return FWHM(hist) / 2.355; }

// Calculate the chi of the overlap between 2 histograms.
// chi2 = sum (hist0_binVal - hist1_binVal)^2/(hist0_binError^2 + hist1_binError^2).
// Sum is over the bins in the range that corresponds to the range (minVal, maxVal).
// minVal, maxVal are not bin numbers - they are the x axis value.
// Currently supports only histograms with same number of bins in the range.
double getHistOverlapChi2(TH1D hist0, TH1D hist1, const double minVal, const double maxVal);

// Set novosibirsk parameters, assuming order is normalization, peak, width, eta.
// Eta is guessed at 0.2
void setNovosibirskParams(TF1 &novosibirskTF1, const TH1 &hist);

inline double xUserToNdc(const double x) {
  gPad->Update(); // this is necessary!
  return (x - gPad->GetX1()) / (gPad->GetX2() - gPad->GetX1());
}

inline double yUserToNdc(const double y) {
  gPad->Update(); // this is necessary!
  return (y - gPad->GetY1()) / (gPad->GetY2() - gPad->GetY1());
}

inline double xNdcToUser(const double x) {
  gPad->Update(); // this is necessary!
  return x * (gPad->GetX2() - gPad->GetX1()) + gPad->GetX1();
}

inline double yNdcToUser(const double y) {
  gPad->Update(); // this is necessary!
  return y * (gPad->GetY2() - gPad->GetY1()) + gPad->GetY1();
}

template <class T>
void drawTObjects(T &object) {
  object.Draw();
}
template <class T, class... Ts>
void drawTObjects(T &first, Ts &... others) {
  first.Draw();
  drawTObjects(others...);
}

inline std::vector<double> getXvalues(TH1 &hist) {

  std::vector<double> xValues;
  xValues.reserve(hist.GetXaxis()->GetNbins());

  for (uint binNum = 0; binNum < xValues.size(); ++binNum) {
    xValues.push_back(hist.GetBinCenter(binNum));
  }

  return xValues;
}

inline std::vector<double> getXvalues(TGraph &graph) { return std::vector<double>(graph.GetX(), graph.GetX() + graph.GetN()); }

// Wait in current pad for double click.
// If there is no pad, or if pad is closed while waiting, terminates the root TApplication.
void waitForDoubleClick();

// Produce TGraph from yValues and set title.
// xValues will be 0, 1, 2,...
template <class T>
TGraph getGraph(const std::vector<T> &xValues, const std::vector<T> yValues, const std::string &title) {
  if (xValues.size() != yValues.size()) {
    throw std::invalid_argument("histFuncs::getGraph: xValues.size() = " + std::to_string(xValues.size()) +
                                ", yValues.size() = " + std::to_string(yValues.size()));
  }
  TGraph graph(yValues.size(), xValues.data(), yValues.data());
  graph.SetTitle(title.c_str());
  
  return graph;
}

// Produce TGraph from yValues and set title.
// xValues will be 0, 1, 2,...
template <class T>
TGraph getGraph(const std::vector<T> &yValues, const std::string &title) {
  std::vector<T> xValues(yValues.size());
  std::iota(xValues.begin(), xValues.end(), T(0));
  return getGraph(xValues, yValues, title);
}

} // namespace myFuncs
