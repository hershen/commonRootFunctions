#pragma once

// STL
#include <memory>
#include <iostream>

// Root
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

// Mine
#include "histFuncs.h"
#include "mathFuncs.h"

namespace myFuncs {

template <class topType>
class ResidualCanvas {
public:
  ResidualCanvas(TCanvas& canvas, topType& topObject, TF1& function)
      : m_canvas(canvas), m_topObject(topObject), m_function(function) {

    if(!function.IsValid()) {
      std::cerr << "ResidualCanvas::ResidualCanvas: function isn't valid. Will probably crash soon...\n";
    }

    prepareCanvas();

    const auto residuals(myFuncs::calcResiduals(m_topObject, m_function));

    m_residualGraph = TGraph(myFuncs::getResidualsGraph(myFuncs::getXvalues(topObject), residuals));

    prepareObjects();
  }

  void draw() {
    m_topPad->cd();
    m_topObject.Draw("AP");
    m_function.Draw("Same");

    m_bottomPad->cd();
    m_residualGraph.Draw("AP");
  }

  TPad* getTopPad() const { return m_topPad; }
  TPad* getBottomPad() const { return m_bottomPad; }
  TGraph getResidualGraph() const { return m_residualGraph; }

private:
  TCanvas& m_canvas;
  topType& m_topObject;
  TF1& m_function;
  TGraph m_residualGraph;
  TPad* m_topPad;
  TPad* m_bottomPad;

  void prepareCanvas() {

    const double bottomPadYpercentage = 0.22;
    const double bottomTopSeperation = 0.05;

    m_canvas.cd();
    m_topPad = new TPad((m_canvas.GetName() + std::string("_topPad")).data(), "", 0, bottomPadYpercentage, 1, 1);
    m_topPad->SetNumber(1);
    m_topPad->Draw();

    // can't set margins with m_topPad->Set__Margin() for some reason. Have to go through m_canvas.cd(x)...
    m_canvas.cd(1)->SetBottomMargin(0.01);
    // Change to canvas before creating second pad
    m_canvas.cd();

    m_bottomPad = new TPad((m_canvas.GetName() + std::string("_bottomPad")).data(), "", 0, 0, 1, bottomPadYpercentage);
    m_bottomPad->SetNumber(2);

    m_bottomPad->Draw();
    m_canvas.cd(2)->SetTopMargin(bottomTopSeperation);
    m_canvas.cd(2)->SetBottomMargin(gStyle->GetPadBottomMargin() * 1.5);
    m_bottomPad->SetGridy();

    m_canvas.cd();
  }

  // Prepare objects for drawing in a prettyResidualGraph
  void prepareObjects() {
    m_topObject.GetXaxis()->SetLabelSize(0.0); // Don't draw x labels on top object
    // m_topObject.GetXaxis()->SetTitleSize(0.0); // Don't draw x title on top object

    // Set axis label and Title size to absolute
    m_topObject.GetYaxis()->SetLabelFont(43);     // Absolute font size in pixel (precision 3)
    m_residualGraph.GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    m_residualGraph.GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    // Title
    m_topObject.GetYaxis()->SetTitleFont(43);     // Absolute font size in pixel (precision 3)
    m_residualGraph.GetXaxis()->SetTitleFont(43); // Absolute font size in pixel (precision 3)
    m_residualGraph.GetYaxis()->SetTitleFont(43); // Absolute font size in pixel (precision 3)

    // Set x + y axis label size
    const double labelSize = std::min(0.03 * m_canvas.cd(1)->GetWh(), 0.03 * m_canvas.cd(1)->GetWw());
    m_topObject.GetYaxis()->SetLabelSize(labelSize);

    m_residualGraph.GetYaxis()->SetLabelSize(labelSize);
    // x axis
    m_residualGraph.GetXaxis()->SetLabelSize(labelSize);

    // Set axis title sizes
    const double titleSize = std::min(0.03 * m_canvas.cd(1)->GetWh(), 0.03 * m_canvas.cd(1)->GetWw());
    m_residualGraph.GetXaxis()->SetTitleSize(titleSize);
    m_residualGraph.GetYaxis()->SetTitleSize(titleSize);
    m_topObject.GetYaxis()->SetTitleSize(titleSize);

    // Set title offsets
    m_residualGraph.GetXaxis()->SetTitleOffset(3.75);

    // Set bottom x title
    m_residualGraph.GetXaxis()->SetTitle(m_topObject.GetXaxis()->GetTitle());
    // Set y title
    m_residualGraph.GetYaxis()->SetTitle("Pull (#sigma)");

    // Set residual y axis divisions
    const auto maxResidual =
        std::abs(*std::max_element(m_residualGraph.GetY(), m_residualGraph.GetY() + m_residualGraph.GetN() - 1,
                                   [](const double residual1, const double residual2) { // find max absolute value residual
                                     return std::abs(residual1) < std::abs(residual2);
                                   }));
    m_residualGraph.SetMaximum(std::ceil(maxResidual));
    m_residualGraph.SetMinimum(-std::ceil(maxResidual));
    const int maxDivisions = std::min(5., std::ceil(maxResidual));
    m_residualGraph.GetYaxis()->SetNdivisions(maxDivisions, false); // false - no optimization - forces current value

    // Set marker size
    double markerSize = myFuncs::linearInterpolate(100, 17500, 1.2, 0.1, m_residualGraph.GetN());
    // markerSize = std::min(static_cast<double>(gStyle->GetMarkerSize()), markerSize);
    m_residualGraph.SetMarkerSize(markerSize);
  }
};

} // namespace myFuncs
