#pragma once

// STL
#include <memory>

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
  ResidualCanvas(TCanvas &canvas, topType &topObject, TF1 &function)
      : m_canvas(canvas), m_topObject(topObject), m_function(function) {

    prepareCanvas();

    const auto residuals(myFuncs::calcResiduals(m_topObject, m_function));

    m_residualGraph = TGraph(myFuncs::getResidualsGraph(myFuncs::getXvalues(topObject), residuals));

    prepareObjects();
  }

  void draw() {
    m_canvas.cd(1);
    m_topPad->SetTopMargin(0.7);
    m_topObject.Draw();

    m_canvas.cd(2);
    m_residualGraph.Draw("AP");
  }

  TPad *getTopPad() const { return m_topPad; }
  TPad *getBottomPad() const { return m_bottomPad; }
  TGraph getResidualGraph() const { return m_residualGraph; }

private:
  TCanvas &m_canvas;
  topType &m_topObject;
  TF1 &m_function;
  TGraph m_residualGraph;
  TPad *m_topPad;
  TPad *m_bottomPad;

  void prepareCanvas() {
    const double topBottomSpace = 0.045;

    m_canvas.cd();
    m_topPad = new TPad((m_canvas.GetName() + std::string("_topPad")).data(), "", 0, 0.25, 1, 1);
    m_topPad->SetNumber(1);
    m_topPad->SetTopMargin(0.7); // gStyle->GetPadTopMargin());
    // std::cout << m_topPad->GetTopMargin() << std::endl;
    m_topPad->SetBottomMargin(0.01); // For some reason 0.0 doesn't look good
    m_topPad->Draw();

    // Change to canvas before creating second pad
    m_canvas.cd();

    m_bottomPad = new TPad((m_canvas.GetName() + std::string("_bottomPad")).data(), "", 0, 0.0, 1, 0.25);
    m_bottomPad->SetNumber(2);
    m_bottomPad->SetTopMargin(topBottomSpace);
    m_bottomPad->Draw();
    m_bottomPad->SetBottomMargin(0.7); // So that title isn't cut out
    m_bottomPad->SetGridy();

    m_canvas.cd();
  }

  // Prepare objects for drawing in a prettyResidualGraph
  void prepareObjects() {
    m_topObject.GetXaxis()->SetLabelSize(0.0); // Don't draw x labels on top object

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
    m_residualGraph.GetXaxis()->SetTitleOffset(3.25);

    // Set bottom x title
    m_residualGraph.GetXaxis()->SetTitle(m_topObject.GetXaxis()->GetTitle());
    // Set y title
    m_residualGraph.GetYaxis()->SetTitle("Residual (#sigma)");
  }
};

} // namespace myFuncs
