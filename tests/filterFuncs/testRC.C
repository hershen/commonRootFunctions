#include "TGraph.h"
#include "filterFuncs.h"
#include "histFuncs.h"
#include "mathFuncs.h"
#include "myRootStyle.h"
#include <iostream>

const double dt = 0.1;
const double tau = 1;

void testRC() {
  myFuncs::setMyRootStyle();

  std::vector<double> zeros(0, 0);
  std::vector<double> ones(100, 1);
  std::vector<double> xs;
  xs.insert(xs.end(), zeros.begin(), zeros.end());
  xs.insert(xs.end(), ones.begin(), ones.end());

  std::vector<double> nom;
  std::vector<double> denom;

  std::tie(nom, denom) = myFuncs::DSP::getCR_RCnCoefficients(1, tau, 1/dt);
  const std::vector<double> filtered_RC = myFuncs::DSP::filter(nom, denom, xs);

  std::tie(nom, denom) = myFuncs::DSP::getCR_RCnCoefficients(2, tau, 1/dt);
  const std::vector<double> filtered_RC2 = myFuncs::DSP::filter(nom, denom, xs);


  std::tie(nom, denom) = myFuncs::DSP::getCR_RCnCoefficients(3, tau, 1/dt);
  const std::vector<double> filtered_RC3 = myFuncs::DSP::filter(nom, denom, xs);

  std::tie(nom, denom) = myFuncs::DSP::getCR_RCnCoefficients(4, tau, 1/dt);
  const std::vector<double> filtered_RC4 = myFuncs::DSP::filter(nom, denom, xs);

  // Draw
  TCanvas *canvas = new TCanvas("canvas", "", 0, 0, 1200, 900);
  std::vector<double> times;
  times.reserve(xs.size());
  for (auto i = 0; i < xs.size(); ++i)
    times.push_back(i * dt);

  TGraph *input = new TGraph(xs.size(), times.data(), xs.data());
  input->SetMarkerColor(kBlue);
  //
  TGraph *output = new TGraph(xs.size(), times.data(), filtered_RC.data());
  output->SetMarkerStyle(kFullCircle);
  output->SetMarkerColor(kGreen);
  TGraph *output2 = new TGraph(xs.size(), times.data(), filtered_RC2.data());
  output2->SetMarkerStyle(kFullSquare);
  output2->SetMarkerColor(kGreen + 1);
  TGraph *output4 = new TGraph(xs.size(), times.data(), filtered_RC4.data());
  output4->SetMarkerStyle(kFullCross);
  output4->SetMarkerColor(kGreen + 3);
  //

  // for(int i = 0; i < filtered.size(); ++i)
  // std::cout << "i = " << i << ", filtered = " << filtered[i] << std::endl;
  // std::cout << "size xs = " << xs.size() << ", size times = " << times.size() << ", size filtere = " <<
  // filtered.size() << std::endl;
  //
  TMultiGraph *multi = new TMultiGraph;
  // multi->Add(input);
  multi->Add(output);
  multi->Add(output2);
  multi->Add(output4);
  multi->SetTitle(";Time (arbitrary);Voltage (arbitrary)");
  multi->Draw("AP");

  myFuncs::Legend *legend = new myFuncs::Legend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(output, "CR-RC", "p");
  legend->AddEntry(output2, "CR-RC^{2}", "p");
  legend->AddEntry(output4, "CR-RC^{4}", "p");
  legend->Draw();

  // TCanvas *analyticalCanvas = new TCanvas("analyticalCanvas", "", 0, 0, 1200, 900);
  for (int n = 1; n <= 4; ++n) {
    if (n == 3)
      continue;
    TF1 *cr_rcn = new TF1(myFuncs::analytical_RC_CRn(n, tau, 1., -dt / 2.0, 300 * dt));
    // cr_rcn->SetLineColor(kGreen + n - 1);
    cr_rcn->SetNpx(1000);
    cr_rcn->Draw("SAME");
  }

  myFuncs::PaveText* pt = new myFuncs::PaveText(0.5);
  pt->AddText( ("Shaping of unit step with CR-RC^{n}, #tau=" + (boost::format("%g") %tau).str()).data() );
  pt->Draw();

  myFuncs::mySaveCanvas(canvas, "CR_RCnShaping");
}
