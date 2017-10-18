#pragma once

#include "TStyle.h"

namespace myFuncs {

TStyle *MyRootStyle();

void setMyRootStyle();

void setRootStyle() {

  gStyle->SetOptStat("ourme");

  gStyle->SetOptFit(1);

  // Set ticks on upper and right axis (not just left and bottom)
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Colors for scales (for example in 2d histograms)
  gStyle->SetPalette(56);

  //   //Set default border outline and tick line width
  //   gStyle->SetLineWidth(2.);
  //
  //   //Set default border outline and tick line width
  //   //Didn't work for some reason
  //   gStyle->SetHistLineWidth(30.);
  //
  //   gStyle->SetStatBorderSize(1.);
}

} // namespace myFuncs
