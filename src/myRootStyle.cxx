#include "myRootStyle.h"

//std
#include <iostream>

//Root
#include "TROOT.h"

void SetBelle2Style();

namespace myFuncs {

TStyle* MyRootStyle()
{
  TStyle *myRootStyle = new TStyle("myRootStyle","My Root style");

	// use plain black on white colors
  Int_t icol=0; // WHITE
  myRootStyle->SetFrameBorderMode(icol);
  myRootStyle->SetFrameFillColor(icol);
  myRootStyle->SetCanvasBorderMode(icol);
  myRootStyle->SetCanvasColor(icol);
  myRootStyle->SetPadBorderMode(icol);
  myRootStyle->SetPadColor(icol);
  myRootStyle->SetStatColor(icol);
	
  // set margin sizes
  myRootStyle->SetPadTopMargin(0.05);
  myRootStyle->SetPadRightMargin(0.05);
  myRootStyle->SetPadBottomMargin(0.16);
  myRootStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  myRootStyle->SetTitleXOffset(1.4);
  myRootStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  myRootStyle->SetTextFont(font);

  myRootStyle->SetTextSize(tsize);
  myRootStyle->SetLabelFont(font,"x");
  myRootStyle->SetTitleFont(font,"x");
  myRootStyle->SetLabelFont(font,"y");
  myRootStyle->SetTitleFont(font,"y");
  myRootStyle->SetLabelFont(font,"z");
  myRootStyle->SetTitleFont(font,"z");
  
  myRootStyle->SetLabelSize(tsize,"x");
  myRootStyle->SetTitleSize(tsize,"x");
  myRootStyle->SetLabelSize(tsize,"y");
  myRootStyle->SetTitleSize(tsize,"y");
  myRootStyle->SetLabelSize(tsize,"z");
  myRootStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  myRootStyle->SetMarkerStyle(20);
  myRootStyle->SetMarkerSize(1.2);
  myRootStyle->SetHistLineWidth(2.);
  myRootStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //myRootStyle->SetErrorX(0.001);
  // get rid of error bar caps
  myRootStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  myRootStyle->SetOptTitle(0);
  
  myRootStyle->SetOptStat("ourme");
  myRootStyle->SetOptFit(1111);
  
  // put tick marks on top and RHS of plots
  myRootStyle->SetPadTickX(1);
  myRootStyle->SetPadTickY(1);
	
	//Set default color of fitted functions
	myRootStyle->SetFuncColor(kRed);
	
	
	myRootStyle->SetMarkerColor(kBlue);
	
	
	myRootStyle->SetLineColor(kBlue);
	
  return myRootStyle;

}

void setMyRootStyle()
{
  static TStyle* myRootStyle = nullptr;
  std::cout << "\nApplying myRootStyle\n" << std::endl ;
  if ( !myRootStyle ) myRootStyle = MyRootStyle();
  gROOT->SetStyle("myRootStyle");
  gROOT->ForceStyle();
}

}//namespace