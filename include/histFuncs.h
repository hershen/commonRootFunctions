#pragma once

//std
#include <iostream>

//Boost
#include "boost/format.hpp"

//Root
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TLegend.h"

namespace myFuncs
{
  
//Draw histograms so that the one with the largest maximum is drawn first and the rest afterwards on the same canvas
void drawHistograms_highestFirst(const std::vector<TH1*> &histVector);

TCanvas* newCanvas(std::string canvasName = "");

//Draw histogram on new canvas. Canvas name is <histname>Canvas
TCanvas* drawNewCanvas(TH1* hist);

//Draw graph on new canvas. Canvas name is <graph>Canvas. Draw option is "AP"
TCanvas* drawNewCanvas(TGraph* graph);

//Add a label above bin number bin
//Takes care of logy, but not of logx pads
void putLabelAboveBin(const TH1* hist, const size_t bin, const std::string& text, const double textHeight = 0.0, const double yOffset = 0.0);


//Compute the integral of a TGraph. Linear interpolation is assumed between the points.
//First sorts the points in ascending order, so they don't have to be sorted.
double integrateTGraph(const TGraph& graph);

//Normalize the TGraph so that the integral = area.
//If area == 0.0 do nothing
//If integral of the tgraph == 0.0 do nothing
void normalize(TGraph& graph, const double area = 1.0);

//Normalize histogram to area
//If area == 0.0 do nothing
//If integral of histogram with options == 0.0 do nothing
void normalize(TH1D& hist, const double area = 1.0, const std::string& options = "width");


//Get bin width with the specified precision".
inline const std::string getBinWidth(const TH1& hist, const std::string& precision = ".3g") { return boost::str(boost::format("%1$" + precision)% hist.GetBinWidth(1) ); }

/**
* Extract a root object rootObjName from a root file file. The file is assumed to be not zombie. If rootObjName doesn't exist in file, print error.
* rootClass should be a pointer (unless you figure out how to use it otherwise).
*/
template <class rootClass> rootClass getRootObjectFromFile(TFile& file, const std::string& rootObjName)
{
	rootClass rootObj = static_cast<rootClass>( file.Get(rootObjName.data() ) );
	if (!rootObj) 
	{
		std::string filename = file.GetName();
		std::cout << "myFuncs::getRootObjectFromFile: Could not find " << rootObjName << " in " << filename << std::endl;;
	}
	return rootObj;
}

//Bin a TGraph into a hist. Hist should already have the required bins.
//The content of each bin will be the average y values of the points corresponding to that bin.
//If errorOnMean is false, the y error for each point will be the population STD in that bin.
//If errorOnMean is true, the y error for each point will be the (population STD in that bin) / sqrt(number of entries in the bin).
void binGraph(const TGraph& graph, TH1& hist, const bool errorOnMean = false);

//Convert a histogram to a TGraphErrors.
//Error on x is the bin width.
//Error on y is the bin error.
TGraphErrors histToGraph(const TH1& hist);

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


// Overwrites files if they already exis
inline void mySaveCanvas(const TCanvas* canvas, const std::string& filename)
{
	//Create figure dump dir. If already exists, or there's a problem, it returns -1.
	gSystem->mkdir("figureDump");
	
	canvas->SaveAs(("figureDump/" + filename + ".png").data());
	canvas->SaveAs(("figureDump/" + filename + ".pdf").data());
	canvas->SaveAs(("figureDump/" + filename + ".C").data());
}


inline void mySaveCanvas(const TCanvas* canvas) {
	mySaveCanvas(canvas, canvas->GetName());	
}

class PaveText : public TPaveText
{
public:
	PaveText (const Double_t x1, const Double_t y1, const Double_t x2, const Double_t y2, Option_t *option="NDCNB"):
	TPaveText(x1, y1, x2, y2, option)
	{
		this->SetTextFont(22);
//     this->SetFillColor(kWhite);
		this->SetFillStyle(0);  //Make fill color transparent
		this->SetBorderSize(0); //No border
	}
	
	//Constructor which puts text at the top of a histogram. Just need to set x2
	//Numbers for x1,y1,y2 assume bottom and left pad margin = 0.16, top pad margin = 0.05.
	//TODO - think if might be accessed with gPad
	PaveText (const Double_t x2, Option_t *option="NDCNB"):
	PaveText (0.16, 0.95, x2, 1, option) {}
};

class myLegend : public TLegend
{
public:
	myLegend (const Double_t x1, const Double_t y1, const Double_t x2, const Double_t y2):
	TLegend(x1, y1, x2, y2)
	{
		this->SetTextFont(22);
		this->SetBorderSize(1); //No shadow
	}
};


template <typename rootObject>
TCanvas& drawWithResiduals(TCanvas& canvas, rootObject object, const TF1& function)
{
	canvas.cd();
	TPad* topPad = new TPad( (canvas.GetName() + std::string("_topPad")).data(), "", 0, 0.2, 1, 1);
	topPad->SetTopMargin(0.03);
	topPad->SetBottomMargin(0.);
	topPad->Draw();
	topPad->cd();
	object.Draw();
	topPad->Update();
	
	//Change to canvas before creating second pad
	canvas.cd();
	
	TPad* bottomPad = new TPad( (canvas.GetName() + std::string("_bottomPad")).data(), "", 0, 0, 1, 0.2);
	topPad->SetBottomMargin(0);
	bottomPad->Draw();
	bottomPad->cd();
	//Draw TGraph
	
	return canvas;
}

double FWHM(const TH1& hist);

inline double FWHM_2355(const TH1& hist) { return FWHM(hist) / 2.355; }
} //namespace