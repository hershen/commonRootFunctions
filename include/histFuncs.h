#pragma once

#include <iostream>

#include "boost/format.hpp"


#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"

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

/**
* Extract a root object rootObjName from a root file file. The file is assumed to be not zombie. If rootObjName doesn't exist in file, print error.
* rootClass should be a pointer (unless you figure out how to use it otherwise).
*/
void binGraph(const TGraph& graph, TH1& hist)
{
	
}

} //namespace