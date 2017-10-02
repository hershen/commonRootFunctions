#include "histFuncs.h"


#include <algorithm>

#include "TRandom3.h"
#include "TMath.h"

#include "TPaveText.h"
#include "TFrame.h"
#include "TF1.h"


bool compareHistMaximum(const TH1* hist1, const TH1* hist2) { return hist1->GetBinContent(hist1->GetMaximumBin()) < hist2->GetBinContent(hist2->GetMaximumBin()); }

namespace myFuncs
{       
void drawHistograms_highestFirst(const std::vector<TH1*> &histVector)
{
  
  if(histVector.size() < 1) return; 
  
  //Find highest histogram
  auto highestHist = std::max_element(histVector.begin(), histVector.end(), compareHistMaximum );
  if( highestHist == histVector.end()) return; //histVector is empty
  
  //Draw highest histogram
  (*highestHist)->Draw();
  
  //Draw rest of histograms
  for(auto hist : histVector)
  {
    if(hist == *highestHist) continue;
    hist->Draw("SAME");
  }
}
  
  
TCanvas* newCanvas(std::string canvasName) 
{
  if(canvasName != "") return new TCanvas(canvasName.data(),canvasName.data(),10,10,1920,985);
    
  TRandom3 rand3(0);
  double randNum = rand3.Uniform();
  canvasName = "canvas" + std::to_string( int(randNum * 1000) );
  return new TCanvas(canvasName.data(),"",10,10,1920,985);
}

TCanvas* drawNewCanvas(TH1* hist)
{
  std::string canvasName = std::string(hist->GetName()) + "Canvas";
  TCanvas* canvas = new TCanvas(canvasName.data(),canvasName.data(),10,10,1920,985);
  hist->Draw();
  return canvas;
}

TCanvas* drawNewCanvas(TGraph* graph)
{
  std::string canvasName = std::string(graph->GetName()) + "Canvas";
  TCanvas* canvas = new TCanvas(canvasName.data(),canvasName.data(),10,10,1920,985);
  graph->Draw("AP");
  return canvas;
}

void putLabelAboveBin(const TH1* hist, const size_t bin, const std::string& text, const double textHeight, const double yOffset)
{
  gPad->Update();
  double binWidth = hist->GetBinWidth(bin);
  double epsilon = binWidth*0.05;
  double xMin = hist->GetBinCenter(bin) - binWidth/2. + epsilon;
  double xMax = hist->GetBinCenter(bin) + binWidth/2. - epsilon;
  
  double yMin = hist->GetBinContent(bin) + yOffset;
  if( gPad->GetLogy() && yMin > 0) yMin = TMath::Log10(yMin - yOffset) + yOffset;

  double yMax = 0.0;
  if( std::abs(textHeight) > 1e-9 )
    yMax = yMin + yOffset + textHeight;
  else //height = top border - max bin height
  {
    double gPadHeight = gPad->GetFrame()->GetY2();
    double histMax = hist->GetMaximum();
    if( gPad->GetLogy() ) histMax = TMath::Log10( histMax );
    double height = gPadHeight - histMax;
    yMax = yMin + yOffset + height - (height * 0.05);
//     std::cout << "yMax = " << yMax << ", histMax = " << histMax << ", gPadHeight = " << gPadHeight << ", height = " << height << std::endl;
  }

  if( gPad->GetLogy() )
  {
    yMin = std::pow(10,yMin);
    yMax = std::pow(10,yMax);
  }
  
  TPaveText* paveText = new TPaveText(xMin,yMin,xMax,yMax,"NB");
  paveText->AddText(text.data());
  paveText->SetFillColorAlpha(kWhite,0); //0 - fully transparent. Not working for some reason...
  paveText->Draw();
  
}

double integrateTGraph(const TGraph& graph)
{
	//Graph needs at least 2 points.
	if(graph.GetN() < 2) return 0.0;
	
	//Create copy so we can sort
	TGraph sortedGraph = graph;
	
	//Sort points
	sortedGraph.Sort();
	
	//Calculate integral, assuming linear interpolaion between points
	//Trapezoid area = (a+b)h/2
	double sum = 0.0;
	double xi = 0;
	double yi = 0;
	double xiP1 = 0;
	double yiP1 = 0;
	
	//Prevent warning of unsigned integer expression because TGraph::GetN returns int isntead of unsigned int
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare"  
	for(size_t iPoint = 0; iPoint < sortedGraph.GetN() - 1; ++iPoint)
	#pragma GCC diagnostic pop
	{
		
		sortedGraph.GetPoint(iPoint, xi, yi);
		sortedGraph.GetPoint(iPoint + 1, xiP1, yiP1);
		double h = xiP1 - xi;
		
		//guard
		if(h <= 0.0) return 0.0;
		
		sum += (yiP1 + yi) * h / 2.0;
	}
	
	return sum;
}

void normalize(TGraph& graph, const double area)
{
	if(area == 0.0) return;
	
	//get area
	double integral = integrateTGraph(graph);
	
	if(integral == 0.0)
	{
		std::cout << "histFuncs::normalize: Integral of graph == 0. Returning without scaling!!!" << std::endl;
		return;
	}
	
	//Normalize each point
	
	const double scaleFactor = area/integral;
	
	//Prevent warning of unsigned integer expression because TGraph::GetN returns int isntead of unsigned int
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare" 
	for(size_t iPoint = 0; iPoint < graph.GetN(); ++iPoint)
	#pragma GCC diagnostic pop
	{
		double x = 0;
	  double y = 0;
		graph.GetPoint(iPoint, x, y);
		graph.SetPoint(iPoint, x, y * scaleFactor);
	}	
}

void normalize(TH1D& hist, const double area, const std::string& options) 
{ 
	if(area == 0.0) return;
	
	const double integral = hist.Integral(options.data());
	if(integral == 0.0) 
	{
		std::cout << "histFuncs::normalize: Integral of histogram == 0. Returning without scaling!!!" << std::endl;
		return;
	}

	hist.Scale(area/integral); 	
} 

void binGraph(const TGraph& graph, TH1& hist, const bool errorOnMean)
{
	//Histogram to store number of points in each bin.
	TH1* numEntriesHist = static_cast<TH1*>( hist.Clone() );
 	numEntriesHist->FillN(graph.GetN(), graph.GetX(), nullptr);

	//Fill histogram with points from graph. 
	//Each bin contains the sum of all the y values in that bin.
	hist.FillN(graph.GetN(), graph.GetX(), graph.GetY() );

	//Prevent warning of unsigned integer expression because TGraph::GetN returns int isntead of unsigned int
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare" 
	for(size_t iBin = 0; iBin <= hist.GetNbinsX(); ++iBin)
	#pragma GCC diagnostic pop
	{
		double content = hist.GetBinContent(iBin); // sum of y values in the bin.
		double error = hist.GetBinError(iBin);     // sum of y^2
		long entries = numEntriesHist->GetBinContent(iBin) ; //num of points in the bin
		
		if(entries > 0) 
		{
			hist.SetBinContent(iBin, content / entries ); //Set content to average
			error = std::sqrt(std::pow(error, 2) / entries - std::pow(content / entries,2) ); //Calculate std in bin
			
			if(errorOnMean) error /= std::sqrt( entries ); //Set error to error mean (std/sqrt(n))
		
			hist.SetBinError(iBin, error); //Set error.
		}
		
	}
	delete numEntriesHist;
}

TGraphErrors histToGraph(const TH1& hist)
{
	TGraphErrors graph;
	
	//Prevent warning of unsigned integer expression because TGraph::GetN returns int isntead of unsigned int
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare" 
	for(size_t iBin = 1; iBin <= hist.GetNbinsX(); ++iBin)
	#pragma GCC diagnostic pop
	{
		size_t iPoint = graph.GetN();
		graph.SetPoint(iPoint, hist.GetBinCenter(iBin), hist.GetBinContent(iBin) );
		graph.SetPointError( iPoint, hist.GetBinWidth(iBin), hist.GetBinError(iBin) );
	}	
	return graph;
	
}

double FWHM(const TH1& hist) {
	const double halfMaximum = hist.GetMaximum() / 2.;
  const auto firstBin = hist.FindFirstBinAbove(halfMaximum);
  const auto lastBin = hist.FindLastBinAbove(halfMaximum);
  return hist.GetBinCenter(lastBin) - hist.GetBinCenter(firstBin);
}

void prepareCanvasResiduals(TCanvas& canvas) {
	const double topBottomSpace = 0.045;
	
	canvas.cd();
	TPad* topPad = new TPad( (canvas.GetName() + std::string("_topPad")).data(), "", 0, 0.25, 1, 1);
	topPad->SetNumber(1);
	topPad->SetTopMargin(0.03);
	topPad->SetBottomMargin(0.01); //For some reason 0.0 doesn't look good
	topPad->Draw();
	topPad->cd();
	topPad->Update();
	
	//Change to canvas before creating second pad
	canvas.cd();
	
	TPad* bottomPad = new TPad( (canvas.GetName() + std::string("_bottomPad")).data(), "", 0, 0, 1, 0.25);
	bottomPad->SetNumber(2);
	bottomPad->SetTopMargin(topBottomSpace);
	bottomPad->SetBottomMargin(0.25); //So that title isn't cut out
	bottomPad->SetGridy();
	bottomPad->Draw();
	canvas.cd();
}

double getHistOverlapChi2(TH1D hist0, TH1D hist1, const double minVal, const double maxVal) {
	
	const int hist0MinBin = hist0.FindBin(minVal);
	const int hist1MinBin = hist1.FindBin(minVal);
	const int hist0MaxBin = hist0.FindBin(maxVal);
	const int hist1MaxBin = hist1.FindBin(maxVal);
	
	const int hist0Range = hist0MaxBin - hist0MinBin;
	const int hist1Range = hist1MaxBin - hist1MinBin;
	
	
	//Check that hists have same number of bins.
	if(hist0Range != hist1Range) {
		std::cout << "histFuncs::getHistOverlapChi2 : Currently supports only histograms with same number of bins in the range (" << minVal << "," << maxVal << ".\n Hist0 has " << hist0Range << " bins and hist 1 has " << hist1Range << " bins." << std::endl;
		return 0.0;
	}
	
	double chi2 = 0.0;
	for(int iHist0Bin = hist0MinBin, iHist1Bin = hist1MinBin; iHist0Bin <= hist0MaxBin; ++iHist0Bin, ++iHist1Bin) {
		const double diff = (hist0.GetBinContent(iHist0Bin) - hist1.GetBinContent(iHist1Bin));
		const double error2 = std::pow(hist0.GetBinError(iHist0Bin),2) + std::pow(hist1.GetBinError(iHist1Bin),2);
		chi2 += std::pow(diff,2)/error2;
	}
	
	return chi2;
}


void setNovosibirskParams(TF1& novosibirskTF1, const TH1& hist) {
	const double norm = hist.GetMaximum();
	const double peak = hist.GetBinCenter( hist.GetMaximumBin() );
	const double width = myFuncs::FWHM(hist);
	
	//Set parameters
	novosibirskTF1.SetParameters(norm, peak, width, 0.0);
	
	//Set width limit
	novosibirskTF1.SetParLimits(2, 0.0, width * 10);
}

} //namespace histFuncs
