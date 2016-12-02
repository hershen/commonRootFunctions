#include "histFuncs.h"


#include <algorithm>

#include "TRandom3.h"
#include "TMath.h"

#include "TPaveText.h"
#include "TFrame.h"



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
  if( abs(textHeight) > 1e-9 )
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

} //namespace histFuncs
