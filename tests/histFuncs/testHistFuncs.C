#include "histFuncs.h"

void testbinGraph()
{
	
	TGraph* graph = new TGraph();
	TRandom3 rand(0);
	
	double maxX = 3.0;
	TF1 f("f","TMath::Exp(x)",0,maxX);
	for(size_t i = 0; i < 1e4; ++i) 
	{
		double x = maxX/1e4 * double(i);
		double y = f.Eval(x);
		graph->SetPoint(i, x, y + rand.Gaus(0,maxX/1e4 * double(i)) );
		
	}
	
	TH1D* hist = new TH1D("hist","",100,0,maxX);
	myFuncs::binGraph(*graph, *hist);
	
	TCanvas* c = new TCanvas("c","",10,10,1000,750);
	hist->Draw("h");
	
	TCanvas* cg = new TCanvas("cg","",10,10,1000,750);
	graph->Draw("AP*");
	
	for(int i = 0; i <= hist->GetNbinsX(); ++i) std::cout << "bin " << i << " has " << hist->GetBinContent(i) << " entries and error " << hist->GetBinError(i) << std::endl;
	
	std::cout << "sum W = " << hist->GetSumOfWeights() << std::endl;
	
	std::cout << std::endl;
	
// 	for(int i = 0; i < graph.GetN() ; ++i) std::cout << "point " << i << " =(" << graph.GetX()[i] << "," << graph.GetY()[i] << ")" << std::endl;
	
}

void testHistFuncs()
{
	testbinGraph();
	
	
}