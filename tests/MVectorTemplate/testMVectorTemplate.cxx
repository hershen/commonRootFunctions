#include <iostream>

#include "MVectorTemplate.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TH1D.h"
#include "mathFuncs.h"
#include "TStyle.h"
#include "TRandom3.h"

#define cline std::cout << "line = " << __LINE__ << std::endl;

size_t numPoints = 100;
size_t numPreZeros = 20;
double dx = 2;
double STD = 0.00005;

TRandom3 ran(0);

double pedestalMin = 0;
double pedestalMax = 10;
double ampMin = 0.1;
double ampMax = 10;
double timeShiftMin = -10;
double timeShiftMax = 10;
  
std::vector<double> fillHistWithValues(TH1D* hist, size_t numPoints, TF1 function, double STD)
{
  std::vector<double> vec;
  double value;
  double pedestal = function.GetParameter(2);
  for(unsigned int idx = 1; idx <= numPoints; ++idx)
  {    
    if(idx-1 < numPreZeros) value = pedestal;
    else value = function.Eval( (idx-1-numPreZeros) * dx);
    
    vec.push_back(value);
    
    if(hist)
    {
      hist->SetBinContent(idx,value);
      hist->SetBinError(idx,STD);
    }
  }
  
  return vec;
}
  

void noName()
{

  
//   TF1 realFunc("realFunc","[0]*(0.5*x*x-[1]*x)+[2]",-30,30);
//   TF1 realFunc("realFunc","[0]*TMath::Landau(x - 4 - [1],1,8,0) + [2]",-10,100);
  TF1 realFunc("realFunc","[0]*TMath::GammaDist( (x-[1]<20?20:x - [1]),2,20,10) + [2]",-10,100);
  double pedestal = 0;
  double timeShift = 0;
  double amplitude = 1.;
  
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  
  
  
  std::vector<double> vec = std::vector<double>(numPoints), vec2 = std::vector<double>(numPoints), vec3 = std::vector<double>(numPoints);
  
  TH1D *hist = new TH1D("hist","hist",numPoints,-dx/2.,numPoints*dx-dx/2.);
  TH1D *hist2 = new TH1D("hist2","hist2",numPoints,-dx/2.,numPoints*dx-dx/2.);
  
  vec = fillHistWithValues(hist, numPoints, realFunc, STD);
  
  
  pedestal = 0;
  timeShift = 4.1;
  amplitude = 1.;
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  
  vec2 = fillHistWithValues(hist2, numPoints, realFunc, STD);
  
  pedestal = 0;
  timeShift = -3;
  amplitude = 1.;
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  vec3 = fillHistWithValues(0, numPoints, realFunc, STD);
  
  MVectorTemplate * vt= new MVectorTemplate();
  vt->setDebugLevel(10);
  vt->setDx(dx);
  
  for(int idx = 0 ; idx < 1000; ++idx)
  {
    pedestal = ran.Uniform(pedestalMin,pedestalMax);
    timeShift = ran.Uniform(timeShiftMin,timeShiftMax);
    amplitude = ran.Uniform(ampMin,ampMax);
//     std::cout << "idx = " << idx << ", pedestal = " << pedestal << ", timeShift = " << timeShift << ", amplitude = " << amplitude << std::endl;
    realFunc.SetParameters(amplitude, timeShift, pedestal);
    
    vec = fillHistWithValues(0, numPoints, realFunc, STD);
    vt->addVector(vec,STD);
  }
  TF1 *templateFunc = vt->getTF1();
  
  cline
  
  //Timing resolution
  TH1D *timingHist = new TH1D("timingHist","timingHist",50,-120,20.);
  TH1D *timingErrHist = new TH1D("timingErrHist","timingErrHist",100,1.2,1.7);
  TH1D *pedestalHist = new TH1D("pedestalHist","pedestalHist",100,-0.2,0.2);
  TH1D *ampHist = new TH1D("ampHist","ampHist",50,-15,5.);
  
//   for(int idx = 0 ; idx < 1000; ++idx)
//   {
//     pedestal = ran.Uniform(pedestalMin,pedestalMax);
//     timeShift = ran.Uniform(timeShiftMin,timeShiftMax);
//     amplitude = ran.Uniform(ampMin,ampMax);
//     realFunc.SetParameters(amplitude, timeShift, pedestal);
//     vec = fillHistWithValues(hist, numPoints, realFunc, STD);
//     
//     hist->Fit(templateFunc,"QME");
//     
//     double fittedAmp = templateFunc->GetParameter(0);
//     double fittedPedestal= templateFunc->GetParameter(1);
//     double fittedTime = templateFunc->GetParameter(2);
//     double fittedTimeErr = templateFunc->GetParError(2);
//     
//     std::cout << "idx = " << idx << ", pedestal = " << pedestal << ", timeShift = " << timeShift << ", amplitude = " << amplitude << std::endl;
//     std::cout << "fitted time = " << fittedTime << std::endl;
//     
//     timingHist->Fill(timeShift - fittedTime);
//     timingErrHist->Fill(fittedTimeErr);
//     ampHist->Fill(fittedAmp / amplitude);
//     pedestalHist->Fill(pedestal - fittedPedestal);    
//   }
//   
//   cline
//   
//   TCanvas *timingCan = new TCanvas("timingCan","timingCan",10,10,700,700);
//   timingHist->Draw();
//   timingCan->SaveAs("timingHist.png");
//   
//   TCanvas *timingErrCan = new TCanvas("timingErrCan","timingErrCan",10,10,700,700);
//   timingErrHist->Draw();
//   timingErrCan->SaveAs("timingErrHist.png");
//   
//   TCanvas *ampCan = new TCanvas("ampCan","ampCan",10,10,700,700);
//   ampHist->Draw();
//   ampCan->SaveAs("ampHist.png");
//   
//   TCanvas *pedestalCan = new TCanvas("pedestalCan","pedestalCan",10,10,700,700);
//   pedestalHist->Draw();
//   pedestalCan->SaveAs("pedestalHist.png");
  
  
  TCanvas *can2 = new TCanvas("can2","can2",10,10,700,700);
  can2->cd();
  
  pedestal = ran.Uniform(pedestalMin,pedestalMax);
  timeShift = ran.Uniform(timeShiftMin,timeShiftMax);
  amplitude = ran.Uniform(ampMin,ampMax);
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  vec = fillHistWithValues(hist, numPoints, realFunc, STD);
  
  std::cout << "pedestal = " << pedestal << ", timeShift = " << timeShift << ", amplitude = " << amplitude << std::endl;
  
  hist->Fit(templateFunc,"QM");
  hist->Draw("");
//   hist->SetLineColor(kBlue);
//   hist2->Draw();
//   hist->Draw("SAME");
//   hist->Fit(templateFunc,"M");
//   templateFunc->Draw("");
//   templateFunc->SetNpx(10000);
  can2->SaveAs("can2.png");

  
}



int main(/*int argc, char *argv[]*/)
{
//   loadSingle
  gStyle->SetOptStat("ourme");
  gStyle->SetOptFit(1111);
  noName();
  return 0;
}

void testMVectorTemplate()
{
//   gROOT->ProcessLine(".L /media/sf_PhD/Root/commonRootFunctions/include/MVectorTemplate.h");
//   gROOT->ProcessLine(".L /media/sf_PhD/Root/commonRootFunctions/src/MVectorTemplate.cxx");
  noName();
}
  

