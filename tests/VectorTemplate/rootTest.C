#include <iostream>

#include "MVectorTemplate.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TH1D.h"
// #include "mathFuncs.h"
// #include "TStyle.h"
#include "TROOT.h"

void rootTest()
{
//   gROOT->ProcessLine(".L /media/sf_PhD/Root/commonRootFunctions/include/MVectorTemplate.h");
//   gROOT->ProcessLine(".L /media/sf_PhD/Root/commonRootFunctions/src/MVectorTemplate.cxx");
  size_t numPoints = 50;
  size_t numPreZeros = 2;
  double dx = 2.;
  
//   TF1 realFunc("realFunc","[0]*(0.5*x*x-[1]*x)+[2]",-30,30);
//   TF1 realFunc("realFunc","[0]*TMath::Landau(x - 4 - [1],1,8,0) + [2]",-10,100);
  TF1 realFunc("realFunc","[0]*TMath::GammaDist( (x-[1]<0?0:x - [1]),2,0,10) + [2]",-10,100);
  double pedestal = 0;
  double timeShift = 0;
  double amplitude = 1.;
  
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  
//   TCanvas can2("can2","can2",10,10,700,700);
  
  std::vector<double> vec = std::vector<double>(numPoints), vec2 = std::vector<double>(numPoints-10);
  
  for(size_t idx = 0; idx < numPoints; ++idx) 
  {
    if(idx < numPreZeros) vec[idx] = pedestal;
    else vec[idx] = realFunc.Eval( (idx-numPreZeros) * dx);
  }
  
  //   pedestal = 10;
  timeShift = -5;
//   amplitude = 5.;
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  for(size_t idx = 0; idx < numPoints-10; ++idx) 
  {
    if(idx < numPreZeros) vec2[idx] = pedestal;
    else vec2[idx] = realFunc.Eval( (idx-numPreZeros) * dx);
    std::cout << "vec2[ " << idx << "] = " << vec2[idx] << std::endl;
  }
  
  
  MVectorTemplate *vt = new MVectorTemplate();
  vt->setDebugLevel(100);
  vt->setDx(dx);
  vt->enablePedestalFit(false);
  vt->enableAmplitudeFit(false);
  vt->addVector(vec,0.0005);
//   vt->addVector(vec2,0.0005);
  TF1 *templateFunc = vt->getTF1();
  
  TH1D *hist = new TH1D("hist","hist",vec2.size(),-dx/2.,vec2.size()*2.-dx/2.);
  for(unsigned int idx = 1; idx <= vec2.size(); ++idx)
  {
    hist->SetBinContent(idx,vec2[idx-1]);
//     fitHist.SetBinError(idx,std);
  }
  TCanvas * c3 = new TCanvas("c3","c3",10,10,700,700);
//   hist->Draw();
  templateFunc->SetNpx(10000);
//   templateFunc->SetParameters(1,0,10);
  
  hist->Fit(templateFunc,"SVM0");
  hist->Draw();
  templateFunc->SetParameters(0.03661,0,-5);
  templateFunc->Draw("SAME");
  
  std::cout << "f(40) = " << templateFunc->Eval(40) << std::endl;
  
  
  TCanvas * c6 = new TCanvas("c6","c6",10,10,700,700);
  
  
  timeShift = 0;
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  realFunc.DrawClone();
  
  timeShift = -5;
  realFunc.SetParameters(amplitude, timeShift, pedestal);
  realFunc.DrawClone("SAME");
  
}