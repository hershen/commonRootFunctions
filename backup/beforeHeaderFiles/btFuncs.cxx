#pragma once

#ifndef btFunctions_cxx
#define btFunctions_cxx

#include <iostream>
#include "TString.h"

//Enumurate TOF channels in V1730
enum TOFchannels {eTOFupstream = 4, eTOFdownstream0 = 12, eTOFdownstream1 = 13, eTOFscint = 6};



namespace btFuncs
{
  TString rootFilesDir = "/media/sf_PhD/Testbeam2015/convertedMidasFiles/";
  TString fileTemplate = "convertedbeamtest_00000%d_%s.mid.gz.root";
  TString treeName = "TimingAndWaveform";
  TString crystalChannelBranchTemplate = "ch%d";
  TString TOFtimeBranchTemplate = "TOFchannel%d";
  TString TOFdofBranchTemplate = "chi2DOFchannel%d";
  
  TString TOFupstreamChBranch = TString::Format( TOFtimeBranchTemplate.Data(), eTOFupstream);
  TString TOFdownstream0ChBranch = TString::Format( TOFtimeBranchTemplate.Data(), eTOFdownstream0);
  TString TOFdownstream1ChBranch = TString::Format( TOFtimeBranchTemplate.Data(), eTOFdownstream1);
  TString TOFscintChBranch = TString::Format( TOFtimeBranchTemplate.Data(), eTOFscint);
  
  const double sampleTime = 2.;
  
  const double ukrainianADC2MeV = 24.66/1056.;
  const double BelleADC2MeV = 24.66/79.51;
  //--------------------------------------------------------------------------------------------
  //TOF_3way
  //********************************************************************************************
  //Simple calculation of TOF using only the 3 TOF microchannels.
  //Downstream time - Upstream time
  inline double TOF_3way(double TOFupStream, double TOFdownStream0, double TOFdownStream1) {return (TOFdownStream0 + TOFdownStream1)/2. - TOFupStream;}
  
  
  template <typename Type>
  double getSimpleAmp(std::vector<Type> *waveform, unsigned int idxForPedestal)
  {
    Type sum = 0;
    Type maxVal = 0;
    for(size_t idx = 0; idx < waveform->size(); ++idx) 
    {
      if( idx < idxForPedestal) sum += (*waveform)[idx];
      if( (*waveform)[idx] > maxVal) maxVal = (*waveform)[idx];
    }
    double amp = double(maxVal) - double(sum)/double(idxForPedestal);
    return amp;
  }
    
  
  //--------------------------------------------------------------------------------------------
  //analytical_RC_CRn
  //********************************************************************************************
  //Analytical response of an RC-(CR)^n filter to a step function starting at time t = 0.
  //Function is 1/n! * (t/tau)^n * exp(-t/tau)
  //tau is the shaping time.
  //Function is shifted in time so that it starts at time startTime.
  //amplitude controls the amplitude.
  //max time is an approximate time for the function to decrease below exp(-expCoefficient)
  //--------------------------------------------------------------------------------------------
//   TF1 analytical_RC_CRn(int n, double tau = 500., double amplitude = 1., double startTime = 0., double maxTime = 0.)
//   {
//     double expCoefficient = 5.;
//     if( abs(maxTime) < 1e-9) maxTime = -tau*(TMath::Log( (amplitude<0?(-amplitude):amplitude) ) + expCoefficient-n*TMath::Log(tau)+TMath::Log(TMath::Factorial(n)) );
//     TF1 analytical("analytical","(x-[3])<0?0:[2]*pow((x-[3])/[1],[0])*TMath::Exp(-((x-[3])/[1]))/TMath::Factorial([0])",0,maxTime);
//     analytical.FixParameter(0,n);
//     analytical.FixParameter(1,tau);
//     analytical.FixParameter(2,amplitude);
//     analytical.SetParameter(3,startTime);
//     
//     return analytical;
//   }
  
} //namespace btFunctions

#endif