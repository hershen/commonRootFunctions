#include "MVectorTemplate.h"

#include <iostream>

#include "TF1.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include "mathFuncs.h"
#include "fileFuncs.h"

#define cline std::cout << "line = " << __LINE__ << std::endl;

//For TMiuit2
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"

MVectorTemplate::MVectorTemplate():
  rNumAveragedFuncs(0), 
  rDx(0.), 
  rTemplateValues(std::vector<double>()),
  bAmplitudeFitEnabled(true),
  bPedestalFitEnabled(true),
  rTF1(TF1()),
  rPeakIdx(0),
  rDebugLevel(0),
  rUseAmplitudeLimits(false),
  rMinAmplitudeLimit(0),    
  rMaxAmplitudeLimit(-1),
  rMinChi2Limit(-1),
  rMaxChi2Limit(DBL_MAX),
  rUseXshiftLimits(false),
  rMinXshiftLimit(0),
  rMaxXshiftLimit(0), 
  rUsePedestalLimits(false),
  rMinPedestalLimit(0),
  rMaxPedestalLimit(0),
  rTreeName("vectorTemplate"),
  rTemplateValuesBranchName("rTemplateValues"),
  rDxBranchName("rDx"),
  rNumAveragedFuncsBranchName("rNumAveragedFuncs"),
  rPeakIdxBranchName("rPeakIdx"),
  rDoubleNumbersEqualThershold(1e-20)
{

}

MVectorTemplate::MVectorTemplate(const std::vector<double>& newVec, const double newDx):
  MVectorTemplate()
{
  setDx(newDx);
  addVector(newVec,0);
}

void MVectorTemplate::setTF1ParNames()
{
  rTF1.SetParName(0,"Amplitude");
  rTF1.SetParName(1,"Pedestal");
  rTF1.SetParName(2,"TimeShift");
}


void MVectorTemplate::resetTemplateValues()
{
  rTemplateValues = std::vector<double>(); 
  rNumAveragedFuncs = 0;
  rTF1 = TF1();
}

void MVectorTemplate::resetTemplateRange()
{
  if(rTF1.IsValid()) rTF1.SetRange(0,rDx*(rTemplateValues.size()-1) );
}

void MVectorTemplate::setDx(const double newDx)
{
  rDx = newDx;
  resetTemplateRange();
}

int MVectorTemplate::addVector(const std::vector<double>& newVector, const double std)
{
  //TODO implement with TMinuit
  //TODO average using a weighted average
  // ------------------------------------------------------------------
  //   //Sanity checks
  // ------------------------------------------------------------------
  if(newVector.empty() )
  {
    std::cerr << "MVectorTemplate::addVector : Error - trying to add empty vector. Not adding. " << std::endl;
    return -1;    
  }
  
  if(std < 0 )
  {
    std::cerr << "MVectorTemplate::addVector : Error - negative standard deviation (" << std << "). Not adding. " << std::endl;
    return -2;    
  }
  
  
  // ------------------------------------------------------------------
  //Set first vector as template, normalize, etc.
  // ------------------------------------------------------------------
  if(rNumAveragedFuncs == 0) //First vector added
  {

    // ------------------------------------------------------------------
    //Take pedestal as first entry
    // ------------------------------------------------------------------
    double pedestal = newVector[0];

    // ------------------------------------------------------------------
    //Search for maximum (and what idx it is at - for later guessing the time shift) of newVector so peak can be normalized
    // ------------------------------------------------------------------
    double maximum = -DBL_MAX;
    for(size_t idx = 0; idx < newVector.size(); ++idx)
    {
      if(newVector[idx] > maximum)
      {
	maximum = newVector[idx];
	rPeakIdx = idx;
      }
    }
    double maxMinusPedestal = maximum - pedestal;
    // ------------------------------------------------------------------
    //Populate rTemplateValues, normalizing and subtracting pedestal
    // ------------------------------------------------------------------
    rTemplateValues = std::vector<double>(newVector.size());
    for(size_t idx = 0; idx < newVector.size(); ++idx)
    {
      rTemplateValues[idx] = (newVector[idx] - pedestal)/maxMinusPedestal;
    }
    
    if (rDebugLevel > 20)
    {
      std::cout << "MVectorTemplate::addVector : pedestal = " << pedestal << ", maximum = " << maximum << std::endl;
      for(size_t iTemplateIdx = 0; iTemplateIdx < rTemplateValues.size(); ++iTemplateIdx)
	std::cout << "MVectorTemplate::addVector : temp[" << iTemplateIdx << "] = " << rTemplateValues[iTemplateIdx] << std::endl;
    }
    
    // ------------------------------------------------------------------
    //Create TF1
    // ------------------------------------------------------------------
    rTF1 = TF1("templateTF1",this,&MVectorTemplate::TF1Eval,0,rDx*(rTemplateValues.size()-1),3,"MVectorTemplate","TF1Eval");
    setTF1ParNames();
  }//First vector added
  
  else //from second vector
  {
    
    //Standard deviation vector
    std::vector<double> stdVector(newVector.size(),std);
    
    //Histrogram of with entries of newVector - used to fit the template
    //Center of each bin corresponds to correct x value (that's why there's a -rDx/2. term)
    TH1D fitHist("fitHist","fitHist", newVector.size(), -rDx/2. , newVector.size() * rDx - rDx/2.);
    for(unsigned int idx = 1; idx <= newVector.size(); ++idx)
    {
      fitHist.SetBinContent(idx,newVector[idx-1]);
      fitHist.SetBinError(idx,std);
    }
    
    // ------------------------------------------------------------------
    //Initial guesses for fit
    // ------------------------------------------------------------------
    //Amplitude initial guess is the value of newVector where the template maximum is. This will not work for fast varying functions, or if the xShift is large
    double amplitudeGuess = newVector[rPeakIdx] / rTemplateValues[rPeakIdx]; 
    
    //Make sure amplitude guess is not out of limits
    if(rUseAmplitudeLimits) 
    {
      amplitudeGuess = std::min(amplitudeGuess, rMaxAmplitudeLimit);
      amplitudeGuess = std::max(amplitudeGuess, rMinAmplitudeLimit);
    }
    double pedestalGuess = newVector[0];
    double xShiftGuess = 0;
    
    setTF1Parameters(amplitudeGuess, pedestalGuess, xShiftGuess);
    
//     std::cout << "amplitudeGuess = " << amplitudeGuess << ", pedestalGuess = " << pedestalGuess << ", xShiftGuess = " <<  xShiftGuess << std::endl;
    
    
    // ------------------------------------------------------------------
    //Fit!
    // ------------------------------------------------------------------
    TFitResultPtr fitResult = fitHist.Fit(&rTF1,"QSN");
    //"Q" Quiet mode (minimum printing)
    //"M" More. Improve fit results. It uses the IMPROVE command of TMinuit (see TMinuit::mnimpr). This algorithm attempts to improve the found local minimum by searching for a better one.
    //"N" Do not store the graphics function, do not draw
    //"S" The result of the fit is returned in the TFitResultPtr (see below Access to the Fit Result)
    // ------------------------------------------------------------------
    //Fit result checks
    // ------------------------------------------------------------------
    
    //TODO implement!
    //How to check if first vector added is corrupt (i.e. flat). Maybe have a cut on maximum - pedestal.
    //Amplitude, pedestal, xShift limits/cuts
    //Chi2 cut?
    
    
    int fitStatus = fitResult;
    if( fitStatus )
    {
      if(fitStatus != 4000)
      {
	if (rDebugLevel > 0) std::cerr << "MVectorTemplate::addVector : fit failed, fit status = " << fitStatus << ". Not adding vector!" << std::endl;
	return -10;
      }
      else if(rDebugLevel > 0) std::cout << "MVectorTemplate::addVector : MINOS falied, fit status = " << fitStatus << ". Vector MIGHT still be added!!!" << std::endl;
    }
      
    // ------------------------------------------------------------------
    //Positive fittedXshift means newVector is forward in time relative to rTemplateValues, i.e. newVector[idx] = rTemplateValues[idx-fittedXshift/rDx]
    // ------------------------------------------------------------------
    double fittedAmplitude = fitResult->Parameter(0);
    double fittedPedestal = fitResult->Parameter(1);
    double fittedXshift = fitResult->Parameter(2);
    
    double chi2_NDF = fitResult->Chi2();
       
    if(rDebugLevel > 10) std::cout << "MVectorTemplate::addVector : fittedAmplitude = " << fittedAmplitude << ", fittedPedestal = " << fittedPedestal << ", fittedXshift = " << fittedXshift << ", chi2 / NDF = " << chi2_NDF << std::endl;
    
    
    if( (abs(fittedAmplitude - rMinAmplitudeLimit) < rDoubleNumbersEqualThershold || abs(fittedAmplitude - rMaxAmplitudeLimit) < rDoubleNumbersEqualThershold)  && rUseAmplitudeLimits)
    {
      if(rDebugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached amplitude limit. fittedAmplitude = " << fittedAmplitude << ". limits are (" << rMinAmplitudeLimit << "," << rMaxAmplitudeLimit << "). Not adding vector!" << std::endl;
      return -11;
    }
    
    if( (abs(fittedXshift - rMinXshiftLimit) < rDoubleNumbersEqualThershold || abs(fittedXshift - rMaxXshiftLimit) < rDoubleNumbersEqualThershold) && rUseXshiftLimits)
    {
      if(rDebugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached xShift limit. fittedXshift = " << fittedXshift << ". limits are (" << rMinXshiftLimit << "," << rMaxXshiftLimit << "). Not adding vector!" << std::endl;
      return -12;
    }
    
    if( (abs(fittedPedestal - rMinPedestalLimit) < rDoubleNumbersEqualThershold || abs(fittedPedestal - rMaxPedestalLimit) < rDoubleNumbersEqualThershold) && rUsePedestalLimits)
    {
      if(rDebugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached pedestal limit. fittedPedestal = " << fittedPedestal << ". limits are (" << rMinPedestalLimit << "," << rMaxPedestalLimit << "). Not adding vector!" << std::endl;
      return -13;
    }
    
    if( chi2_NDF < rMinChi2Limit || chi2_NDF > rMaxChi2Limit ) 
    {
      if(rDebugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, Chi2/NDF of fit = " << chi2_NDF << " is too large. limits are (" << rMinChi2Limit << "," << rMaxChi2Limit << "). Not adding vector!" << std::endl;
      return -30; 
    }
    
//     if(chi2_NDF > 2) std::cout << "chi2 / NDF = " << chi2_NDF << std::endl;

    // ------------------------------------------------------------------
    //Average template with new waveform
    //x corresponding to newVector[0] is: - fittedXshift (with minus!!)
    //x corresponding to newVector[size] is size*rDx - fittedXshift
    // ------------------------------------------------------------------
    
    //Find which are the first and last indexes in the template which needs to be averaged.
    //firstTemplateIdx,lastTemplatedIdx will be averaged with the interpolated value from newVector
    size_t firstTemplateIdx = 0, lastTemplateIdx = 0;
    if(fittedXshift > 0)
    {
      //template:                      Idx0---------Idx1---------Idx2          Idx(size-4)-------------Idx(size-3)-------------Idx(size-2)--
      //newVector: Idx0---------Idx1---------Idx2---------Idx3     Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
      //In this example both vectors are same length and fittedXshift = 1.5rDx
      //firstTemplateIdx = 0, lastTemplateIdx = newVector.size()-3
      firstTemplateIdx = 0; 
      lastTemplateIdx = newVector.size() - 1 - size_t(fittedXshift / rDx) - 1; 
    }
    else
    {
      //template:  Idx0--------Idx1--------Idx2          Idx(size-2)-------------Idx(size-1)
      //newVector:                   Idx0--------Idx1                Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
      //In this example both vectors are same length and fittedXshift = 1.5rDx
      //firstTemplateIdx = 2, lastTemplateIdx = rTemplateValues.size()-1
      firstTemplateIdx = size_t(-fittedXshift / rDx) + 1;
      lastTemplateIdx = rTemplateValues.size() - 1; //last element
    }

    // ------------------------------------------------------------------
    //Average newVector with template
    // ------------------------------------------------------------------
    
    //Create template of newVector
    MVectorTemplate newVectorTemplate(newVector, rDx);
    TF1 * newVectorTF1 = newVectorTemplate.getTF1();
    
    double originalNewVectorMax = newVector[ newVectorTemplate.getPeakIdx()];
    double newVectorPedestal = (newVector[0] - fittedPedestal)/(originalNewVectorMax-newVector[0]);;
    
    double newVectorAmplitude = (originalNewVectorMax - newVector[0]) / fittedAmplitude;
    
    newVectorTF1->SetParameters(newVectorAmplitude, newVectorPedestal, -fittedXshift);
    
    //Average
    double peakVal = -DBL_MAX;
    for(size_t iTemplateIdx = firstTemplateIdx; iTemplateIdx <= lastTemplateIdx; ++iTemplateIdx)
    {
      double newVectorValue = newVectorTF1->Eval(iTemplateIdx * rDx);
      
      //average vectors
      rTemplateValues[iTemplateIdx] = (rTemplateValues[iTemplateIdx]*double(rNumAveragedFuncs) + newVectorValue)/double(rNumAveragedFuncs+1);
      
      //Store maximum value (used to normalize later) and its index
      if(rTemplateValues[iTemplateIdx] > peakVal)
      {
	peakVal = rTemplateValues[iTemplateIdx];
	rPeakIdx = iTemplateIdx;
      }
      
      if (rDebugLevel > 20)
	std::cout << "MVectorTemplate::addVector : temp[" << iTemplateIdx << "] = " << rTemplateValues[iTemplateIdx] << std::endl;
    }
    
    // ------------------------------------------------------------------
    //Clip ends of template
    //x corresponding to newVector[0] is - fittedXshift (with minus!!)
    //x corresponding to newVector[size] is size*rDx - fittedXshift
    // ------------------------------------------------------------------
    if(fittedXshift > 0)
    {
      rTemplateValues.resize( lastTemplateIdx + 1); //lastTemplateIdx + 1 == size of new template vector;
      resetTemplateRange();
    }
    else
    {
      if(firstTemplateIdx> 0 ) rTemplateValues.erase( rTemplateValues.begin() , rTemplateValues.begin() + firstTemplateIdx - 1 );
      resetTemplateRange();
    }
    
    // ------------------------------------------------------------------
    //re-Normalize vector - keep peak == 1
    //Make sure first entry is zero
    // ------------------------------------------------------------------
    double pedestal = rTemplateValues[0];
    for(size_t iTemplateIdx = 0; iTemplateIdx < rTemplateValues.size(); ++iTemplateIdx)
    {
      rTemplateValues[iTemplateIdx] = (rTemplateValues[iTemplateIdx] - pedestal )/ (peakVal-pedestal);
    }
    
  }//add, from second vector
  
  ++rNumAveragedFuncs;
  rTF1.SetParameters(1,0,0);
  return 0;
}

void MVectorTemplate::enableAmplitudeFit(const bool &enableAmplitudeFit)
{
  bAmplitudeFitEnabled = enableAmplitudeFit;
  setTF1Parameters();
}

// double MVectorTemplate::eval(const double &x, const double &amplitude, const double &pedestal, const double &xShift) const
// {
//   if( !rTF1.IsValid())
//   {
//     std::cerr << "MVectorTemplate::eval : TF1 not valid" << std::endl;
//     return 0.;
//   }
//   
//   double originalAmplitude = rTF1.GetParameter(0);
//   double originalPedestal = rTF1.GetParameter(1);
//   double originalXshift = rTF1.GetParameter(2);
//   
//   rTF1.SetParameters(amplitude,pedestal,xShift);
//   double returnVal = rTF1.Eval(x);
//   rTF1.SetParameters(originalAmplitude,originalPedestal,originalXshift);
//   
//   return returnVal;
//   
// }

void MVectorTemplate::enablePedestalFit(const bool &enablePedestalFit) 
{
  bPedestalFitEnabled = enablePedestalFit;
  setTF1Parameters();
}

double MVectorTemplate::TF1Eval(double *var, double *params) 
{
//   double amplitude = params[0];
//   double pedistal = params[1];
//   double xShift = params[2];
  double effectiveX = var[0] - params[2]; 
  
//         cout << "amplitude = " << amplitude << ", pedistal = " << pedistal << ", xShift = " << xShift << ", x = " << var[0] << ", effectiveX = " << effectiveX << ". ";
  
  if(effectiveX <= 0) return params[1] + params[0]*rTemplateValues[0];
  if(effectiveX >= rDx * double(rTemplateValues.size() - 1)) return params[1] + params[0]*rTemplateValues.back();
  
  int idx = int(effectiveX / rDx);

  //Linear interpolation
  double x1 = double(idx)*rDx, x2 = double(idx + 1) * rDx;
//     cout << "vec[idx + 1] = " << vec[idx + 1];
  double y1 = rTemplateValues[idx], y2 = rTemplateValues[idx + 1];
  
//     std::cout << "idx = " << idx << ", x1 = " << x1 << ", y1 = " << y1 << ". x2 = " << x2 << ", y2 = " << y2 << std::endl;
//   double effY = (y2 - y1)/(x2-x1)*(effectiveX-x1) + y1;
  
  double returnVal = params[1] + params[0]*((y2 - y1)/(x2-x1)*(effectiveX-x1) + y1);
  
//     std::cout << "returnVal = " << returnVal << std::endl;
  return returnVal;    
}

void MVectorTemplate::setTF1Parameters(const double &amplitude, const double &pedestal, const double &xShift)
{
  if( ! rTF1.IsValid() ) return;
  
  if(bAmplitudeFitEnabled) 
  {
    rTF1.SetParameter(0, amplitude);
    rTF1.SetParLimits(0, rMinAmplitudeLimit,rMaxAmplitudeLimit);
  }
  else rTF1.FixParameter(0, amplitude);
  
  if(bPedestalFitEnabled)
  {
    rTF1.SetParameter(1, pedestal);
    rTF1.SetParLimits(1, rMinPedestalLimit,rMaxPedestalLimit);
  }
  else rTF1.FixParameter(1, pedestal);
  
  rTF1.SetParameter(2, xShift);
  rTF1.SetParLimits(2,rMinXshiftLimit, rMaxXshiftLimit);
}

void MVectorTemplate::setTF1Parameters()
{
  if( ! rTF1.IsValid() ) return;
  setTF1Parameters(rTF1.GetParameter(0), rTF1.GetParameter(1), rTF1.GetParameter(2));
}

void MVectorTemplate::setAmplitudeLimits(const double &newMinAmplitudeLimit, const double &newMaxAmplitudeLimit)
{
  if( abs(newMinAmplitudeLimit - newMaxAmplitudeLimit) > rDoubleNumbersEqualThershold) rUseAmplitudeLimits = true;
  else rUseAmplitudeLimits = false;
  rMinAmplitudeLimit = newMinAmplitudeLimit; 
  rMaxAmplitudeLimit = newMaxAmplitudeLimit;
  setTF1Parameters();
}

void MVectorTemplate::setXshiftLimits(const double &newMinXshiftLimit, const double &newMaxXshiftLimit)
{
  if( abs(newMinXshiftLimit - newMaxXshiftLimit) > rDoubleNumbersEqualThershold) rUseXshiftLimits = true;
  else rUseXshiftLimits = false;
  rMinXshiftLimit = newMinXshiftLimit; 
  rMaxXshiftLimit = newMaxXshiftLimit;
  setTF1Parameters();
}

void MVectorTemplate::setPedestalLimits(const double &newMinPedestalLimit, const double &newMaxPedestalLimit)
{
  if( abs(newMinPedestalLimit - newMaxPedestalLimit) > rDoubleNumbersEqualThershold) rUsePedestalLimits = true;
  else rUsePedestalLimits = false;
  rMinPedestalLimit = newMinPedestalLimit; 
  rMaxPedestalLimit = newMaxPedestalLimit;
  setTF1Parameters();  
}

int MVectorTemplate::saveTemplateToTFile(const std::string &fullFileName, const std::string &treeDescription)
{
  
  // ------------------------------------------------------------------
  //Create TFile
  // ------------------------------------------------------------------
  TFile outputFile(fullFileName.data(),"RECREATE");
  
  if( outputFile.IsZombie() )
  {
    std::cerr << "MVectorTemplate::saveTemplateToTfile : could not open file " << fullFileName.data() << ". Aborting" << std::endl;
    return -1;
  }
  
  // ------------------------------------------------------------------
  //Create TTree
  // ------------------------------------------------------------------
  TTree tree(rTreeName.data(), treeDescription.data());
  
  // ------------------------------------------------------------------
  //Create pointers to variables that will be saved
  // ------------------------------------------------------------------
  std::vector<double> * templateValuesPointer = &rTemplateValues;
  
  // ------------------------------------------------------------------
  //Create branches
  // ------------------------------------------------------------------
  tree.Branch(rTemplateValuesBranchName.data(), "std::vector<double>", &templateValuesPointer);
  tree.Branch(rDxBranchName.data(), &rDx);
  tree.Branch(rPeakIdxBranchName.data(), &rPeakIdx,"rPeakIdx/l");
  tree.Branch(rNumAveragedFuncsBranchName.data(), &rNumAveragedFuncs,"rNumAveragedFuncs/l");

  // ------------------------------------------------------------------
  //Save
  // ------------------------------------------------------------------
  tree.Fill();
  outputFile.Write();
  outputFile.Close();
  return 0;
  
}


int MVectorTemplate::loadTemplateFromTFile(const std::string &fullFileName)
{
  
  std::vector<std::string> branchNamesV = {rTemplateValuesBranchName.data(), rDxBranchName.data(), rPeakIdxBranchName.data(), rNumAveragedFuncsBranchName.data()};
  std::vector<double> * rTemplateValuesPointer = 0;  

  std::vector<void*> pointerV = {&rTemplateValuesPointer, &rDx, &rPeakIdx, &rNumAveragedFuncs};

  TChain * inputChain = myFuncs::openChain_setBranch(fullFileName.data(), rTreeName.data(), branchNamesV, pointerV);
  
  if ( !inputChain)
  {
    std::cout << "MVectorTemplate::loadTemplateFromTFile : Error - Could not load file " << fullFileName.data() << ". Aborting" << std::endl;
    return -1;
  }
    
  if ( inputChain->GetEntries() < 1)
  {
    std::cout << "MVectorTemplate::loadTemplateFromTFile : Error - No entries in file " << fullFileName.data() << ". Aborting" << std::endl;
    return -2;
  }

  inputChain->GetEntry(0);
  
  rTemplateValues = *rTemplateValuesPointer;

  delete rTemplateValuesPointer;
  delete inputChain;

	rTF1 = TF1("templateTF1",this,&MVectorTemplate::TF1Eval,0,rDx*(rTemplateValues.size()-1),3,"MVectorTemplate","TF1Eval");
// 	rTF1 = TF1("templateTF1",
// 				[&](double*var, double *p)
// 				{   
// 					const double effectiveX = var[0] - p[2]; 
// 					if(effectiveX <= 0) return p[1] + p[0]*rTemplateValues[0];
// 					if(effectiveX >= rDx * double(rTemplateValues.size() - 1)) return p[1] + p[0]*rTemplateValues.back();
// 					const int idx = int(effectiveX / rDx);
// 
// 					//Linear interpolation
// 					const double x1 = double(idx)*rDx;
// 					const double x2 = double(idx + 1) * rDx;
// 					const double y1 = rTemplateValues[idx];
// 					const double y2 = rTemplateValues[idx + 1];
// 					double returnVal = p[1] + p[0]*((y2 - y1)/(x2-x1)*(effectiveX-x1) + y1);
// 					
// 					return returnVal;    
// 				}//Lambda
// 				
// 				, 0, rDx*(rTemplateValues.size()-1), 3);
// 	setTF1Parameters(1.,0.,0.);
// 	setTF1ParNames();
// 		
		
	return 0;
}