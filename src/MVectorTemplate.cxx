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
	m_tF1(TF1()),
	m_treeName("vectorTemplate"),
  m_templateValuesBranchName("m_templateValues"),
  m_dxBranchName("rDx"),
  m_numAveragedFuncsBranchName("rNumAveragedFuncs"),
  m_peakIdxBranchName("m_peakIdx"),
  m_numAveragedFuncs(0), 
  m_peakIdx(0),
	m_dx(0.), 
  m_minAmplitudeLimit(0),    
  m_maxAmplitudeLimit(-1),
  m_minChi2Limit(-1),
  m_maxChi2Limit(DBL_MAX),
  m_minXshiftLimit(0),
  m_maxXshiftLimit(0), 
  m_minPedestalLimit(0),
  m_maxPedestalLimit(0),
  m_doubleNumbersEqualThershold(1e-20),
  m_debugLevel(0), 
	m_amplitudeFitEnabled(true),
  m_pedestalFitEnabled(true),
  m_useAmplitudeLimits(false),
  m_useXshiftLimits(false),
  m_usePedestalLimits(false)
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
  m_tF1.SetParName(0,"Amplitude");
  m_tF1.SetParName(1,"Pedestal");
  m_tF1.SetParName(2,"TimeShift");
}


void MVectorTemplate::resetTemplateValues()
{
  m_templateValues.clear(); 
  m_numAveragedFuncs = 0;
  m_tF1 = TF1();
}

void MVectorTemplate::resetTemplateRange()
{
  if(m_tF1.IsValid()) m_tF1.SetRange(0,m_dx*(m_templateValues.size()-1) );
}

void MVectorTemplate::setDx(const double newDx)
{
  m_dx = newDx;
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
  if(m_numAveragedFuncs == 0) //First vector added
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
				m_peakIdx = idx;
      }
    }
    double maxMinusPedestal = maximum - pedestal;
    // ------------------------------------------------------------------
    //Populate m_templateValues, normalizing and subtracting pedestal
    // ------------------------------------------------------------------
    m_templateValues = std::vector<double>(newVector.size());
    for(size_t idx = 0; idx < newVector.size(); ++idx)
    {
      m_templateValues[idx] = (newVector[idx] - pedestal)/maxMinusPedestal;
    }
    
    if (m_debugLevel > 20)
    {
      std::cout << "MVectorTemplate::addVector : pedestal = " << pedestal << ", maximum = " << maximum << std::endl;
      for(size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
				std::cout << "MVectorTemplate::addVector : temp[" << iTemplateIdx << "] = " << m_templateValues[iTemplateIdx] << std::endl;
    }
    
    // ------------------------------------------------------------------
    //Create TF1
    // ------------------------------------------------------------------
    m_tF1 = TF1("templateTF1",this,&MVectorTemplate::TF1Eval,0,m_dx*(m_templateValues.size()-1),3,"MVectorTemplate","TF1Eval");
    setTF1ParNames();
  }//First vector added
  
  else //from second vector
  {
    
    //Standard deviation vector
    std::vector<double> stdVector(newVector.size(),std);
    
    //Histrogram of with entries of newVector - used to fit the template
    //Center of each bin corresponds to correct x value (that's why there's a -m_dx/2. term)
    TH1D fitHist("fitHist","fitHist", newVector.size(), -m_dx/2. , newVector.size() * m_dx - m_dx/2.);
    for(unsigned int idx = 1; idx <= newVector.size(); ++idx)
    {
      fitHist.SetBinContent(idx,newVector[idx-1]);
      fitHist.SetBinError(idx,std);
    }
    
    // ------------------------------------------------------------------
    //Initial guesses for fit
    // ------------------------------------------------------------------
    //Amplitude initial guess is the value of newVector where the template maximum is. This will not work for fast varying functions, or if the xShift is large
    double amplitudeGuess = newVector[m_peakIdx] / m_templateValues[m_peakIdx]; 
    
    //Make sure amplitude guess is not out of limits
    if(m_useAmplitudeLimits) 
    {
      amplitudeGuess = std::min(amplitudeGuess, m_maxAmplitudeLimit);
      amplitudeGuess = std::max(amplitudeGuess, m_minAmplitudeLimit);
    }
    double pedestalGuess = newVector[0];
    double xShiftGuess = 0;
    
    setTF1Parameters(amplitudeGuess, pedestalGuess, xShiftGuess);
    
//     std::cout << "amplitudeGuess = " << amplitudeGuess << ", pedestalGuess = " << pedestalGuess << ", xShiftGuess = " <<  xShiftGuess << std::endl;
    
    
    // ------------------------------------------------------------------
    //Fit!
    // ------------------------------------------------------------------
    TFitResultPtr fitResult = fitHist.Fit(&m_tF1,"QSN");
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
				if (m_debugLevel > 0) std::cerr << "MVectorTemplate::addVector : fit failed, fit status = " << fitStatus << ". Not adding vector!" << std::endl;
				return -10;
      }
      else if(m_debugLevel > 0) std::cout << "MVectorTemplate::addVector : MINOS falied, fit status = " << fitStatus << ". Vector MIGHT still be added!!!" << std::endl;
    }
      
    // ------------------------------------------------------------------
    //Positive fittedXshift means newVector is forward in time relative to m_templateValues, i.e. newVector[idx] = m_templateValues[idx-fittedXshift/m_dx]
    // ------------------------------------------------------------------
    double fittedAmplitude = fitResult->Parameter(0);
    double fittedPedestal = fitResult->Parameter(1);
    double fittedXshift = fitResult->Parameter(2);
    
    double chi2_NDF = fitResult->Chi2();
       
    if(m_debugLevel > 10) std::cout << "MVectorTemplate::addVector : fittedAmplitude = " << fittedAmplitude << ", fittedPedestal = " << fittedPedestal << ", fittedXshift = " << fittedXshift << ", chi2 / NDF = " << chi2_NDF << std::endl;
    
    
    if( (std::abs(fittedAmplitude - m_minAmplitudeLimit) < m_doubleNumbersEqualThershold || std::abs(fittedAmplitude - m_maxAmplitudeLimit) < m_doubleNumbersEqualThershold)  && m_useAmplitudeLimits)
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached amplitude limit. fittedAmplitude = " << fittedAmplitude << ". limits are (" << m_minAmplitudeLimit << "," << m_maxAmplitudeLimit << "). Not adding vector!" << std::endl;
      return -11;
    }
    
    if( (std::abs(fittedXshift - m_minXshiftLimit) < m_doubleNumbersEqualThershold || std::abs(fittedXshift - m_maxXshiftLimit) < m_doubleNumbersEqualThershold) && m_useXshiftLimits)
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached xShift limit. fittedXshift = " << fittedXshift << ". limits are (" << m_minXshiftLimit << "," << m_maxXshiftLimit << "). Not adding vector!" << std::endl;
      return -12;
    }
    
    if( (std::abs(fittedPedestal - m_minPedestalLimit) < m_doubleNumbersEqualThershold || std::abs(fittedPedestal - m_maxPedestalLimit) < m_doubleNumbersEqualThershold) && m_usePedestalLimits)
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached pedestal limit. fittedPedestal = " << fittedPedestal << ". limits are (" << m_minPedestalLimit << "," << m_maxPedestalLimit << "). Not adding vector!" << std::endl;
      return -13;
    }
    
    if( chi2_NDF < m_minChi2Limit || chi2_NDF > m_maxChi2Limit ) 
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, Chi2/NDF of fit = " << chi2_NDF << " is too large. limits are (" << m_minChi2Limit << "," << m_maxChi2Limit << "). Not adding vector!" << std::endl;
      return -30; 
    }
    
//     if(chi2_NDF > 2) std::cout << "chi2 / NDF = " << chi2_NDF << std::endl;

    // ------------------------------------------------------------------
    //Average template with new waveform
    //x corresponding to newVector[0] is: - fittedXshift (with minus!!)
    //x corresponding to newVector[size] is size*m_dx - fittedXshift
    // ------------------------------------------------------------------
    
    //Find which are the first and last indexes in the template which needs to be averaged.
    //firstTemplateIdx,lastTemplatedIdx will be averaged with the interpolated value from newVector
    size_t firstTemplateIdx = 0, lastTemplateIdx = 0;
    if(fittedXshift > 0)
    {
      //template:                      Idx0---------Idx1---------Idx2          Idx(size-4)-------------Idx(size-3)-------------Idx(size-2)--
      //newVector: Idx0---------Idx1---------Idx2---------Idx3     Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
      //In this example both vectors are same length and fittedXshift = 1.5m_dx
      //firstTemplateIdx = 0, lastTemplateIdx = newVector.size()-3
      firstTemplateIdx = 0; 
      lastTemplateIdx = newVector.size() - 1 - size_t(fittedXshift / m_dx) - 1; 
    }
    else
    {
      //template:  Idx0--------Idx1--------Idx2          Idx(size-2)-------------Idx(size-1)
      //newVector:                   Idx0--------Idx1                Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
      //In this example both vectors are same length and fittedXshift = 1.5m_dx
      //firstTemplateIdx = 2, lastTemplateIdx = m_templateValues.size()-1
      firstTemplateIdx = size_t(-fittedXshift / m_dx) + 1;
      lastTemplateIdx = m_templateValues.size() - 1; //last element
    }

    // ------------------------------------------------------------------
    //Average newVector with template
    // ------------------------------------------------------------------
    
    //Create template of newVector
    MVectorTemplate newVectorTemplate(newVector, m_dx);
    TF1 * newVectom_tF1 = newVectorTemplate.getTF1();
    
    double originalNewVectorMax = newVector[ newVectorTemplate.getPeakIdx()];
    double newVectorPedestal = (newVector[0] - fittedPedestal)/(originalNewVectorMax-newVector[0]);;
    
    double newVectorAmplitude = (originalNewVectorMax - newVector[0]) / fittedAmplitude;
    
    newVectom_tF1->SetParameters(newVectorAmplitude, newVectorPedestal, -fittedXshift);
    
    //Average
    double peakVal = -DBL_MAX;
    for(size_t iTemplateIdx = firstTemplateIdx; iTemplateIdx <= lastTemplateIdx; ++iTemplateIdx)
    {
      double newVectorValue = newVectom_tF1->Eval(iTemplateIdx * m_dx);
      
      //average vectors
      m_templateValues[iTemplateIdx] = (m_templateValues[iTemplateIdx]*double(m_numAveragedFuncs) + newVectorValue)/double(m_numAveragedFuncs+1);
      
      //Store maximum value (used to normalize later) and its index
      if(m_templateValues[iTemplateIdx] > peakVal)
      {
	peakVal = m_templateValues[iTemplateIdx];
	m_peakIdx = iTemplateIdx;
      }
      
      if (m_debugLevel > 20)
	std::cout << "MVectorTemplate::addVector : temp[" << iTemplateIdx << "] = " << m_templateValues[iTemplateIdx] << std::endl;
    }
    
    // ------------------------------------------------------------------
    //Clip ends of template
    //x corresponding to newVector[0] is - fittedXshift (with minus!!)
    //x corresponding to newVector[size] is size*m_dx - fittedXshift
    // ------------------------------------------------------------------
    if(fittedXshift > 0)
    {
      m_templateValues.resize( lastTemplateIdx + 1); //lastTemplateIdx + 1 == size of new template vector;
      resetTemplateRange();
    }
    else
    {
      if(firstTemplateIdx> 0 ) m_templateValues.erase( m_templateValues.begin() , m_templateValues.begin() + firstTemplateIdx - 1 );
      resetTemplateRange();
    }
    
    // ------------------------------------------------------------------
    //re-Normalize vector - keep peak == 1
    //Make sure first entry is zero
    // ------------------------------------------------------------------
    double pedestal = m_templateValues[0];
    for(size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
    {
      m_templateValues[iTemplateIdx] = (m_templateValues[iTemplateIdx] - pedestal )/ (peakVal-pedestal);
    }
    
  }//add, from second vector
  
  ++m_numAveragedFuncs;
  m_tF1.SetParameters(1,0,0);
  return 0;
}

void MVectorTemplate::enableAmplitudeFit(const bool enableAmplitudeFit)
{
  m_amplitudeFitEnabled = enableAmplitudeFit;
  setTF1Parameters();
}

// double MVectorTemplate::eval(const double &x, const double &amplitude, const double &pedestal, const double &xShift) const
// {
//   if( !m_tF1.IsValid())
//   {
//     std::cerr << "MVectorTemplate::eval : TF1 not valid" << std::endl;
//     return 0.;
//   }
//   
//   double originalAmplitude = m_tF1.GetParameter(0);
//   double originalPedestal = m_tF1.GetParameter(1);
//   double originalXshift = m_tF1.GetParameter(2);
//   
//   m_tF1.SetParameters(amplitude,pedestal,xShift);
//   double returnVal = m_tF1.Eval(x);
//   m_tF1.SetParameters(originalAmplitude,originalPedestal,originalXshift);
//   
//   return returnVal;
//   
// }

void MVectorTemplate::enablePedestalFit(const bool enablePedestalFit) 
{
  m_pedestalFitEnabled = enablePedestalFit;
  setTF1Parameters();
}

double MVectorTemplate::TF1Eval(double *var, double *params) 
{
//   double amplitude = params[0];
//   double pedistal = params[1];
//   double xShift = params[2];
  double effectiveX = var[0] - params[2]; 
  
//         cout << "amplitude = " << amplitude << ", pedistal = " << pedistal << ", xShift = " << xShift << ", x = " << var[0] << ", effectiveX = " << effectiveX << ". ";
  
  if(effectiveX <= 0) return params[1] + params[0]*m_templateValues[0];
  if(effectiveX >= m_dx * double(m_templateValues.size() - 1)) return params[1] + params[0]*m_templateValues.back();
  
  int idx = int(effectiveX / m_dx);

  //Linear interpolation
  double x1 = double(idx)*m_dx, x2 = double(idx + 1) * m_dx;
//     cout << "vec[idx + 1] = " << vec[idx + 1];
  double y1 = m_templateValues[idx], y2 = m_templateValues[idx + 1];
  
//     std::cout << "idx = " << idx << ", x1 = " << x1 << ", y1 = " << y1 << ". x2 = " << x2 << ", y2 = " << y2 << std::endl;
//   double effY = (y2 - y1)/(x2-x1)*(effectiveX-x1) + y1;
  
  double returnVal = params[1] + params[0]*((y2 - y1)/(x2-x1)*(effectiveX-x1) + y1);
  
//     std::cout << "returnVal = " << returnVal << std::endl;
  return returnVal;    
}

void MVectorTemplate::setTF1Parameters(const double amplitude, const double pedestal, const double xShift)
{
  if( !m_tF1.IsValid() ) return;
  
  if(m_amplitudeFitEnabled) 
  {
    m_tF1.SetParameter(0, amplitude);
    m_tF1.SetParLimits(0, m_minAmplitudeLimit,m_maxAmplitudeLimit);
  }
  else m_tF1.FixParameter(0, amplitude);
  
  if(m_pedestalFitEnabled)
  {
    m_tF1.SetParameter(1, pedestal);
    m_tF1.SetParLimits(1, m_minPedestalLimit,m_maxPedestalLimit);
  }
  else m_tF1.FixParameter(1, pedestal);
  
  m_tF1.SetParameter(2, xShift);
  m_tF1.SetParLimits(2,m_minXshiftLimit, m_maxXshiftLimit);
}

void MVectorTemplate::setTF1Parameters()
{
  if( !m_tF1.IsValid() ) return;
  setTF1Parameters(m_tF1.GetParameter(0), m_tF1.GetParameter(1), m_tF1.GetParameter(2));
}

void MVectorTemplate::setAmplitudeLimits(const double newMinAmplitudeLimit, const double newMaxAmplitudeLimit)
{
  if( std::abs(newMinAmplitudeLimit - newMaxAmplitudeLimit) > m_doubleNumbersEqualThershold) m_useAmplitudeLimits = true;
  else m_useAmplitudeLimits = false;
  m_minAmplitudeLimit = newMinAmplitudeLimit; 
  m_maxAmplitudeLimit = newMaxAmplitudeLimit;
  setTF1Parameters();
}

void MVectorTemplate::setXshiftLimits(const double newMinXshiftLimit, const double newMaxXshiftLimit)
{
  if( std::abs(newMinXshiftLimit - newMaxXshiftLimit) > m_doubleNumbersEqualThershold) m_useXshiftLimits = true;
  else m_useXshiftLimits = false;
  m_minXshiftLimit = newMinXshiftLimit; 
  m_maxXshiftLimit = newMaxXshiftLimit;
  setTF1Parameters();
}

void MVectorTemplate::setPedestalLimits(const double newMinPedestalLimit, const double newMaxPedestalLimit)
{
  if( std::abs(newMinPedestalLimit - newMaxPedestalLimit) > m_doubleNumbersEqualThershold) m_usePedestalLimits = true;
  else m_usePedestalLimits = false;
  m_minPedestalLimit = newMinPedestalLimit; 
  m_maxPedestalLimit = newMaxPedestalLimit;
  setTF1Parameters();  
}

int MVectorTemplate::saveTemplateToTFile(const std::string& fullFileName, const std::string& treeDescription)
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
  TTree tree(m_treeName.data(), treeDescription.data());
  
  // ------------------------------------------------------------------
  //Create pointers to variables that will be saved
  // ------------------------------------------------------------------
  std::vector<double>* templateValuesPointer = &m_templateValues;
  
  // ------------------------------------------------------------------
  //Create branches
  // ------------------------------------------------------------------
  tree.Branch(m_templateValuesBranchName.data(), "std::vector<double>", &templateValuesPointer);
  tree.Branch(m_dxBranchName.data(), &m_dx);
  tree.Branch(m_peakIdxBranchName.data(), &m_peakIdx,"m_peakIdx/l");
  tree.Branch(m_numAveragedFuncsBranchName.data(), &m_numAveragedFuncs,"m_numAveragedFuncs/l");

  // ------------------------------------------------------------------
  //Save
  // ------------------------------------------------------------------
  tree.Fill();
  outputFile.Write();
  outputFile.Close();
  return 0;
  
}


int MVectorTemplate::loadTemplateFromTFile(const std::string& fullFileName)
{
  
  std::vector<std::string> branchNamesV = {m_templateValuesBranchName.data(), m_dxBranchName.data(), m_peakIdxBranchName.data(), m_numAveragedFuncsBranchName.data()};
  std::vector<double>* m_templateValuesPointer = 0;  

  std::vector<void*> pointerV = {&m_templateValuesPointer, &m_dx, &m_peakIdx, &m_numAveragedFuncs};

  TChain* inputChain = myFuncs::openChain_setBranch(fullFileName.data(), m_treeName.data(), branchNamesV, pointerV);
  
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
  
	//Prevent clang static analyzer from thinking this derefrences null pointer - makes this line invisible to the analyzer
	#ifndef __clang_analyzer__
  m_templateValues = *m_templateValuesPointer;
	#endif
	
  delete m_templateValuesPointer;
  delete inputChain;

	m_tF1 = TF1("templateTF1", this, &MVectorTemplate::TF1Eval, 0, m_dx*(m_templateValues.size()-1), 3, "MVectorTemplate", "TF1Eval");
// 	m_tF1 = TF1("templateTF1",
// 				[&](double*var, double *p)
// 				{   
// 					const double effectiveX = var[0] - p[2]; 
// 					if(effectiveX <= 0) return p[1] + p[0]*m_templateValues[0];
// 					if(effectiveX >= m_dx * double(m_templateValues.size() - 1)) return p[1] + p[0]*m_templateValues.back();
// 					const int idx = int(effectiveX / m_dx);
// 
// 					//Linear interpolation
// 					const double x1 = double(idx)*m_dx;
// 					const double x2 = double(idx + 1) * m_dx;
// 					const double y1 = m_templateValues[idx];
// 					const double y2 = m_templateValues[idx + 1];
// 					double returnVal = p[1] + p[0]*((y2 - y1)/(x2-x1)*(effectiveX-x1) + y1);
// 					
// 					return returnVal;    
// 				}//Lambda
// 				
// 				, 0, m_dx*(m_templateValues.size()-1), 3);
// 	setTF1Parameters(1.,0.,0.);
// 	setTF1ParNames();
// 		
		
	return 0;
}