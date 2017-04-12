#include "MVectorTemplate.h"

#include <iostream>

#include "TF1.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

#include "mathFuncs.h"
#include "fileFuncs.h"

#define cline std::cout << "line = " << __LINE__ << std::endl;

//For TMiuit2
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"

namespace myFuncs {
	
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
  m_minChi2_NdfLimit(-1),
  m_maxChi2_NdfLimit(3),
  m_minXshiftLimit(0),
  m_maxXshiftLimit(0), 
  m_minPedestalLimit(0),
  m_maxPedestalLimit(0),
  m_doubleNumbersEqualThershold(1e-20),
  m_debugLevel(0), 
	m_amplitudeFitEnabled(true),
  m_pedestalFitEnabled(true),
  m_xShiftFitEnabled(true),
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

void MVectorTemplate::setTF1ParNames() {
  m_tF1.SetParName(0,"Amplitude");
  m_tF1.SetParName(1,"Pedestal");
  m_tF1.SetParName(2,"xShift");
}


void MVectorTemplate::resetTemplateValues() {
  m_templateValues.clear(); 
  m_numAveragedFuncs = 0;
  m_tF1 = TF1();
}

void MVectorTemplate::resetTemplateRange() {
  if(m_tF1.IsValid()) m_tF1.SetRange(0,m_dx*(m_templateValues.size()-1) );
}

void MVectorTemplate::setDx(const double newDx) {
  m_dx = newDx;
  resetTemplateRange();
}


double MVectorTemplate::calcSimplePedestal(const std::vector<double>& vector, const double percentage) const {
	
	const long numElementsToAverage = std::distance(vector.begin(), vector.begin() + (vector.size() * percentage) );
	
	if( numElementsToAverage <= 1 ) return vector[0];
	
	return std::accumulate(vector.begin(), vector.begin() + numElementsToAverage, 0.0) / numElementsToAverage;
}


void MVectorTemplate::addFirstVector(const std::vector<double>& newVector) {
	// ------------------------------------------------------------------
	//Take pedestal as first 5% of entries
	// ------------------------------------------------------------------
	const double pedestal = calcSimplePedestal(newVector);

	// ------------------------------------------------------------------
	//Search for maximum (and what idx it is at - for later guessing the time shift) of newVector so peak can be normalized
	// ------------------------------------------------------------------
	const auto maxElementIt = std::max_element(newVector.begin(), newVector.end());
	const double maximum = *maxElementIt;
	m_peakIdx = std::distance(newVector.begin(), maxElementIt);
	
	const double maxMinusPedestal = maximum - pedestal;
	// ------------------------------------------------------------------
	//Populate m_templateValues, normalizing and subtracting pedestal
	// ------------------------------------------------------------------
	m_templateValues.clear();
	m_templateValues.reserve(newVector.size());
	
	const double ov_maxMinusPedestal = 1.0 / maxMinusPedestal;
	
	for(size_t idx = 0; idx < newVector.size(); ++idx)
		m_templateValues.push_back( (newVector[idx] - pedestal) * ov_maxMinusPedestal );
	
	if (m_debugLevel > 20)
	{
		std::cout << "MVectorTemplate::addVector : pedestal = " << pedestal << ", maximum = " << maximum << ", maxMinusPedestal = " << maxMinusPedestal << std::endl;
		for(size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
			std::cout << "MVectorTemplate::addVector : temp[" << iTemplateIdx << "] = " << m_templateValues[iTemplateIdx] << std::endl;
	}
	
	// ------------------------------------------------------------------
	//Create TF1
	// ------------------------------------------------------------------
	m_tF1 = TF1("templateTF1", this, &MVectorTemplate::TF1Eval, 0, m_dx*(m_templateValues.size()-1), 3, "MVectorTemplate", "TF1Eval");
	setTF1ParNames();
	m_tF1.SetParameters(maxMinusPedestal, pedestal, 0.0);
}

TFitResult MVectorTemplate::addVector(const std::vector<double>& newVector, const double std) {
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
  if(m_numAveragedFuncs == 0) { 
  
    addFirstVector(newVector);
		++m_numAveragedFuncs;
		return TFitResult();			
  }//First vector added
  
  else { //from second vector
  
    //Standard deviation vector
    std::vector<double> stdVector(newVector.size(), std);
    
    //Histrogram with entries of newVector - used to fit the template
    //Center of each bin corresponds to correct x value (that's why there's a -m_dx/2. term)
    TH1D fitHist("fitHist","fitHist", newVector.size(), -m_dx/2.0 , newVector.size() * m_dx - m_dx/2.0);
    for(size_t idx = 1; idx <= newVector.size(); ++idx)
    {
      fitHist.SetBinContent(idx, newVector[idx-1]);
      fitHist.SetBinError(idx, std);
    }
    
    // ------------------------------------------------------------------
    //Initial guesses for fit
    // ------------------------------------------------------------------
    
    const double pedestalGuess = calcSimplePedestal(newVector);
		
		//Search for maximum around m_peakIdx
		const int elementsBeforeAfter = newVector.size() * 0.1;
		
		long firstElement = m_peakIdx - elementsBeforeAfter;
		if(firstElement < 0) firstElement = 0;
		
		long lastElement = m_peakIdx + elementsBeforeAfter;
		if(lastElement > static_cast<long>(newVector.size()) - 1) lastElement = newVector.size() - 1;
		
		const double maxElement = *std::max_element(newVector.begin() + firstElement, newVector.begin() + lastElement);
    //Amplitude initial guess is the value of newVector where the template maximum is. This will not work for fast varying functions, or if the xShift is large
    double amplitudeGuess = (maxElement - pedestalGuess);
		    
    //Make sure amplitude guess is not out of limits
    if(m_useAmplitudeLimits) 
    {
      amplitudeGuess = std::min(amplitudeGuess, m_maxAmplitudeLimit);
      amplitudeGuess = std::max(amplitudeGuess, m_minAmplitudeLimit);
    }
    
    if(m_debugLevel > 15)
			std::cout << "amplitudeGuess = " << amplitudeGuess << ", pedestalGuess = " << pedestalGuess << std::endl;
    
		
    setTF1Parameters(amplitudeGuess, pedestalGuess, 0.0); //Amplitude, pedestal, xShift
    
//     std::cout << "amplitudeGuess = " << amplitudeGuess << ", pedestalGuess = " << pedestalGuess << ", xShiftGuess = " <<  xShiftGuess << std::endl;
    
    
    // ------------------------------------------------------------------
    //Fit!
    // ------------------------------------------------------------------
		TFitResultPtr fitResult = fitTemplate(fitHist);
		
		
		// ------------------------------------------------------------------
		//Positive fittedXshift means newVector is forward in time relative to m_templateValues, i.e. newVector[idx] = m_templateValues[idx-fittedXshift/m_dx]
		// ------------------------------------------------------------------
		//TODO implement!
		//How to check if first vector added is corrupt (i.e. flat). Maybe have a cut on maximum - pedestal.
		//Amplitude, pedestal, xShift limits/cuts
		//Chi2 cut?
		const double fittedAmplitude = fitResult->Parameter(0);
    const double fittedPedestal = fitResult->Parameter(1);
    const double fittedXshift = fitResult->Parameter(2);
		
		const int fitStatus = fitResult;
		if( fitStatus ) {
			if(fitStatus != 4000) {
				if (m_debugLevel > 0) std::cerr << "MVectorTemplate::addVector : fit failed, fit status = " << fitStatus << ". Not adding vector!" << std::endl;
				return -10;
			}
			else if(m_debugLevel > 0) std::cout << "MVectorTemplate::addVector : MINOS falied, fit status = " << fitStatus << ". Vector MIGHT still be added!!!" << std::endl;
		}
		
		const double chi2_ndf = fitResult->Chi2() / fitResult->Ndf();
		
		if( chi2_ndf < m_minChi2_NdfLimit || chi2_ndf > m_maxChi2_NdfLimit ) 
    {
      if(m_debugLevel > 0) std::cerr << "MVectorTemplate::addVector : fit failed, Chi2/NDF of fit = " << chi2_ndf << " is our of limits. limits are (" << m_minChi2_NdfLimit << "," << m_maxChi2_NdfLimit << "). Not adding vector!" << std::endl;
      return -30; 
    }
    
    if(m_debugLevel > 10) std::cout << "MVectorTemplate::addVector : fittedAmplitude = " << fittedAmplitude << ", fittedPedestal = " << fittedPedestal << ", fittedXshift = " << fittedXshift << ", chi2 / NDF = " << chi2_ndf << std::endl;
    
    
    if( m_useAmplitudeLimits and (std::abs(fittedAmplitude - m_minAmplitudeLimit) < m_doubleNumbersEqualThershold || std::abs(fittedAmplitude - m_maxAmplitudeLimit) < m_doubleNumbersEqualThershold) )
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached amplitude limit. fittedAmplitude = " << fittedAmplitude << ". limits are (" << m_minAmplitudeLimit << "," << m_maxAmplitudeLimit << "). Not adding vector!" << std::endl;
      return -11;
    }
    
    if( m_useXshiftLimits and (std::abs(fittedXshift - m_minXshiftLimit) < m_doubleNumbersEqualThershold || std::abs(fittedXshift - m_maxXshiftLimit) < m_doubleNumbersEqualThershold) )
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached xShift limit. fittedXshift = " << fittedXshift << ". limits are (" << m_minXshiftLimit << "," << m_maxXshiftLimit << "). Not adding vector!" << std::endl;
      return -12;
    }
    
    if(m_usePedestalLimits and (std::abs(fittedPedestal - m_minPedestalLimit) < m_doubleNumbersEqualThershold || std::abs(fittedPedestal - m_maxPedestalLimit) < m_doubleNumbersEqualThershold) )
    {
      if(m_debugLevel > 10) std::cerr << "MVectorTemplate::addVector : fit failed, reached pedestal limit. fittedPedestal = " << fittedPedestal << ". limits are (" << m_minPedestalLimit << "," << m_maxPedestalLimit << "). Not adding vector!" << std::endl;
      return -13;
    }

    // ------------------------------------------------------------------
    //Average template with new waveform
    //x corresponding to newVector[0] is: - fittedXshift (with minus!!)
    //x corresponding to newVector[size] is size*m_dx - fittedXshift
    // ------------------------------------------------------------------
    
    //Find which are the first and last indexes in the template which needs to be averaged.
    //firstTemplateIdx,lastTemplatedIdx will be averaged with the interpolated value from newVector
    size_t firstTemplateIdx = 0;
		size_t lastTemplateIdx = 0;
    if(fittedXshift > 0)
    {
      //template:                      Idx0---------Idx1---------Idx2          Idx(size-4)-------------Idx(size-3)-------------Idx(size-2)--
      //newVector: Idx0---------Idx1---------Idx2---------Idx3     Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
      //In this example both vectors are same length and fittedXshift = 1.5m_dx
      //firstTemplateIdx = 0, lastTemplateIdx = newVector.size()-3
      firstTemplateIdx = 0; 
      lastTemplateIdx = newVector.size() - 1 - static_cast<size_t>(fittedXshift / m_dx) - 1; 
    }
    else
    {
      //template:  Idx0--------Idx1--------Idx2          Idx(size-2)-------------Idx(size-1)
      //newVector:                   Idx0--------Idx1                Idx(size-3)-------------Idx(size-2)-------------Idx(size-1)
      //In this example both vectors are same length and fittedXshift = 1.5m_dx
      //firstTemplateIdx = 2, lastTemplateIdx = m_templateValues.size()-1
      firstTemplateIdx = static_cast<size_t>(-fittedXshift / m_dx) + 1;
      lastTemplateIdx = m_templateValues.size() - 1; //last element
    }

    // ------------------------------------------------------------------
    //Average newVector with template
    // ------------------------------------------------------------------
    
    //Create template of newVector
    MVectorTemplate newVectorTemplate(newVector, m_dx);
    TF1* newVector_tF1 = newVectorTemplate.getTF1();
		
		
		//We now want to average the template with the TF1 created from the new vector.
		//We first set the pedestal to zero, because the template has pedestal of zero.
		newVector_tF1->SetParameter(1, 0.0); //set pedestal to 0.
		
		//We renormalize the new TF1. The newVector_tF1->GetParameter(0) term "brings back" the tf1 to it's original height.
		//The / fittedAmplitude term scales it by how much it's similar to the template.
		//I.e. if tF1->GetParameter(0) > fittedAmplitude, this means the new vector is higher than the current template. This will put the peak of the new vector above the current template, as it should.
		if(newVector_tF1->GetParameter(0) != 0.0 ) 
			newVector_tF1->SetParameter(0, newVector_tF1->GetParameter(0) / fittedAmplitude );
		
		TCanvas c("c","c",0,0,1200,900);
		fitHist->Draw();
		newVector_tF1->SetParameter(2, -fittedXshift); 

		
		newVector_tF1->SetNpx(1000);
		
		newVector_tF1->Draw("SAME");
		
		gPad->WaitPrimitive();
		
    const double ov_numAveragedFuncsP1 = 1.0 / static_cast<double>(m_numAveragedFuncs + 1);
		
    //Average
    double peakVal = -DBL_MAX;
		
    for(size_t iTemplateIdx = firstTemplateIdx; iTemplateIdx <= lastTemplateIdx; ++iTemplateIdx)
    {
      const double newVectorValue = newVector_tF1->Eval(iTemplateIdx * m_dx);
      
      //average vectors
      m_templateValues[iTemplateIdx] = (m_templateValues[iTemplateIdx] * static_cast<double>(m_numAveragedFuncs) + newVectorValue) * ov_numAveragedFuncsP1;
      
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
    if(fittedXshift > 0){
			if(m_templateValues.size() > lastTemplateIdx + 1) {
				m_templateValues.resize( lastTemplateIdx + 1); //lastTemplateIdx + 1 == size of new template vector;
				resetTemplateRange();
			}			
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
    const double pedestal = calcSimplePedestal(m_templateValues);
    for(size_t iTemplateIdx = 0; iTemplateIdx < m_templateValues.size(); ++iTemplateIdx)
      m_templateValues[iTemplateIdx] = (m_templateValues[iTemplateIdx] - pedestal )/ (peakVal-pedestal);
        
//     m_tF1.SetParameters(1,0,0);
		
		++m_numAveragedFuncs;
		
		return *fitResult;
  }//add, from second vector
  
  
  
  return TFitResult();
}

TFitResultPtr MVectorTemplate::fitTemplate(TH1D& fitHist) {
	
	//Attempt to fit at current x.
	
	TFitResultPtr fitResult = fitHist.Fit(&m_tF1,"QSN");
	//"Q" Quiet mode (minimum printing)
	//"M" More. Improve fit results. It uses the IMPROVE command of TMinuit (see TMinuit::mnimpr). This algorithm attempts to improve the found local minimum by searching for a better one.
	//"N" Do not store the graphics function, do not draw
	//"S" The result of the fit is returned in the TFitResultPtr (see below Access to the Fit Result)
	
	if ( fitGood(fitResult) )
		return fitResult;
	
	//Fit failed. Move xShift and retry
	double xShiftGuess = m_dx * m_templateValues.size() * 0.05;
	if(m_debugLevel >=20) 
		std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
	m_tF1.SetParameter(2, xShiftGuess);
	fitResult = fitHist.Fit(&m_tF1,"QSN");
	if ( fitGood(fitResult) )
		return fitResult;
	
	//move xShift and retry
	xShiftGuess = -m_dx * m_templateValues.size() * 0.05 ;
	if(m_debugLevel >=20) 
		std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
	m_tF1.SetParameter(2, xShiftGuess);
	fitResult = fitHist.Fit(&m_tF1,"QSN");
	if ( fitGood(fitResult) )
		return fitResult;
	
	xShiftGuess = m_dx * m_templateValues.size() * 0.1;
	if(m_debugLevel >=20) 
		std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
	m_tF1.SetParameter(2, xShiftGuess);
	fitResult = fitHist.Fit(&m_tF1,"QSN");
	if ( fitGood(fitResult) )
		return fitResult;
	
	//move xShift and retry
	xShiftGuess = -m_dx * m_templateValues.size() * 0.1;
	if(m_debugLevel >=20) 
		std::cout << "setting xShift Guess to " << xShiftGuess << std::endl;
	m_tF1.SetParameter(2, xShiftGuess);
	fitResult = fitHist.Fit(&m_tF1,"QSN");
	
	//Return fitResult in any case
	return fitResult;
}

bool MVectorTemplate::fitGood(const TFitResultPtr& fitResult) const {
	
	const int fitStatus = fitResult;
	
	//Fit status
	if( fitStatus != 0 and fitStatus != 4000) //4000 is MINOS problem
			return false;
	
	// chi2/ndf 
	const double chi2_ndf = fitResult->Chi2() / fitResult->Ndf();
	if(chi2_ndf < m_minChi2_NdfLimit or chi2_ndf > m_maxChi2_NdfLimit)
		return false;

	return true;
}

void MVectorTemplate::enableAmplitudeFit(const bool enableAmplitudeFit) {
  m_amplitudeFitEnabled = enableAmplitudeFit;
  setTF1Parameters();
}

void MVectorTemplate::enablePedestalFit(const bool enablePedestalFit) {
  m_pedestalFitEnabled = enablePedestalFit;
  setTF1Parameters();
}

void MVectorTemplate::enableXshiftFit(const bool enableXshiftFit) {
  m_xShiftFitEnabled = enableXshiftFit;
  setTF1Parameters();
}

double MVectorTemplate::TF1Eval(double *var, double *params) {
//   double amplitude = params[0];
//   double pedestal = params[1];
//   double xShift = params[2];
	
  const double effectiveX = var[0] - params[2]; 
  
  if(effectiveX <= 0) return params[1];// + params[0]*m_templateValues[0];
  if(effectiveX >= m_dx * static_cast<double>(m_templateValues.size() - 1)) return params[1] + params[0]*m_templateValues.back();

  const int idx = static_cast<int>(effectiveX / m_dx);

  //Linear interpolation
  const double x1 = static_cast<double>(idx) * m_dx;
	const double x2 = static_cast<double>(idx + 1) * m_dx;
//     cout << "vec[idx + 1] = " << vec[idx + 1];
  const double y1 = m_templateValues[idx];
	const double y2 = m_templateValues[idx + 1];

	return params[1] + params[0]*((y2 - y1)/(x2-x1)*(effectiveX-x1) + y1);
}

void MVectorTemplate::setTF1Parameters(const double amplitude, const double pedestal, const double xShift) {
  if( !m_tF1.IsValid() ) return;
  
  if(m_amplitudeFitEnabled) {
    m_tF1.SetParameter(0, amplitude);
    m_tF1.SetParLimits(0, m_minAmplitudeLimit,m_maxAmplitudeLimit);
  }
  else m_tF1.FixParameter(0, amplitude);
  
  if(m_pedestalFitEnabled) {
    m_tF1.SetParameter(1, pedestal);
    m_tF1.SetParLimits(1, m_minPedestalLimit,m_maxPedestalLimit);
  }
  else m_tF1.FixParameter(1, pedestal);
  
	if(m_xShiftFitEnabled) {
		m_tF1.SetParameter(2, xShift);
		m_tF1.SetParLimits(2, m_minXshiftLimit, m_maxXshiftLimit);
	}
	else m_tF1.FixParameter(2, xShift);
}

void MVectorTemplate::setTF1Parameters() {
  if( !m_tF1.IsValid() ) return;
  setTF1Parameters(m_tF1.GetParameter(0), m_tF1.GetParameter(1), m_tF1.GetParameter(2));
}

void MVectorTemplate::setAmplitudeLimits(const double newMinAmplitudeLimit, const double newMaxAmplitudeLimit) {
  if( std::abs(newMinAmplitudeLimit - newMaxAmplitudeLimit) > m_doubleNumbersEqualThershold) m_useAmplitudeLimits = true;
  else m_useAmplitudeLimits = false;
  m_minAmplitudeLimit = newMinAmplitudeLimit; 
  m_maxAmplitudeLimit = newMaxAmplitudeLimit;
  setTF1Parameters();
}

void MVectorTemplate::setXshiftLimits(const double newMinXshiftLimit, const double newMaxXshiftLimit) {
  if( std::abs(newMinXshiftLimit - newMaxXshiftLimit) > m_doubleNumbersEqualThershold) m_useXshiftLimits = true;
  else m_useXshiftLimits = false;
  m_minXshiftLimit = newMinXshiftLimit; 
  m_maxXshiftLimit = newMaxXshiftLimit;
  setTF1Parameters();
}

void MVectorTemplate::setPedestalLimits(const double newMinPedestalLimit, const double newMaxPedestalLimit) {
  if( std::abs(newMinPedestalLimit - newMaxPedestalLimit) > m_doubleNumbersEqualThershold) m_usePedestalLimits = true;
  else m_usePedestalLimits = false;
  m_minPedestalLimit = newMinPedestalLimit; 
  m_maxPedestalLimit = newMaxPedestalLimit;
  setTF1Parameters();  
}

int MVectorTemplate::saveTemplateToTFile(const std::string& fullFileName, const std::string& treeDescription) {
  
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


int MVectorTemplate::loadTemplateFromTFile(const std::string& fullFileName) {
  
  std::vector<std::string> branchNamesV = {m_templateValuesBranchName.data(), m_dxBranchName.data(), m_peakIdxBranchName.data(), m_numAveragedFuncsBranchName.data()};
  std::vector<double>* m_templateValuesPointer = 0;  

  std::vector<void*> pointerV = {&m_templateValuesPointer, &m_dx, &m_peakIdx, &m_numAveragedFuncs};

  std::unique_ptr<TChain> inputChain (myFuncs::openChain_setBranch(fullFileName.data(), m_treeName.data(), branchNamesV, pointerV));
  
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

} //namespace myFuncs