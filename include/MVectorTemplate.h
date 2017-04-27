#pragma once

#include <vector>
#include "TF1.h"

class TFitResultPtr;
class TFitResult;
class TH1D;

namespace myFuncs{
	
class MVectorTemplate
{

  public:
    MVectorTemplate();
    MVectorTemplate(const std::vector<double>& newVec, const double newDx);
    
    
    inline unsigned int getNumAveragedFuncs() const {return m_numAveragedFuncs;}
    //Should be used with care
    inline void setNumAveragedFuncs(const unsigned int numAveraged) {m_numAveragedFuncs = numAveraged;}
    
    inline double getDx() const {return m_dx;}
    void setDx(const double newDx);
    
    //inline void setVectorValues(const std::vector<double> &newValues) {m_templateValues = newValues; m_numAveragedFuncs = 0;}
    inline std::vector<double> getTemplateValues() const {return m_templateValues;}
    
    void resetTemplateValues(); 
    
    // ------------------------------------------------------------------
    //enableAmplitudeFit
    //Currently no way of actually setting amplitude of fit. It is guessed when a vector is added.
    // ------------------------------------------------------------------
    void enableAmplitudeFit(const bool enableAmplitudeFit);
    inline bool isAmplitudeFitEnabled() const {return m_amplitudeFitEnabled;}
    
    // ------------------------------------------------------------------
    //enablePedestalFit
    //Currently no way of actually setting pedestal of fit. It is guessed when a vector is added.
    // ------------------------------------------------------------------
    void enablePedestalFit(const bool enablePedestalFit);
    inline bool isPedestalFitEnabled() const {return m_pedestalFitEnabled;}
    
    void enableXshiftFit(const bool enableXshiftFit);
    inline bool isXshiftFitEnabled() const {return m_xShiftFitEnabled;}
    
    inline unsigned int getTemplateSize() const {return m_templateValues.size();}
    
    
    //Adds a new vector to the template. The values are averaged (taking into account the relative wight of the new vector according to the number of previous vectors added.
    //m_tF1 parameters reset at the end (otherwise they keep values of last vector added.
    TFitResult addVector(const std::vector<double>& newVector, const double std);
    
    inline TF1* getTF1() {return &m_tF1;}
    
    inline int getDebugLevel() const {return m_debugLevel;}
    inline void setDebugLevel(const int newDebugLevel) {m_debugLevel = newDebugLevel;}
    
    inline void getAmplitudeLimits(double& minAmplitude, double& maxAmplitude) const {minAmplitude = m_minAmplitudeLimit; maxAmplitude = m_maxAmplitudeLimit;}
    void setAmplitudeLimits(const double newMinAmplitudeLimit, const double newMaxAmplitudeLimit);
    
    inline void getChi2_NdfLimits(double& minChi2_Ndf, double& maxChi2_Ndf) const {minChi2_Ndf = m_minChi2_NdfLimit; maxChi2_Ndf = m_maxChi2_NdfLimit;}
    void setChi2_NdfLimits(const double newMinChi2_NdfLimit, const double newMaxChi2_NdfLimit) { m_minChi2_NdfLimit = newMinChi2_NdfLimit; m_maxChi2_NdfLimit = newMaxChi2_NdfLimit;}
    
    inline void getXshiftLimits(double& minXshift, double& maxXshift) const {minXshift = m_minXshiftLimit; maxXshift = m_maxXshiftLimit;}
    void setXshiftLimits(const double newMinXshiftLimit, const double newMaxXshiftLimit);
    
    inline void getPedestalLimits(double& minPedestal, double& maxPedestal) const {minPedestal = m_minPedestalLimit; maxPedestal = m_maxPedestalLimit;}
    void setPedestalLimits(const double newMinPedestalLimit, const double newMaxPedestalLimit);
    
    inline size_t getPeakIdx() const {return m_peakIdx;}
    
    int saveTemplateToTFile(const std::string& fullFileName, const std::string& treeDescription);
       
    int loadTemplateFromTFile(const std::string& fullFileName);
    
		void setTF1Parameters(const double amplitude, const double pedestal, const double xShift);
		
		//This is used to keep track where the x axis 0 is.
		//It can change if items are removed from the beggining of the template, for example.
		inline double getXvalueOfFirstTemplateEntry() const {return m_xValueOfFirstTemplateEntry;}
		void setXvalueOfFirstTemplateEntry(const double xValueOfFirstTemplateEntry);
		
private:
		
		void setTF1ParNames();
    void resetTemplateRange();
    
    
    
    //Overloaded
    void setTF1Parameters();
		
		//Evaluate function that will be called by TF1 objects
    double TF1Eval(double *var, double *params);
		
		double calcSimplePedestal(const std::vector<double>& vector, const double percentage = 0.05) const;
		
		//Check if fit of template to new vector is successful
		bool fitGood(const TFitResultPtr& fitResult) const;
		
		//Fit the template to fitHist. If fit fails, retry by moving xShift to a few different values around 0.
		TFitResultPtr fitTemplate(TH1D& fitHist);
			
		//Add first vector to template (i.e. make template out of this vector)
		void addFirstVector(const std::vector<double>& newVector);
		
		//Calculate a guess for the amplitude
		double getAmplitudeGuess(const std::vector<double>& vector, const double pedestal) const;
		
		size_t getEffPeakIdx() const;
		
    //TF1 based on the template
    //This is a bit dangerous as we provide the pointer to this TF1. This means that the user can change its properties (range, parameters, etc)
    //Parameters - 0: amplitude, 1: pedestal, 2: xShift
    TF1 m_tF1;
		
		//The vector that holds the actual values
    std::vector<double> m_templateValues;
		
    const std::string m_treeName;
    const std::string m_templateValuesBranchName;
    const std::string m_dxBranchName;
    const std::string m_numAveragedFuncsBranchName;
    const std::string m_peakIdxBranchName;
  
		//Remembers how many functions were averaged to create the template
    size_t m_numAveragedFuncs;
		
    size_t m_peakIdx;

		//Difference between the points
    double m_dx;
		
    double m_minAmplitudeLimit;    
    double m_maxAmplitudeLimit;
    
    double m_minChi2_NdfLimit;
    double m_maxChi2_NdfLimit;
    
    
    double m_minXshiftLimit;
    double m_maxXshiftLimit;
    
    
    double m_minPedestalLimit;
    double m_maxPedestalLimit;
		
		//This is used to keep track where the x axis 0 is.
		//It can change if items are removed from the beggining of the template, for example.
		double m_xValueOfFirstTemplateEntry;
    
    const double m_doubleNumbersEqualThershold;
		
		int m_debugLevel;
		
		bool m_amplitudeFitEnabled;
		bool m_pedestalFitEnabled;
		bool m_xShiftFitEnabled;
		mutable bool m_useAmplitudeLimits;
		mutable bool m_useXshiftLimits;
		mutable bool m_usePedestalLimits;
	  
};

} //namespace myFuncs