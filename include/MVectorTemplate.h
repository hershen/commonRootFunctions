#pragma once

#ifndef MVectorTemplate_cxx
#define MVectorTemplate_cxx


#include <vector>
#include "TF1.h"

class MVectorTemplate
{

  public:
    MVectorTemplate();
    MVectorTemplate(const std::vector<double> &newVec, const double &newDx);
    
    
    inline unsigned int getNumAveragedFuncs() const {return rNumAveragedFuncs;}
    
    inline double getDx() const {return rDx;}
    void setDx(const double &newDx);
    
    //inline void setVectorValues(const std::vector<double> &newValues) {rTemplateValues = newValues; rNumAveragedFuncs = 0;}
    inline std::vector<double> getTemplateValues() const {return rTemplateValues;}
    
    void resetTemplateValues(); 
    
    // ------------------------------------------------------------------
    //enableAmplitudeFit
    //Currently no way of actually setting amplitude of fit. It is guessed when a vector is added.
    // ------------------------------------------------------------------
    void enableAmplitudeFit(const bool &enableAmplitudeFit);
    inline bool isAmplitudeFitEnabled() const {return bAmplitudeFitEnabled;}
    
    // ------------------------------------------------------------------
    //enablePedestalFit
    //Currently no way of actually setting pedestal of fit. It is guessed when a vector is added.
    // ------------------------------------------------------------------
    void enablePedestalFit(const bool &enablePedestalFit);
    inline bool isPedestalFitEnabled() const {return bPedestalFitEnabled;}
    
    inline unsigned int getTemplateSize() const {return rTemplateValues.size();}
    
    
    //Adds a new vector to the template. The values are averaged (taking into account the relative wight of the new vector according to the number of previous vectors added.
    //rTF1 parameters reset at the end (otherwise they keep values of last vector added.
    int addVector(const std::vector<double> &newVector, const double &std);
    
    inline TF1 *getTF1() {return &rTF1;}
    
    inline int getDebugLevel() const {return rDebugLevel;}
    inline void setDebugLevel(const int &newDebugLevel) {rDebugLevel = newDebugLevel;}
    
    inline void getAmplitudeLimits(double &minAmplitude, double &maxAmplitude) const {minAmplitude = rMinAmplitudeLimit; maxAmplitude = rMaxAmplitudeLimit;}
    void setAmplitudeLimits(const double &newMinAmplitudeLimit, const double &newMaxAmplitudeLimit);
    
    inline void getChi2Limits(double &minChi2, double &maxChi2) const {minChi2 = rMinChi2Limit; maxChi2 = rMaxChi2Limit;}
    void setChi2Limits(const double &newMinChi2Limit, const double &newMaxChi2Limit) { rMinChi2Limit = newMinChi2Limit; rMaxChi2Limit = newMaxChi2Limit;}
    
    inline void getXshiftLimits(double &minXshift, double &maxXshift) const {minXshift = rMinXshiftLimit; maxXshift = rMaxXshiftLimit;}
    void setXshiftLimits(const double &newMinXshiftLimit, const double &newMaxXshiftLimit);
    
    inline void getPedestalLimits(double &minPedestal, double &maxPedestal) const {minPedestal = rMinPedestalLimit; maxPedestal = rMaxPedestalLimit;}
    void setPedestalLimits(const double &newMinPedestalLimit, const double &newMaxPedestalLimit);
    
    inline size_t getPeakIdx() const {return rPeakIdx;}
    
    int saveTemplateToTFile(const std::string &fullFileName, const std::string &treeDescription);
       
    int loadTemplateFromTFile(const std::string &fullFileName);
    
  private:
    
    //Remembers how many functions were averaged to create the template
    size_t rNumAveragedFuncs;
    
    //Difference between the points
    double rDx;
    
    //The vector that holds the actual values
    std::vector<double> rTemplateValues;
    
    bool bAmplitudeFitEnabled;
    
    bool bPedestalFitEnabled;
    
    //TF1 based on the template
    //This is a bit dangerous as we provide the pointer to this TF1. This means that the user can change its properties (range, parameters, etc)
    //Parameters - 0: amplitude, 1: pedestal, 2: xShift
    TF1 rTF1;
    
    size_t rPeakIdx;
    
    int rDebugLevel;
    
    mutable bool rUseAmplitudeLimits;
    double rMinAmplitudeLimit;    
    double rMaxAmplitudeLimit;
    
    double rMinChi2Limit;
    double rMaxChi2Limit;
    
    mutable bool rUseXshiftLimits;
    double rMinXshiftLimit;
    double rMaxXshiftLimit;
    
    mutable bool rUsePedestalLimits;
    double rMinPedestalLimit;
    double rMaxPedestalLimit;
    
    //Evaluate function that will be called by TF1 objects
    double TF1Eval(double *var, double *params);
    
    const std::string rTreeName;
    const std::string rTemplateValuesBranchName;
    const std::string rDxBranchName;
    const std::string rNumAveragedFuncsBranchName;
    const std::string rPeakIdxBranchName;
    
    const double rDoubleNumbersEqualThershold;
 
    void setTF1ParNames();
    void resetTemplateRange();
    
    void setTF1Parameters(const double &amplitude, const double &pedestal, const double &xShift);
    
    //Overloaded
    void setTF1Parameters();
    
};
 

#endif