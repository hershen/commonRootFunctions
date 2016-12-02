#pragma once

#ifndef myMathFunctions_cxx
#define myMathFunctions_cxx

#include <iostream>
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TFitResult.h"

//For TMiuit2
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"

TRandom3 random3;
int randomSeed = 1234;

namespace mathFuncs
{
  
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
  TF1 analytical_RC_CRn(int n, double tau = 500., double amplitude = 1., double startTime = 0., double maxTime = 0.)
  {
    double expCoefficient = 5.;
    if( abs(maxTime) < 1e-9) maxTime = -tau*(TMath::Log( (amplitude<0?(-amplitude):amplitude) ) + expCoefficient-n*TMath::Log(tau)+TMath::Log(TMath::Factorial(n)) );
    TF1 analytical("analytical","(x-[3])<0?0:[2]*pow((x-[3])/[1],[0])*TMath::Exp(-((x-[3])/[1]))/TMath::Factorial([0])",0,maxTime);
    analytical.FixParameter(0,n);
    analytical.FixParameter(1,tau);
    analytical.FixParameter(2,amplitude);
    analytical.SetParameter(3,startTime);
    
    return analytical;
  }
  
  //--------------------------------------------------------------------------------------------
  //calcResiduals
  //********************************************************************************************
  //Function to calculate residuals between a a model function (TF1) and y(x) yValuesV taken at points xValuesV.
  //x values are in vector xValuesV of type xValueType
  //y values are in vector yValuesV of type yValueType
  //--------------------------------------------------------------------------------------------
  template <typename xValType, typename yValType> 
  std::vector<double> calcResiduals(std::vector<xValType> xValuesV, std::vector<yValType> yValuesV, TF1 modelFunc)
  {
    // ------------------------------------------------------------------
    //Sanity checks
    // ------------------------------------------------------------------
    if(xValuesV.size() != yValuesV.size()) 
    {
      std::cout << "calcResiduals: xValuesV.size() != yValuesV.size(). Aborting." << std::endl;
      return std::vector<double>(0,0);
    }
    
    //return vector
    std::vector<double> residualV = std::vector<double>(xValuesV.size(),0);
    
    
    // ------------------------------------------------------------------
    //Calculate residuals
    // ------------------------------------------------------------------
    for(size_t idx = 0; idx < xValuesV.size(); ++idx) residualV[idx] = double(yValuesV[idx]) - modelFunc.Eval(xValuesV[idx]);
    
    
    return residualV;
  }
 
  //--------------------------------------------------------------------------------------------
  //vecMean
  //********************************************************************************************
  //Calculate the mean of a vector
  //y values are in vector yValuesV of type yValueType
  //-------------------------------------------------------------------------------------------- 
  template <typename Type> 
  double vecMean(std::vector<Type> vec)
  {
    Type sum = 0;
    for(size_t idx = 0; idx < vec.size(); ++idx) sum += vec[idx];
    
    return double(sum)/double(vec.size());
  }
  
  //--------------------------------------------------------------------------------------------
  //vecSum2
  //********************************************************************************************
  //Calculate Sum( x_i ^2 )
  //-------------------------------------------------------------------------------------------- 
  template <typename Type> 
  Type vecSum2(std::vector<Type> vec)
  {
    Type sum2 = 0;
    for(size_t idx = 0; idx < vec.size(); ++idx) sum2 += vec[idx]*vec[idx];
    
    return sum2;
  }

  
    
  //--------------------------------------------------------------------------------------------
  //convertArray2TF1Internal
  //********************************************************************************************
  //Internal function for convertVector2TF1.
  
  int vecSize = 0;
  double convertArray2TF1Internal(double *var, double *params)
  {
    int numOfAdditionalParams = 2;
    double timeShift = params[0];
    double dT = params[1];
 
    double x = var[0] - timeShift;
    
    if(x <= 0) return params[numOfAdditionalParams];
    if(x >= dT * double(vecSize - 1)) return params[numOfAdditionalParams + vecSize - 1];
    
    int idx = int(x / dT); // x corrisponds to a time between idx and idx+1
    
    double x1 = double(idx)*dT, x2 = double(idx + 1) *dT;
    double y1 = params[numOfAdditionalParams + idx], y2 = params[numOfAdditionalParams + idx + 1];
    
//     std::cout << "x1 = " << x1 << ", y1 = " << ". x2 = " << x2 << ", y2 = " << y2 << std::endl;
    
    return (y2 - y1)/(x2-x1)*(x-x1) + y1;    
  }

  
  //--------------------------------------------------------------------------------------------
  //convertVector2TF1
  //********************************************************************************************
  //Converts a vector of double values into a TF1.
  //Only parameter 0 should regularly change! Others change the function itself.
  
  
  //It is assumed that the values in vecValues are evenly spaced with spacing dT. I.e. vecValues[0] = f(0), vecValues[1] = f(dT), vecValues[2] = f(2dT)...
  //timeShift can shift the x axis.
  //Values between the ones in vecValues are evaluated using a linear interpolation between the points.
  //Evals for x axis values smaller than the first entry in vecValues return vecValues[0]. Evals for values greater than the last entry in vecValues return vecValues[last].
  //The params array given to convertArray2TF1Internal is made up of {timeShift, dT, vecValues}
  //All parameters are fixed, except timeshift which floats.
  //There is a file scope "global" variable vecSize which holds the number of parameters in vecValues (because the internal function doesn't know this.
  //The 
  //All values of vecValues are given to the internal function
  //-------------------------------------------------------------------------------------------- 
  TF1 convertVector2TF1(double dT, std::vector<double> vecValues, double timeShift = 0.)
  {
    double tMin = timeShift;
    double tMax = (vecValues.size()-1) * dT + timeShift;
    
    vecSize = vecValues.size();
    
    TF1 convertArray2TF1("convertArray2TF1Internal",convertArray2TF1Internal, tMin, tMax, vecSize + 2);
    convertArray2TF1.SetParameter(0,timeShift);
    convertArray2TF1.FixParameter(1,dT);
    //set parameters
    for(int idx = 2; idx < vecSize + 2; ++idx) convertArray2TF1.FixParameter(idx, vecValues[idx-2]);
    return convertArray2TF1;
  }
    
    
  TF1 globalFunc;
  
  double CFDfuncInternal(double *var, double *params)
  {
    double DLY = params[0];
    double fraction = params[1];
    double x = var[0];
    return fraction * globalFunc.Eval(x) - globalFunc.Eval(x - DLY);
  }
  
  
  //--------------------------------------------------------------------------------------------
  //CFD
  //********************************************************************************************
  //Implementing a CFD for TF1
  //-------------------------------------------------------------------------------------------- 
  TF1 CFD(TF1 func, double DLY, double fraction = 0.2)
  {
    globalFunc = func;
    
    TF1 CFDfunc("CFDfunc",CFDfuncInternal, func.GetXmin(), func.GetXmax(), 2);
    CFDfunc.SetParameters(DLY,fraction);
    return CFDfunc;
  }
  
  //Assumes that it is a positive signal
  double CFtime(TF1 func, double DLY, double fraction = 0.2/*, double initialGuess = 0.*/)
  {
    double minMaxPrecision = 1e-3;
    TF1 CFDfunc = CFD(func, DLY, fraction );
    size_t originalNpx = CFDfunc.GetNpx();

    CFDfunc.SetNpx(100);

    // ------------------------------------------------------------------
    //find minimum and maximum of function
    //Incorporates a machinism to make sure that xMax < xMin. This is only correct for positive signals. If xMax > xMin, start looking again from xMin.
    // ------------------------------------------------------------------
    double xStart = CFDfunc.GetXmin();
    double xMin, xMax;
    do
    {
      xMin = CFDfunc.GetMinimumX(xStart,CFDfunc.GetXmax(),minMaxPrecision);
      xMax = CFDfunc.GetMaximumX(xStart,CFDfunc.GetXmax(),minMaxPrecision);
      xStart = xMax;
    } while (xMin < xMax && abs(xMin-xMax) > 1e-9 );
//     std::cout << "xMin = " << xMin << ", xMax = " << xMax << std::endl;
    
    
    // ------------------------------------------------------------------
    //find zero crossing
    // ------------------------------------------------------------------
    //GetX(Double_t fy,Double_t xmin = 0,Double_t	xmax = 0, Double_t epsilon = 1.E-10, Int_t maxiter = 100, Bool_t logx = false)
    double zeroTime = CFDfunc.GetX(0.,xMax,xMin,1e-20,100,false );
    CFDfunc.SetNpx(originalNpx);
    return zeroTime;
    
  }
  

  //--------------------------------------------------------------------------------------------
  //cfdCR_RCnZeroCrossingTime
  //********************************************************************************************
  //Returns the zero crossing time of a CFD function from an analytical CR-(RC)^n output of a step function starting at t=0
  
  double cfdCR_RCnZeroCrossingTime(int n, double fraction, double tau, double DLY)
  {
    double exp = TMath::Exp(DLY/tau/double(n));
    double zeroTime = -DLY*exp / (TMath::Power(fraction,1./double(n)) - exp) ; //TMath::Power(TMath::Exp(DLY/tau)/fraction, 1./double(n) );
    return zeroTime;
  }

  
  //--------------------------------------------------------------------------------------------
  //addGaussianNoise
  //********************************************************************************************
  //Add gaussian noise to an input vector of type Type. Output vector is alway double.
  //If seed is not given than the default filescope seed is used.
  //--------------------------------------------------------------------------------------------
  template <typename Type> 
  std::vector<double> addGaussianNoise(std::vector<Type> inputValues, double mean, double  sigma, int seed = 0)
  {
    if(seed == 0) random3.SetSeed(randomSeed);
    else random3.SetSeed(seed);
    
    std::vector<double> noisyValues = inputValues;
    
    for(size_t idx = 0; idx < noisyValues.size(); ++idx) noisyValues[idx] += random3.Gaus(mean,sigma);
    return noisyValues;
  }
  
  //--------------------------------------------------------------------------------------------
  //getGaussianFit
  //********************************************************************************************
  //Fit a gaussian to hist.
  //If xMinInitial or xMaxInitial are given, the fit is performed only in the range (xMinInitial, xMaxInitial)
  //If xMinFracOfSigma or xMaxFracOfSigma are given, the fit is re-performed in the range (mean - sigma*xMinFracOfSigma, mean + sigma*xMaxFracOfSigma) where mean and sigma are taken from the first fit. I.e., the fit can be run again in the range (mean - 2std, mean + 2std).
  //The function returns the fit result in fitResult
  //--------------------------------------------------------------------------------------------
  TF1 getGaussianFit(TH1D hist, TFitResultPtr &fitResult, double xMinInitial = 0., double xMaxInitial = 0., double xMinFracOfSigma = 0., double xMaxFracOfSigma = 0. )
  {

    TF1 gausFunc("gausFunc", "gaus", xMinInitial,  xMaxInitial);//hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
    
    fitResult = hist.Fit("gausFunc","S", "", xMinInitial, xMaxInitial);
    
    if(fitResult) return TF1(); //fit not succesful

    // ------------------------------------------------------------------
    //Re-run fit in the range (mean - sigma*xMinFracOfSigma, mean + sigma*xMaxFracOfSigma)
    // ------------------------------------------------------------------
    if( abs(xMinFracOfSigma) > 1e-9 || abs(xMaxFracOfSigma) > 1e-9 )
    {
      double mean = fitResult->Parameter(1);
      double std = fitResult->Parameter(2);
      
      double xMin = mean - xMinFracOfSigma*std;
      double xMax = mean + xMaxFracOfSigma*std;
      fitResult = hist.Fit("gausFunc","S","", xMin, xMax);
    
      //Change function range
      gausFunc.SetRange(xMin, xMax);
    }
    
    return gausFunc;

  }
  
  //--------------------------------------------------------------------------------------------
  //getGaussianFitResult - Overloaded
  //********************************************************************************************
  //Same as before, only without the fitResult parameter.
  //--------------------------------------------------------------------------------------------
  inline TF1 getGaussianFit(TH1D hist, double xMinInitial = 0., double xMaxInitial = 0., double xMinFracOfSigma = 0., double xMaxFracOfSigma = 0.)
  {
    TFitResultPtr dummyFitResult;
    return getGaussianFit(hist, dummyFitResult, xMinInitial, xMaxInitial, xMinFracOfSigma, xMaxFracOfSigma);
  }

} //namespace myMathFunctions

#endif