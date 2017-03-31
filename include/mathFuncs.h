#pragma once

//std
#include <iostream>
#include <numeric> //For std::accumulate

//Root
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TMatrixD.h"
#include "TCanvas.h"

//For TMiuit2
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"

namespace myFuncs
{  
	
  //trying to compile with g++
//   int vecSize = 0;
//   TF1 globalFunc;
  
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
  TF1 analytical_RC_CRn(int n, double tau = 500., double amplitude = 1., double startTime = 0., double maxTime = 0.);
    
  //--------------------------------------------------------------------------------------------
  //calcResiduals
  //********************************************************************************************
  //Function to calculate residuals between a a model function (TF1) and y(x) yValues taken at points xValuesV.
	//Residual defined as f(x) - yValue
  //x values are in vector xValues of type xValueType
  //y values are in vector yValues of type yValueType
	//If there are more yValues than xValues, the residuals are calculated only for the xValues given.
	//More xValues than yValues are not allowed.
  //--------------------------------------------------------------------------------------------
	//Not implemented with boost zip itirators because ROOT doesn't compile Boost libraries well.
  template <typename xValType, typename yValType> 
  std::vector<double> calcResiduals(const std::vector<xValType>& xValues, const std::vector<yValType>& yValues, const std::vector<double>& stds, const TF1& modelFunc)
	{
    // ------------------------------------------------------------------
    //Sanity checks
    // ------------------------------------------------------------------
		//Allow cases where there are more y values than x values
    if(yValues.size() < xValues.size()) 
      throw std::invalid_argument("yValuesV.size() < xValuesV.size()");
		
		if(stds.size() < xValues.size()) 
      throw std::invalid_argument("stds.size() < xValuesV.size()");
		
		//Output vector
    std::vector<double> residuals;
		residuals.reserve( xValues.size() );
		
		auto xIt = xValues.begin();
		auto yIt = yValues.begin();
		auto stdIt = stds.begin();
		//Loop on num elements = xValues.size()
		while( xIt != xValues.end() )
		{
// 			std::cout << "xVal = " << *xIt << ", yVal = " << *yIt << std::endl;
			residuals.push_back ( (modelFunc.Eval(*xIt) - static_cast<double>(*yIt)) / *stdIt );
		  ++xIt;
			++yIt;
			++stdIt;
		}
		
		return residuals;
	}
	
	//Overloaded
	template <typename xValType, typename yValType> 
	std::vector<double> calcResiduals(const std::vector<xValType>& xValues, const std::vector<yValType>& yValues, const double std, const TF1& modelFunc)
	{
		return calcResiduals(xValues, yValues, std::vector<double>(xValues.size(), std), modelFunc);
	}
	
	
  //--------------------------------------------------------------------------------------------
  //convertArray2TF1Internal
  //********************************************************************************************
  //Internal function for convertVector2TF1.
  double convertArray2TF1Internal(double *var, double *params);
  
  
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
  TF1 convertVector2TF1(double dT, std::vector<double> vecValues, double timeShift);    
    
  
  
  double CFDfuncInternal(double *var, double *params);
   
  //--------------------------------------------------------------------------------------------
  //CFD
  //********************************************************************************************
  //Implementing a CFD for TF1
  //-------------------------------------------------------------------------------------------- 
  TF1 CFD(TF1 func, double DLY, double fraction);

  
  //Assumes that it is a positive signal
  double CFtime(TF1 func, double DLY, double fraction/*, double initialGuess = 0.*/);
  
  

  //--------------------------------------------------------------------------------------------
  //cfdCR_RCnZeroCrossingTime
  //********************************************************************************************
  //Returns the zero crossing time of a CFD function from an analytical CR-(RC)^n output of a step function starting at t=0
  
  double cfdCR_RCnZeroCrossingTime(int n, double fraction, double tau, double DLY);
  

  
  //--------------------------------------------------------------------------------------------
  //addGaussianNoise
  //********************************************************************************************
  //Add gaussian noise to an input vector of type Type. Output vector is alway double.
  //If seed is not given the default constructor is used. NOTE - the number sequence generated is always the same in this case.
  //If a seed is provided, it is used.
  //--------------------------------------------------------------------------------------------
  template <typename Type> 
  std::vector<double> addGaussianNoise(std::vector<Type> inputValues, double mean, double  sigma, int seed = -1);
  
  
  //--------------------------------------------------------------------------------------------
  //getGaussianFit
  //********************************************************************************************
  //Fit a gaussian to hist.
  //If xMinInitial or xMaxInitial are given, the fit is performed only in the range (xMinInitial, xMaxInitial)
  //If xMinFracOfSigma or xMaxFracOfSigma are given, the fit is re-performed in the range (mean - sigma*xMinFracOfSigma, mean + sigma*xMaxFracOfSigma) where mean and sigma are taken from the first fit. I.e., the fit can be run again in the range (mean - 2std, mean + 2std).
  //The function returns the fit result in fitResult
  //--------------------------------------------------------------------------------------------
  TF1 getGaussianFit(TH1D hist, TFitResultPtr &fitResult, double xMinInitial, double xMaxInitial, double xMinFracOfSigma, double xMaxFracOfSigma);
  
  
  //--------------------------------------------------------------------------------------------
  //getGaussianFitResult - Overloaded
  //********************************************************************************************
  //Same as before, only without the fitResult parameter.
  //--------------------------------------------------------------------------------------------
  inline TF1 getGaussianFit(TH1D hist, double xMinInitial, double xMaxInitial, double xMinFracOfSigma, double xMaxFracOfSigma);
  
  //--------------------------------------------------------------------------------------------
  //getCorrelationMatrix
  //********************************************************************************************
  //Calculate the correlation matrix from the covariance matrix.
  //I.e. M_ij = M_ij / sqrt(M_ii) / sqrt(M_jj)
  //If one of the diagonal elements M_ii is zero, return original matrix
  //--------------------------------------------------------------------------------------------
  TMatrixD getCorrelationMatrix(const TMatrixD &covarianceMatrix);
  
  //--------------------------------------------------------------------------------------------
  //printMatrix
  //********************************************************************************************
  //Print a matrix on screen
  //headings - the heading of each row / column
  //width - each element will be padded to fill width spaces
  //aligh - align to left, right, etc.
  //precision - the precision syntax to be used (as in printf )
  //--------------------------------------------------------------------------------------------
  void printMatrix(const TMatrixD matrix, const std::vector<std::string>& headings, int width = 10, std::ios_base& align( std::ios_base& str ) = std::left, const std::string& precision = ".3g" );
  
  //--------------------------------------------------------------------------------------------
  //drawMatrix
  //********************************************************************************************
  //Print a matrix on screen
  //headings - the heading of each row / column
  //precision - the precision syntax to be used (as in printf )
  //zMin, zMax - used to set z axis limits
  //--------------------------------------------------------------------------------------------
  TCanvas* drawMatrix(const TMatrixD matrix, std::string title, const std::vector<std::string>& xAxisHeadings, const std::vector<std::string>& yAxisHeadings, const double zMin = 0.0, const double zMax = 0.0, const std::string& precision = ".3g" );
	
	
	//--------------------------------------------------------------------------------------------
  //sumVector
  //********************************************************************************************
	//Returns sum(x[i])
  //wrapper around std::accumulate with inital value of 0.0 so not to confuse it with 0 in which case everything is converted to an int.
  //--------------------------------------------------------------------------------------------
	template <typename T>
	T sumVector(const std::vector<T>& vector, const T& initialValue = 0.0) {
		return std::accumulate(vector.begin(), vector.end(), initialValue);
	}
	
	//--------------------------------------------------------------------------------------------
  //sumVectorSquared
  //********************************************************************************************
  //Returns sum(x[i]*x[i])
  //--------------------------------------------------------------------------------------------
	template <typename T>
	T sumVectorSquared(const std::vector<T>& vector, const T& initialValue = 0.0) {
		auto sum = initialValue;
		
		for(auto& element : vector)
			sum += element*element;
		
		return sum;
	}

} //namespace myMathFunctions

