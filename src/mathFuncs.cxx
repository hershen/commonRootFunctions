#include "mathFuncs.h"
#include "histFuncs.h"

//std
#include <iostream>
#include <iomanip>  //for setw

//Boost
#include "boost/format.hpp"

//Root
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TH2D.h"
#include "TStyle.h"

//Mine
#include "histFuncs.h"

namespace myFuncs
{  
	
  TF1 globalFunc;
  int vecSize = 0;
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
  TF1 analytical_RC_CRn(int n, double tau, double amplitude, double startTime, double maxTime)
  {
    
    if( std::abs(maxTime) < 1e-9) 
    {
      double expCoefficient = 5.;
      maxTime = -tau*(TMath::Log( (amplitude<0?(-amplitude):amplitude) ) + expCoefficient-n*TMath::Log(tau)+TMath::Log(TMath::Factorial(n)) );
    }
    TF1 analytical("analytical","(x-[3])<0?0:[2]*pow((x-[3])/[1],[0])*TMath::Exp(-((x-[3])/[1]))/TMath::Factorial([0])",0,maxTime);
    analytical.FixParameter(0,n);
    analytical.FixParameter(1,tau);
    analytical.FixParameter(2,amplitude);
    analytical.SetParameter(3,startTime);
    
    return analytical;
  }

  
  //--------------------------------------------------------------------------------------------
  //convertArray2TF1Internal
  //********************************************************************************************
  //Internal function for convertVector2TF1.
  
  
  double convertArray2TF1Internal(double *var, double *params)
  {
    int numOfAdditionalParams = 2;
    double timeShift = params[0];
    double dT = params[1];
 
    double x = var[0] - timeShift;
    
    if(x <= 0) return params[numOfAdditionalParams];
    if(x >= dT * double(myFuncs::vecSize - 1)) return params[numOfAdditionalParams + myFuncs::vecSize - 1];
    
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
    
    myFuncs::vecSize = vecValues.size();
    
    TF1 convertArray2TF1("convertArray2TF1Internal",convertArray2TF1Internal, tMin, tMax, myFuncs::vecSize + 2);
    convertArray2TF1.SetParameter(0,timeShift);
    convertArray2TF1.FixParameter(1,dT);
    //set parameters
    for(int idx = 2; idx < myFuncs::vecSize + 2; ++idx) convertArray2TF1.FixParameter(idx, vecValues[idx-2]);
    return convertArray2TF1;
  }
    
    
  
  
  double CFDfuncInternal(double *var, double *params)
  {
    double DLY = params[0];
    double fraction = params[1];
    double x = var[0];
    return fraction * myFuncs::globalFunc.Eval(x) - myFuncs::globalFunc.Eval(x - DLY);
  }
  
  
  //--------------------------------------------------------------------------------------------
  //CFD
  //********************************************************************************************
  //Implementing a CFD for TF1
  //-------------------------------------------------------------------------------------------- 
  TF1 CFD(TF1 func, double DLY, double fraction = 0.2)
  {
    myFuncs::globalFunc = func;
    
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
    } while (xMin < xMax && std::abs(xMin-xMax) > 1e-9 );
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
  //--------------------------------------------------------------------------------------------
  template <typename Type> 
  std::vector<double> addGaussianNoise(std::vector<Type> inputValues, double mean, double  sigma, int seed)
  {
    TRandom3 random3;
    if(seed >=0) random3.SetSeed(seed);
    
    std::vector<double> noisyValues = inputValues;
    
    for(size_t idx = 0; idx < noisyValues.size(); ++idx) noisyValues[idx] += random3.Gaus(mean,sigma);
    return noisyValues;
  }
  
  //Declare all usefule implementations in order to avoid linker problems!
  template std::vector<double> addGaussianNoise<double>(std::vector<double> inputValues, double mean, double  sigma, int seed = -1);
//   template std::vector<double> addGaussianNoise<int>(std::vector<int> inputValues, double mean, double  sigma, int seed = -1);
  
  
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
    if( std::abs(xMinFracOfSigma) > 1e-9 || std::abs(xMaxFracOfSigma) > 1e-9 )
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

  
  
  //--------------------------------------------------------------------------------------------
  
  
  //********************************************************************************************
  //getCorrelationMatrix
  //********************************************************************************************
  TMatrixD getCorrelationMatrix(const TMatrixD &covarianceMatrix)
  {
    
    //Check if matrix is square
    if( covarianceMatrix.GetNcols() != covarianceMatrix.GetNrows() ) 
    {
      std::cout << "mathFuncs::getCorrelationMatrixvv : #rows != # columns. Returning original covariance matrix" << std::endl;
      return covarianceMatrix;
    }
    
    TMatrixD correlationMatrix = covarianceMatrix;
    
    for( auto i = 0; i < covarianceMatrix.GetNcols(); ++i)
    {
      if( covarianceMatrix(i,i) == 0.0 )
      {
	std::cout << "mathFuncs::getCorrelationMatrixvv : covarianceMatrix(" << i << "," << i << ") is zero. Returning original covariance matrix" << std::endl;
	return covarianceMatrix;
      }
      
      for (auto j = 0; j < covarianceMatrix.GetNcols(); ++j)
      {
	if( covarianceMatrix(j,j) == 0.0 )
	{
	  std::cout << "mathFuncs::getCorrelationMatrixvv : covarianceMatrix(" << j << "," << j << "). Returning original covariance matrix" << std::endl;
	  return covarianceMatrix;
	}
	correlationMatrix(i,j) = correlationMatrix(i,j) / TMath::Sqrt( covarianceMatrix(i,i) ) / TMath::Sqrt( covarianceMatrix(j,j) ) ;
      }
    }
    return correlationMatrix;
  }
  
  
  //********************************************************************************************
  //printMatrix
  //********************************************************************************************
  void printMatrix(const TMatrixD matrix, const std::vector<std::string>& headings, int width, std::ios_base& align( std::ios_base& str ), const std::string& precision)
  {
    std::string formatArg = "%1$" + precision;
    
    std::cout << std::setw(width) << align << "";
    
    //Print upper row headers
    for (size_t idx = 0; idx < headings.size(); ++idx) std::cout << std::setw(width) << align <<  headings[idx];
    std::cout << std::endl;
    
    //Print values
		//Prevent warning of unsigned integer expression because TMatrixD::GetNrows returns int isntead of unsigned int
		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wsign-compare" 
    for (size_t i = 0; i < matrix.GetNrows(); ++i)
		#pragma GCC diagnostic pop
    {
      //Print heading at begining of line
      if( i < headings.size() ) std::cout << std::setw(width) << align <<  headings[i];
      
      for (auto j = 0; j < matrix.GetNcols(); ++j) std::cout << std::setw(width) << align << str( boost::format(formatArg) % matrix(i,j) );
      
      std::cout << std::endl;
    }
  }
  
  //********************************************************************************************
  //drawMatrix
  //********************************************************************************************
  TCanvas* drawMatrix(const TMatrixD matrix, std::string title, const std::vector<std::string>& xAxisHeadings, const std::vector<std::string>& yAxisHeadings, const double zMin, const double zMax, const std::string& precision )
  {
		//Prevent warning of unsigned integer expression because TMatrixD::GetNcols returns int isntead of unsigned int
		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wsign-compare" 
    if( xAxisHeadings.size() < matrix.GetNcols() )
		#pragma GCC diagnostic pop
    {
      std::cout << "mathFuncs::drawMatrix : matrix has " << matrix.GetNcols() << "columns, but only "<< xAxisHeadings.size() << " headings. Abroting. " << std::endl;
      return new TCanvas();
    }
    
    //Prevent warning of unsigned integer expression because TMatrixD::GetNcols returns int isntead of unsigned int
		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wsign-compare" 
    if( yAxisHeadings.size() < matrix.GetNcols() )
		#pragma GCC diagnostic pop
    {
      std::cout << "mathFuncs::drawMatrix : matrix has " << matrix.GetNrows() << "rows, but only "<< yAxisHeadings.size() << " headings. Abroting. " << std::endl;
      return new TCanvas();
    }
      
    TH2D* hist = new TH2D(title.data(), title.data(), matrix.GetNcols(), 0, matrix.GetNcols(), matrix.GetNrows(), 0, matrix.GetNrows());
    
    for( auto i = 0; i < matrix.GetNcols(); ++i)
      for(auto j = matrix.GetNrows() - 1; j >= 0 ; --j)
	hist->Fill( xAxisHeadings[i].data(), yAxisHeadings[j].data(), matrix(i,j) );

    //Set minimum and maximum z axis (for plotting)
    if ( std::abs(zMin) > 1e-6 || std::abs(zMax) > 1e-6)
    {
      hist->SetMinimum(zMin);
      hist->SetMaximum(zMax);
    }
    
    //Set precision format for text labels
    const char* originalPaintTextFormat = gStyle->GetPaintTextFormat();
    gStyle->SetPaintTextFormat( precision.data() );
    
    
    hist->LabelsOption( "d" ); //Draw x axis labels pointing down
    hist->SetStats(false);     //No stat box
    
    TCanvas* canvas = myFuncs::newCanvas( title + "Canvas" );
    canvas->SetGrid();
    canvas->SetTicks();
    
    gStyle->SetPalette(67); //67 - kAurora, 84 - kGreenPink
    
    hist->Draw("colz"); // color pads   
    canvas->Update();

    hist->Draw("textsame");  // add text
    canvas->Update();
    
    gStyle->SetPaintTextFormat( originalPaintTextFormat );
    
    return canvas;
  }
    
} //namespace myMathFunctions
