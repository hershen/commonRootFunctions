#include "fftFuncs.h"

//STD
#include <iostream>

//ROOT
#include "TVirtualFFT.h"

namespace myFuncs {      
	
std::vector<TComplex> fftR2C(const std::vector<double> input, std::string options) {
	
	int inputSize = input.size(); //Should be const, but can't because TVirtualFFT::FFT needs to accept non-const pointer
	
	if(inputSize == 0)
	{
		std::cerr <<"PulseTools::RealToComplexFFT - ERROR! empty pulse passed into this function" << std::endl;
		return std::vector<TComplex>();
	}

	//Do the FFT
	std::unique_ptr<TVirtualFFT> fftr2c(TVirtualFFT::FFT(1, &inputSize, ("R2C " + options).data() ));
	fftr2c->SetPoints(input.data());
	fftr2c->Transform();
	
	//Output temperary vectors
	
	//vectors have to be right size, because we're going to overwrite their internal arrays.
	//There are inputSize/2 + 1 complex numbers in the output (fraction rounded down).
	std::vector<double> real(inputSize/2 + 1);
	std::vector<double> imag(inputSize/2 + 1);
	
	//Read transformed points
	fftr2c->GetPointsComplex(real.data(),imag.data());
	
	
	//Create output vector
	std::vector<TComplex> output;
	output.reserve( real.size() );
	
	//Take care to normalize by 1/sqrt(inputSize) because root's fft doesn't do this.
	const double one_inputSize = 1.0 / inputSize;
	
	for(unsigned int iPoint = 0; iPoint < real.size(); ++iPoint) {
		output.push_back( TComplex(real[iPoint] * one_inputSize, imag[iPoint] * one_inputSize) );
	}
	
	return output;
}
	
} //namespace histFuncs
