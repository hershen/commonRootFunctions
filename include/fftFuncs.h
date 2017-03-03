#pragma once

#include <vector>
#include <string>

//ROOT
#include "TComplex.h"

namespace myFuncs {

	//--------------------------------------------
	//Insperation taken from CDMS's PulseTools::RealToComplexFFT
	//--------------------------------------------
	std::vector<TComplex> fftR2C(const std::vector<double> input, std::string options = "P");

}//myFuncs namespace