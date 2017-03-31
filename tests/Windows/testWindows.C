#include "Windows.h"
#include "TMath.h"
#include <iostream>

bool testHamming() 
{
	bool pass = true;
	
	const myFuncs::Hamming hamming(10);
	
	if( std::abs(hamming.eval(0)-0.08) > 1e-9) {
		std::cout << "testHamming: expecting hamming(0) = 0.0, obtained " << hamming.eval(0) << std::endl;
		pass = false;
	}
	
	if( std::abs(hamming.eval(10)-0.08) > 1e-9) {
		std::cout << "testHamming: expecting hamming(10) = 0.0, obtained " << hamming.eval(10) << std::endl;
		pass = false;
	}
	
	std::vector<double> v(11,1.0);
	
	const auto windowed = hamming.windowAnInput(v);
	
	for(unsigned int i = 0; i < windowed.size(); ++i) {
		if( std::abs(windowed[i] - (0.54-0.46*TMath::Cos(2.0*TMath::Pi()*i/10.0)) ) > 1e-9 ) {
			std::cout << "testHamming: expecting windowed[" << i << "] = " << 0.54-0.46*TMath::Cos(2.0*TMath::Pi()*i/10.0) << ", obtained " << windowed[i] << std::endl;
			pass = false;
		}
	}
	
	return pass;
}

bool testHann() 
{
	bool pass = true;
	
	const myFuncs::Hann hann(10);
	
	if( std::abs(hann.eval(0)-0.0) > 1e-9) {
		std::cout << "testHann: expecting hann(0) = 0.0, obtained " << hann.eval(0) << std::endl;
		pass = false;
	}
	
	if( std::abs(hann.eval(10)-0.0) > 1e-9) {
		std::cout << "testHann: expecting hann(10) = 0.0, obtained " << hann.eval(10) << std::endl;
		pass = false;
	}
	
	std::vector<double> v(11,1.0);
	
	const auto windowed = hann.windowAnInput(v);
	
	for(unsigned int i = 0; i < windowed.size(); ++i) {
		if( std::abs(windowed[i] - (0.5-0.5*TMath::Cos(2.0*TMath::Pi()*static_cast<double>(i)/10.0) ) ) > 1e-9 ) {
			std::cout << "testHann: expecting windowed[" << i << "] = " << 0.5-0.5*TMath::Cos(2.0*TMath::Pi()*static_cast<double>(i)/10.0) << ", obtained " << windowed[i] << std::endl;
			pass = false;
		}
	}
	
	return pass;
}

void testWindows()
{
	auto result = testHann();
	std::cout << "testHann test result = " << result << std::endl;
	
	result = testHamming();
	std::cout << "testHamming test result = " << result << std::endl;
}

int main()
{
	testWindows();
}