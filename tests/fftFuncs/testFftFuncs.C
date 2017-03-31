#include "fftFuncs.h"
#include "mathFuncs.h"


bool testRealSequence2psd(const size_t N) {
	
	std::vector<double> x;
	x.reserve(N);
	
	TF1 cosFunc("cosFunc","1 + TMath::Cos(2.0*TMath::Pi()/2.0 * x) +  TMath::Cos(2.0*TMath::Pi()/5.0 * x) + TMath::Cos(2.0*TMath::Pi()/10.0 * x)",0, N);
		
	//Create cosine samples
	for (auto i = 0; i < N; ++i) {
		x.push_back( cosFunc.Eval(i) );
	}
	
	const auto psd = myFuncs::realSequence2psd(x);
	
	if(psd.size() != (N/2 + 1)) {
		std::cerr << "testRealSequence2psd: Failed - size of psd is " << psd.size() << ", expected " << N/2 + 1 << std::endl;
		return false;
	}
	
	const double xSquared = myFuncs::sumVectorSquared(x);
	const double psdSum = myFuncs::sumVector(psd);
	
	if( std::abs(xSquared - psdSum) > 1e-9) {
		std::cerr << "testRealSequence2psd: Perseval's theorem failed. sum(x[i]*x[i]) =  " << xSquared << ", sum(psd[i]) (already squared) =  " << psdSum << std::endl;
		return false;
	}
	
	
	//Check individual entries
	if(N==100) {
		for(size_t i = 0; i < psd.size(); ++i) {
			
			if(i == 0) {
					if(std::abs(psd[i] - 100) > 1e-9) {
						std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 100 " << std::endl;
						return false;
				}
			}
			else if(i == 10) {
					if( std::abs(psd[i] - 50) > 1e-9) {
						std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 50 " << std::endl;
						return false;
				}
			}
			else if(i == 20) {
					if(std::abs(psd[i] - 50) > 1e-9) {
						std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 50 " << std::endl;
						return false;
				}
			}
			else if(i == 50) {
					if(std::abs(psd[i] - 100) > 1e-9) {
						std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 100 " << std::endl;
						return false;
				}
			}
			else if(std::abs(psd[i]) > 1e-9) {
				std::cerr << "testRealSequence2psd: psd[" << i << "] =  " << psd[i] << ", expecting 0 " << std::endl;
				return false;
			}
			
		}
	}
	
	return true;
		
}

void testFftFuncs() {
	
	//The 100 is important because we check individual results of the psd
	if (!testRealSequence2psd(100)) {
		std::cerr << "testRealSequence2psd with input 100 Failed!!!" << std::endl;
	}
	
	if (!testRealSequence2psd(101)) {
		std::cerr << "testRealSequence2psd with input 101 Failed!!!" << std::endl;
	}
}